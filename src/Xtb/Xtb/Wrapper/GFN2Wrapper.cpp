/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Internal Includes */
#include "Xtb/Wrapper/GFN2Wrapper.h"
#include "Xtb/Wrapper/XtbSettings.h"

/* External Include */
#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/CalculatorBasics/ResultsAutoCompleter.h>
#include <Utils/GeometricDerivatives/NumericalHessianCalculator.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <Utils/Solvation/ImplicitSolvation.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <xtb.h>
#include <boost/exception/diagnostic_information.hpp>

namespace Scine {
namespace Xtb {

std::mutex GFN2Wrapper::_mtx;

GFN2Wrapper::GFN2Wrapper() {
  _settings.modifyString(Utils::SettingsNames::method, this->method());
}

const Scine::Utils::Results& GFN2Wrapper::calculate(std::string /* dummy */) {
  if (!_settings.valid()) {
    _settings.throwIncorrectSettings();
  }
  verifyPesValidity();
#if defined(_OPENMP)
  const int nCores = _settings.getInt(Utils::SettingsNames::externalProgramNProcs);
  omp_set_dynamic(0); // Explicitly disable dynamic teams
  omp_set_num_threads(nCores);
#endif
  // Prepare Data
  const int natoms = _structure->size();
  auto elements = _structure->getElements();
  Eigen::VectorXi attyp(natoms);
  for (int i = 0; i < natoms; i++) {
    attyp[i] = Scine::Utils::ElementInfo::Z(elements[i]);
  }
  const double charge = _settings.getInt(Utils::SettingsNames::molecularCharge); // double because xtb wants double
  const int uhf = _settings.getInt(Utils::SettingsNames::spinMultiplicity) - 1;
  auto coord = _structure->getPositions();

  // Prepare XTB classes
  xtb_TEnvironment env = xtb_newEnvironment();
  xtb_TCalculator calc = xtb_newCalculator();
  xtb_TResults res = xtb_newResults();
  xtb_TMolecule mol = xtb_newMolecule(env, &natoms, attyp.data(), coord.data(), &charge, &uhf, nullptr, nullptr);
  if (xtb_checkEnvironment(env) != 0) {
    xtb_showEnvironment(env, nullptr);
    _cleanDataStructures(env, calc, res, mol);
    throw std::runtime_error("XTB molecule setup failed.");
  }

  // Setup XTB model
  _mtx.lock();
  xtb_loadGFN2xTB(env, mol, calc, nullptr);
  _mtx.unlock();
  if (xtb_checkEnvironment(env) != 0) {
    xtb_showEnvironment(env, nullptr);
    _cleanDataStructures(env, calc, res, mol);
    throw std::runtime_error("XTB method setup failed.");
  }

  // Apply settings
  double acc = _settings.getDouble(Utils::SettingsNames::selfConsistenceCriterion) / 1e-6; // to arrive at Xtb accuracy
                                                                                           // value
  xtb_setAccuracy(env, calc, acc);
  xtb_setMaxIter(env, calc, _settings.getInt(Utils::SettingsNames::maxScfIterations));
  xtb_setElectronicTemp(env, calc, _settings.getDouble(Utils::SettingsNames::electronicTemperature));
  xtb_setVerbosity(env, _settings.getInt("print_level"));

  if (Utils::Solvation::ImplicitSolvation::solvationNeededAndPossible(_availableSolvationModels, _settings)) {
    std::string solvent = _settings.getString(Utils::SettingsNames::solvent);
    std::for_each(solvent.begin(), solvent.end(), [](char& c) { c = ::tolower(c); });
    std::vector<std::string> availableSolvents = {"acetone", "acetonitrile", "benzene", "ch2cl2",   "chcl3",
                                                  "cs2",     "dmso",         "ether",   "methanol", "toluene",
                                                  "thf",     "water",        "h2o"};
    if (std::find(availableSolvents.begin(), availableSolvents.end(), solvent) == availableSolvents.end()) {
      _cleanDataStructures(env, calc, res, mol);
      throw std::runtime_error("The given solvent is not available for implicit solvation within GFN2.");
    }
    double temp = _settings.getDouble("temperature");
    int state = 3;  // 1 bar of ideal gas and 1 mol/L of liquid solution
    int grid = 230; // n_grid_points, xtb default value
    xtb_setSolvent(env, calc, &solvent[0], &state, &temp, &grid);
  }

  // Run XTB singlepoint
  try {
    xtb_singlepoint(env, mol, calc, res);
  }
  catch (...) {
    _cleanDataStructures(env, calc, res, mol);
    throw Core::UnsuccessfulCalculationException("Xtb calculation failed:\n" +
                                                 boost::current_exception_diagnostic_information());
  }
  if (xtb_checkEnvironment(env) != 0) {
    // necessary raw pointers for xtb wrapper
    const int buffersize = 512;
    char error[buffersize] = "";
    xtb_getError(env, &error[0], &buffersize);
    std::string errorMessage(error);
    xtb_showEnvironment(env, nullptr);
    _cleanDataStructures(env, calc, res, mol);
    throw Core::UnsuccessfulCalculationException("Xtb calculation failed:\n" + errorMessage);
  }

  // Parse output
  this->_results = Scine::Utils::Results();
  // - Energy
  double energy = 0.0;
  xtb_getEnergy(env, res, &energy);
  if (xtb_checkEnvironment(env) != 0) {
    xtb_showEnvironment(env, nullptr);
    this->_results.set<Scine::Utils::Property::SuccessfulCalculation>(false);
    _cleanDataStructures(env, calc, res, mol);
    throw Core::UnsuccessfulCalculationException("Could not read XTB energy.");
  }
  this->_results.set<Scine::Utils::Property::Energy>(energy);
  // - Gradients
  if (_requiredProperties.containsSubSet(Scine::Utils::Property::Gradients)) {
    Utils::GradientCollection grad = Utils::GradientCollection::Zero(natoms, 3);
    xtb_getGradient(env, res, grad.data());
    if (xtb_checkEnvironment(env) != 0) {
      xtb_showEnvironment(env, nullptr);
      this->_results.set<Scine::Utils::Property::SuccessfulCalculation>(false);
      _cleanDataStructures(env, calc, res, mol);
      throw Core::UnsuccessfulCalculationException("Could not read XTB gradients.");
    }
    this->_results.set<Scine::Utils::Property::Gradients>(grad);
  }
  // - Bond orders
  if (_requiredProperties.containsSubSet(Scine::Utils::Property::BondOrderMatrix)) {
    Eigen::MatrixXd wbo = Eigen::MatrixXd::Zero(natoms, natoms);
    xtb_getBondOrders(env, res, wbo.data());
    if (xtb_checkEnvironment(env) != 0) {
      xtb_showEnvironment(env, nullptr);
      this->_results.set<Scine::Utils::Property::SuccessfulCalculation>(false);
      _cleanDataStructures(env, calc, res, mol);
      throw Core::UnsuccessfulCalculationException("Could not read XTB bond orders.");
    }
    Scine::Utils::BondOrderCollection bos(natoms);
    bos.setMatrix(wbo.sparseView(1e-12, 1.0));
    this->_results.set<Scine::Utils::Property::BondOrderMatrix>(bos);
  }
  // - Partial charges
  std::vector<double> q(natoms, 0.0);
  xtb_getCharges(env, res, q.data());
  if (xtb_checkEnvironment(env) != 0) {
    xtb_showEnvironment(env, nullptr);
    this->_results.set<Scine::Utils::Property::SuccessfulCalculation>(false);
    _cleanDataStructures(env, calc, res, mol);
    throw Core::UnsuccessfulCalculationException("Could not read XTB partial charges.");
  }
  this->_results.set<Scine::Utils::Property::AtomicCharges>(q);
  // - Occupation
  auto occupation = Scine::Utils::LcaoUtils::ElectronicOccupation();
  if (uhf == 0) {
    occupation.fillLowestRestrictedOrbitalsWithElectrons(attyp.sum() - static_cast<int>(charge));
  }
  else {
    int alpha = (attyp.sum() - static_cast<int>(charge) + uhf) / 2;
    int beta = (attyp.sum() - static_cast<int>(charge) - uhf) / 2;
    occupation.fillLowestUnrestrictedOrbitals(alpha, beta);
  }
  this->_results.set<Scine::Utils::Property::ElectronicOccupation>(occupation);
  //  - Dipole
  // Eigen::VectorXd dipole = Eigen::VectorXd::Zero(3);
  // xtb_getDipole(env, res, dipole.data());
  // Calculate Hessian
  if (_requiredProperties.containsSubSet(Scine::Utils::Property::Hessian) or
      _requiredProperties.containsSubSet(Scine::Utils::Property::Thermochemistry)) {
    Utils::NumericalHessianCalculator hessianCalculator(*this);
    auto numericalResult = hessianCalculator.calculate();
    this->_results.set<Utils::Property::Hessian>(numericalResult.take<Utils::Property::Hessian>());
  }

  // set successful to be able to autocomplete thermochemistry
  this->_results.set<Scine::Utils::Property::SuccessfulCalculation>(true);
  _settings.modifyString(Utils::SettingsNames::spinMode,
                         Utils::SpinModeInterpreter::getStringFromSpinMode(Utils::SpinMode::RestrictedOpenShell));
  this->_results.set<Scine::Utils::Property::ProgramName>(program);

  // - Thermochemistry
  if (_requiredProperties.containsSubSet(Scine::Utils::Property::Hessian) or
      _requiredProperties.containsSubSet(Scine::Utils::Property::Thermochemistry)) {
    Scine::Utils::ResultsAutoCompleter completer(*_structure);
    completer.setTemperature(_settings.getDouble(Utils::SettingsNames::temperature));
    completer.setPressure(_settings.getDouble(Utils::SettingsNames::pressure));
    completer.setMolecularSymmetryNumber(_settings.getInt(Utils::SettingsNames::symmetryNumber));
    completer.addOneWantedProperty(Scine::Utils::Property::Thermochemistry);
    completer.generateProperties(this->_results, *_structure);
  }

  // Clean XTB data structures
  _cleanDataStructures(env, calc, res, mol);

  return this->_results;
}

} /* namespace Xtb */
} /* namespace Scine */
