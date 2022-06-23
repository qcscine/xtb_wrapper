/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Internal Includes */
#include "Xtb/Wrapper/XtbCalculatorBase.h"
#include "Xtb/Wrapper/XtbState.h"
#include <algorithm>
#include <cctype>
#include <string>

namespace Scine {
namespace Xtb {

XtbCalculatorBase::XtbCalculatorBase(const XtbCalculatorBase& other) {
  _settings = other._settings;
  _results = other._results;
  _requiredProperties = other._requiredProperties;
  _structure = std::make_unique<Scine::Utils::AtomCollection>(*(other._structure));
}

void XtbCalculatorBase::setStructure(const Scine::Utils::AtomCollection& structure) {
  _structure = std::make_unique<Scine::Utils::AtomCollection>(structure);
  this->_results = Scine::Utils::Results();
}

std::unique_ptr<Scine::Utils::AtomCollection> XtbCalculatorBase::getStructure() const {
  return std::make_unique<Scine::Utils::AtomCollection>(*_structure);
}

void XtbCalculatorBase::modifyPositions(Scine::Utils::PositionCollection newPositions) {
  if (!this->_structure) {
    throw std::runtime_error("Failed to modify non existing structure.");
  }
  _structure->setPositions(newPositions);
  this->_results = Scine::Utils::Results();
}

const Scine::Utils::PositionCollection& XtbCalculatorBase::getPositions() const {
  if (!this->_structure) {
    throw std::runtime_error("Failed to get non existing structure.");
  }
  return _structure->getPositions();
}

void XtbCalculatorBase::setRequiredProperties(const Scine::Utils::PropertyList& requiredProperties) {
  if (!this->possibleProperties().containsSubSet(requiredProperties)) {
    throw std::runtime_error("Unavailable Properties requested.");
  }
  _requiredProperties = requiredProperties;
}

Scine::Utils::PropertyList XtbCalculatorBase::getRequiredProperties() const {
  return _requiredProperties;
}

Scine::Utils::Settings& XtbCalculatorBase::settings() {
  return _settings;
}

const Scine::Utils::Settings& XtbCalculatorBase::settings() const {
  return _settings;
}

Scine::Utils::Results& XtbCalculatorBase::results() {
  return _results;
}

const Scine::Utils::Results& XtbCalculatorBase::results() const {
  return _results;
}

void XtbCalculatorBase::loadState(std::shared_ptr<Scine::Core::State> state) {
  auto castState = std::dynamic_pointer_cast<XtbState>(state);
  if (!castState)
    throw Scine::Core::StateCastingException();
  this->setStructure(castState->system);
}

std::shared_ptr<Scine::Core::State> XtbCalculatorBase::getState() const {
  return std::make_shared<XtbState>(*_structure);
}

void XtbCalculatorBase::verifyPesValidity() {
  if (!_structure) {
    throw std::runtime_error("The " + name() + " calculator does currently not hold a structure");
  }
  const int charge = _settings.getInt(Utils::SettingsNames::molecularCharge);
  const int multiplicity = _settings.getInt(Utils::SettingsNames::spinMultiplicity);

  // Check selected method
  std::string method = _settings.getString(Utils::SettingsNames::method);
  std::transform(method.begin(), method.end(), method.begin(), [](unsigned char c) { return std::tolower(c); });
  std::string model = this->method();
  std::transform(model.begin(), model.end(), model.begin(), [](unsigned char c) { return std::tolower(c); });
  if (method != model && method != "any") {
    throw std::runtime_error("The " + name() + " calculator does not provide the requested method.");
  }
  _settings.modifyString(Utils::SettingsNames::method, model);

  // get n electrons for uncharged species and available AOs
  int nElectrons = 0;
  int nAos = 0;
  for (int i = 0; i < _structure->size(); ++i) {
    const auto ele = _structure->getElement(i);
    // check if bigger element number than the last included one
    if (Utils::ElementInfo::Z(ele) > Utils::ElementInfo::Z(_nElectronsAndAos.rbegin()->first)) {
      throw std::runtime_error(
          "XTB: The structure includes an element that is not supported by the GFN-X method family.");
    }
    auto parameters = _nElectronsAndAos.at(ele);
    nElectrons += parameters.first;
    nAos += parameters.second;
  }

  // check charge
  if (charge >= nElectrons) {
    throw std::runtime_error("XTB: The chosen molecular charge (" + std::to_string(charge) + ") is too positive for " +
                             std::to_string(nElectrons) + " electrons.");
  }
  else if (nElectrons - charge > 2 * nAos) {
    throw std::runtime_error("XTB: Not enough orbitals to accommodate the chosen molecular charge (" +
                             std::to_string(charge) + ").");
  }
  nElectrons -= charge;

  // check multiplicity
  if (multiplicity > nElectrons + 1) {
    throw std::runtime_error(
        "XTB: Molecular charge removes more electrons than are actually present in the calculation.");
  }
  int numberSpotsLeftForElectrons = 2 * nAos - nElectrons;
  if (multiplicity > numberSpotsLeftForElectrons + 1) {
    throw std::runtime_error("XTB: The chosen spin multiplicity (" + std::to_string(multiplicity) +
                             ") is too large (not enough orbitals).");
  }
  // Check that number of electrons and spin multiplicity are compatible
  else if ((multiplicity + nElectrons) % 2 == 0) {
    throw std::runtime_error("XTB: The chosen spin multiplicity (" + std::to_string(multiplicity) +
                             ") is not compatible with the molecular charge (" + std::to_string(charge) + ").");
  }
}

void XtbCalculatorBase::_cleanDataStructures(xtb_TEnvironment& env, xtb_TCalculator& calc, xtb_TResults& res,
                                             xtb_TMolecule& mol) {
  xtb_delResults(&res);
  xtb_delCalculator(&calc);
  xtb_delMolecule(&mol);
  xtb_delEnvironment(&env);
}

} /* namespace Xtb */
} /* namespace Scine */
