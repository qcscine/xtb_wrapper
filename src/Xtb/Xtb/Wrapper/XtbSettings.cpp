/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Internal Includes */
#include "Xtb/Wrapper/XtbSettings.h"
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <xtb.h>

/* External Includes */
#include <cstring>
#if defined(_OPENMP)
#  include <omp.h>
#endif

namespace Scine {
namespace Xtb {

XtbSettings::XtbSettings() : Scine::Utils::Settings("XtbSettings") {
  using namespace Scine::Utils;
  using namespace Scine::Utils::UniversalSettings;

  // Method
  StringDescriptor method("The method to be used.");
  method.setDefaultValue("");
  this->_fields.push_back(SettingsNames::method, method);

  // Multiplicity
  IntDescriptor spin_multiplicity("The multiplicity.");
  spin_multiplicity.setMinimum(1);
  spin_multiplicity.setDefaultValue(1);
  this->_fields.push_back(SettingsNames::spinMultiplicity, spin_multiplicity);

  // Charge
  IntDescriptor molecular_charge("The molecular charge.");
  molecular_charge.setDefaultValue(0);
  this->_fields.push_back(SettingsNames::molecularCharge, molecular_charge);

  // Spin mode
  OptionListDescriptor spin_mode("The spin mode.");
  spin_mode.addOption(Utils::SpinModeInterpreter::getStringFromSpinMode(SpinMode::Any));
  spin_mode.addOption(Utils::SpinModeInterpreter::getStringFromSpinMode(SpinMode::RestrictedOpenShell));
  spin_mode.setDefaultOption(Utils::SpinModeInterpreter::getStringFromSpinMode(SpinMode::RestrictedOpenShell));
  this->_fields.push_back(SettingsNames::spinMode, spin_mode);

  // Print level
  IntDescriptor prlevel("The verbosity level of the XTB output.");
  prlevel.setMinimum(0);
  prlevel.setMaximum(2);
  prlevel.setDefaultValue(0);
  this->_fields.push_back("print_level", prlevel);

  // Accuracy
  DoubleDescriptor scfAcc("The energy accuracy used for XTB calculations. This settings automatically influences "
                          "integral cutoffs and wavefunction accuracy.");
  scfAcc.setDefaultValue(1e-7);
  this->_fields.push_back(SettingsNames::selfConsistenceCriterion, scfAcc);

  // Electronic temperature
  DoubleDescriptor etemp("The temperature that controls the extent of Fermi smearing.");
  etemp.setMinimum(0.0);
  etemp.setDefaultValue(300.0);
  this->_fields.push_back(SettingsNames::electronicTemperature, etemp);

  // Max SCF iterations
  IntDescriptor maxiter("The maximum number of SCF iterations.");
  maxiter.setMinimum(0);
  maxiter.setDefaultValue(100);
  this->_fields.push_back(SettingsNames::maxScfIterations, maxiter);

  // Solvent
  StringDescriptor solvent("The implicit solvent to be used.");
  solvent.setDefaultValue("");
  this->_fields.push_back(SettingsNames::solvent, solvent);

  // Solvation
  StringDescriptor solvation("The solvation model to be used.");
  solvation.setDefaultValue("");
  this->_fields.push_back(SettingsNames::solvation, solvation);

  // Temperature used for thermochemical calculations
  DoubleDescriptor thermTemp("The temperature used for the thermochemical calculation.");
  thermTemp.setMinimum(0.0);
  thermTemp.setDefaultValue(298.15);
  this->_fields.push_back(SettingsNames::temperature, thermTemp);

  // Symmetry number used for thermochemical calculations
  IntDescriptor symmetryNumber("The molecular symmetry number used for the thermochemical calculation.");
  symmetryNumber.setMinimum(1);
  symmetryNumber.setDefaultValue(1);
  this->_fields.push_back(SettingsNames::symmetryNumber, symmetryNumber);

  // Parallel execution
  IntDescriptor parallel("The maximum number of cores to be used.");
#if defined(_OPENMP)
#  pragma omp parallel
  { parallel.setDefaultValue(omp_get_num_threads()); }
#else
  parallel.setDefaultValue(1);
#endif
  this->_fields.push_back(SettingsNames::externalProgramNProcs, parallel);

  this->resetToDefaults();
}

} /* namespace Xtb */
} /* namespace Scine */
