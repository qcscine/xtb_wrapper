/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef XTB_XTBCALCULATORBASE_H_
#define XTB_XTBCALCULATORBASE_H_

/* Internal Includes */
#include "Xtb/Wrapper/XtbSettings.h"

/* External Includes */
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Geometry/ElementTypes.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <xtb.h>

namespace Scine {

namespace Utils {} // namespace Utils

namespace Xtb {

/**
 * @class
 * @brief The SCINE Calculator Base for all Xtb Calculators.
 *
 * This class comprises some of the common mapping features required in order to
 * couple XTB and the SCINE interfaces.
 */
class XtbCalculatorBase
  : public Scine::Utils::CloneInterface<Scine::Utils::Abstract<XtbCalculatorBase>, Scine::Core::Calculator> {
 public:
  /// @brief Default Constructor
  XtbCalculatorBase() = default;
  /// @brief Default Destructor.
  ~XtbCalculatorBase() override = default;
  /// @brief Copy Constructor.
  XtbCalculatorBase(const XtbCalculatorBase& other);
  /**
   * @brief Sets new structure and initializes the underlying method with the
   * parameter given in the settings.
   * @param structure The structure to be assigned.
   */
  void setStructure(const Scine::Utils::AtomCollection& structure) final;
  /**
   * @brief Getter for the underlying element types and positions.
   */
  std::unique_ptr<Scine::Utils::AtomCollection> getStructure() const final;
  /**
   * @brief Allows to modify the positions of the underlying
   * Utils::AtomCollection
   * @param newPositions the new positions to be assigned to the underlying
   * Utils::AtomCollection
   */
  void modifyPositions(Scine::Utils::PositionCollection newPositions) final;
  /**
   * @brief Getter for the coordinates of the underlying Utils::AtomCollection
   */
  const Scine::Utils::PositionCollection& getPositions() const final;
  /**
   * @brief Sets the properties to calculate.
   * @param requiredProperties a Scine::Utils::PropertyList object, contains an
   * enum class that work as a bitset, switching on and off the bits
   * representing a property.
   */
  void setRequiredProperties(const Scine::Utils::PropertyList& requiredProperties) final;
  /**
   * @brief Gets the properties to calculate.
   * @return requiredProperties a Scine::Utils::PropertyList object, contains an
   * enum class that work as a bitset, switching on and off the bits
   * representing a property.
   */
  virtual Scine::Utils::PropertyList getRequiredProperties() const final;
  /**
   * @brief Getter for the name of the underlying calculator.
   * @returns Returns the name of the underlying calculator.
   */
  virtual std::string name() const = 0;
  /**
   * @brief Getter for the name of the underlying method.
   * @returns Returns the name of the underlying method.
   */
  virtual std::string method() const = 0;
  /**
   * @brief Returns a list of all properties this Calculator can produce.
   * @return Scine::Utils::PropertyList A list of all properties this Calculator
   * can produce.
   */
  virtual Scine::Utils::PropertyList possibleProperties() const = 0;
  /**
   * @brief The main function running calculations.
   * @param dummy   A dummy parameter.
   * @return Scine::Utils::Results Return the result of the calculation.
   */
  virtual const Scine::Utils::Results& calculate(std::string dummy) = 0;
  /**
   * @brief Accessor for the Settings used in this method wrapper.
   * @returns Scine::Utils::Settings& The Settings.
   */
  Scine::Utils::Settings& settings() final;
  /**
   * @brief Const accessor for the Settings used in this method wrapper.
   * @returns const Scine::Utils::Settings& The Settings.
   */
  const Scine::Utils::Settings& settings() const final;
  /**
   * @brief Accessor for the Results stored in this method wrapper.
   * @returns Scine::Utils::Results& The results of the previous calculation.
   */
  Scine::Utils::Results& results() final;
  /**
   * @brief Const accessor for the Results used in this method wrapper.
   * @returns const Scine::Utils::Results& The results of the previous
   * calculation.
   */
  const Scine::Utils::Results& results() const final;
  /**
   * @brief Exchange the current state/system for a different one.
   * @param state The new state/system.
   */
  void loadState(std::shared_ptr<Scine::Core::State> state) final;
  /**
   * @brief Get a copy of current state/system.
   * @return std::shared_ptr<Scine::Core::State> The current state/system.
   */
  std::shared_ptr<Scine::Core::State> getState() const final;

  /**
   * @brief Checks charge and spin multiplicity in settings to be a valid input for the Xtb Wrapper
   * @throws std::runtime_error for wrong input of charge or multiplicity
   */
  void verifyPesValidity();
  /**
   * @brief Whether the calculator has no underlying Python code and can therefore
   * release the global interpreter lock in Python bindings
   */
  bool allowsPythonGILRelease() const override {
    return true;
  };

 protected:
  XtbSettings _settings;
  Scine::Utils::Results _results;
  Scine::Utils::PropertyList _requiredProperties;
  std::unique_ptr<Scine::Utils::AtomCollection> _structure;
  std::vector<std::string> _availableSolvationModels = std::vector<std::string>{"gbsa"};
  bool _setExternalCharges = false;
  void _cleanDataStructures(xtb_TEnvironment& env, xtb_TCalculator& calc, xtb_TResults& res, xtb_TMolecule& mol);
  void setExternalCharges(xtb_TEnvironment& env, xtb_TCalculator& calc, xtb_TResults& res, xtb_TMolecule& mol);
  std::map<Utils::ElementType, std::pair<int, int>> _nElectronsAndAos = {
      {Utils::ElementType::H, {1, 1}},   {Utils::ElementType::He, {2, 4}},  {Utils::ElementType::Li, {1, 4}},
      {Utils::ElementType::Be, {2, 4}},  {Utils::ElementType::B, {3, 4}},   {Utils::ElementType::C, {4, 4}},
      {Utils::ElementType::N, {5, 4}},   {Utils::ElementType::O, {6, 4}},   {Utils::ElementType::F, {7, 4}},
      {Utils::ElementType::Ne, {8, 9}},  {Utils::ElementType::Na, {1, 4}},  {Utils::ElementType::Mg, {2, 9}},
      {Utils::ElementType::Al, {3, 9}},  {Utils::ElementType::Si, {4, 9}},  {Utils::ElementType::P, {5, 9}},
      {Utils::ElementType::S, {6, 9}},   {Utils::ElementType::Cl, {7, 9}},  {Utils::ElementType::Ar, {8, 9}},
      {Utils::ElementType::K, {1, 4}},   {Utils::ElementType::Ca, {2, 9}},  {Utils::ElementType::Sc, {3, 9}},
      {Utils::ElementType::Ti, {4, 9}},  {Utils::ElementType::V, {5, 9}},   {Utils::ElementType::Cr, {6, 9}},
      {Utils::ElementType::Mn, {7, 9}},  {Utils::ElementType::Fe, {8, 9}},  {Utils::ElementType::Co, {9, 9}},
      {Utils::ElementType::Ni, {10, 9}}, {Utils::ElementType::Cu, {11, 9}}, {Utils::ElementType::Zn, {2, 4}},
      {Utils::ElementType::Ga, {3, 9}},  {Utils::ElementType::Ge, {4, 9}},  {Utils::ElementType::As, {5, 9}},
      {Utils::ElementType::Se, {6, 9}},  {Utils::ElementType::Br, {7, 9}},  {Utils::ElementType::Kr, {8, 9}},
      {Utils::ElementType::Rb, {1, 4}},  {Utils::ElementType::Sr, {2, 9}},  {Utils::ElementType::Y, {3, 9}},
      {Utils::ElementType::Zr, {4, 9}},  {Utils::ElementType::Nb, {5, 9}},  {Utils::ElementType::Mo, {6, 9}},
      {Utils::ElementType::Tc, {7, 9}},  {Utils::ElementType::Ru, {8, 9}},  {Utils::ElementType::Rh, {9, 9}},
      {Utils::ElementType::Pd, {10, 9}}, {Utils::ElementType::Ag, {11, 9}}, {Utils::ElementType::Cd, {2, 4}},
      {Utils::ElementType::In, {3, 9}},  {Utils::ElementType::Sn, {4, 9}},  {Utils::ElementType::Sb, {5, 9}},
      {Utils::ElementType::Te, {6, 9}},  {Utils::ElementType::I, {7, 9}},   {Utils::ElementType::Xe, {8, 9}},
      {Utils::ElementType::Cs, {1, 4}},  {Utils::ElementType::Ba, {2, 9}},  {Utils::ElementType::La, {3, 9}},
      {Utils::ElementType::Ce, {3, 9}},  {Utils::ElementType::Pr, {3, 9}},  {Utils::ElementType::Nd, {3, 9}},
      {Utils::ElementType::Pm, {3, 9}},  {Utils::ElementType::Sm, {3, 9}},  {Utils::ElementType::Eu, {3, 9}},
      {Utils::ElementType::Gd, {3, 9}},  {Utils::ElementType::Tb, {3, 9}},  {Utils::ElementType::Dy, {3, 9}},
      {Utils::ElementType::Ho, {3, 9}},  {Utils::ElementType::Er, {3, 9}},  {Utils::ElementType::Tm, {3, 9}},
      {Utils::ElementType::Yb, {3, 9}},  {Utils::ElementType::Lu, {3, 9}},  {Utils::ElementType::Hf, {4, 9}},
      {Utils::ElementType::Ta, {5, 9}},  {Utils::ElementType::W, {6, 9}},   {Utils::ElementType::Re, {7, 9}},
      {Utils::ElementType::Os, {8, 9}},  {Utils::ElementType::Ir, {9, 9}},  {Utils::ElementType::Pt, {10, 9}},
      {Utils::ElementType::Au, {11, 9}}, {Utils::ElementType::Hg, {2, 4}},  {Utils::ElementType::Tl, {3, 4}},
      {Utils::ElementType::Pb, {4, 4}},  {Utils::ElementType::Bi, {5, 4}},  {Utils::ElementType::Po, {6, 4}},
      {Utils::ElementType::At, {7, 9}},  {Utils::ElementType::Rn, {8, 9}}};
};

} /* namespace Xtb */
} /* namespace Scine */

#endif /* XTB_XTBCALCULATORBASE_H_ */
