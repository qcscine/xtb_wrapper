/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef XTB_GFN2WRAPPER_H_
#define XTB_GFN2WRAPPER_H_

/* Internal Includes */
#include "Xtb/Wrapper/XtbCalculatorBase.h"
/* External Includes*/
#include <Utils/CalculatorBasics.h>
#include <Utils/Technical/CloneInterface.h>
#include <mutex>

namespace Scine {
namespace Xtb {

/**
 * @brief The SCINE Calculator for molecular GFN2 Calculations.
 */
class GFN2Wrapper : public Scine::Utils::CloneInterface<GFN2Wrapper, XtbCalculatorBase> {
 public:
  static constexpr const char* model = "GFN2";
  static constexpr const char* program = "Xtb";
  /// @brief Default Constructor
  GFN2Wrapper() = default;
  /// @brief Default Destructor.
  ~GFN2Wrapper() = default;
  /// @brief Copy Constructor.
  GFN2Wrapper(const GFN2Wrapper& other) = default;
  /**
   * @brief Getter for the name of the underlying method.
   * @returns Returns the name of the underlying method.
   */
  std::string name() const final {
    return "XtbGFN2Calculator";
  };
  /**
   * @brief Report the possible properties.
   * @return Scine::Utils::PropertyList
   */
  Scine::Utils::PropertyList possibleProperties() const final {
    return Scine::Utils::Property::Energy | Utils::Property::AtomicCharges | Scine::Utils::Property::Gradients |
           Scine::Utils::Property::Hessian | Scine::Utils::Property::BondOrderMatrix | Scine::Utils::Property::Thermochemistry;
  };
  /**
   * @brief Check if the method family is supported by this calculator.
   * @param methodFamily The method family as all caps string.
   * @return true  If it is supported.
   * @return false If it is not supported.
   */
  bool supportsMethodFamily(const std::string& methodFamily) const final {
    return methodFamily == "GFN2";
  }
  /**
   * @brief Run the calculation.
   * @param dummy
   * @return const Scine::Utils::Results& The results.
   */
  const Scine::Utils::Results& calculate(std::string dummy) final;

 private:
  static std::mutex _mtx;
};

} /* namespace Xtb */
} /* namespace Scine */

#endif /* XTB_GFN2WRAPPER_H_ */
