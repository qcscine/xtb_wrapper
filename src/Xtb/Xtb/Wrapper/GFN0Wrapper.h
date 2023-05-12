/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef XTB_GFN0WRAPPER_H_
#define XTB_GFN0WRAPPER_H_

/* Internal Includes */
#include "Xtb/Wrapper/XtbCalculatorBase.h"
/* External Includes*/
#include <Utils/CalculatorBasics.h>
#include <Utils/Technical/CloneInterface.h>
#include <mutex>

namespace Scine {
namespace Xtb {

/**
 * @brief The SCINE Calculator for molecular GFN0 Calculations.
 */
class GFN0Wrapper : public Scine::Utils::CloneInterface<GFN0Wrapper, XtbCalculatorBase, Scine::Core::Calculator> {
 public:
  static constexpr const char* model = "GFN0";
  static constexpr const char* program = "Xtb";
  /// @brief Default Constructor
  GFN0Wrapper();
  /// @brief Default Destructor.
  ~GFN0Wrapper() = default;
  /// @brief Copy Constructor.
  GFN0Wrapper(const GFN0Wrapper& other) = default;
  /**
   * @brief Getter for the name of the underlying method.
   * @returns Returns the name of the underlying method.
   */
  std::string method() const final {
    return GFN0Wrapper::model;
  };
  /**
   * @brief Getter for the name of the underlying calculator.
   * @returns Returns the name of the underlying calculator.
   */
  std::string name() const final {
    return "XtbGFN0Calculator";
  };
  /**
   * @brief Report the possible properties.
   * @return Scine::Utils::PropertyList
   */
  Scine::Utils::PropertyList possibleProperties() const final {
    return Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian |
           Utils::Property::SuccessfulCalculation | Utils::Property::Thermochemistry;
  };
  /**
   * @brief Check if the method family is supported by this calculator.
   * @param methodFamily The method family as all caps string.
   * @return true  If it is supported.
   * @return false If it is not supported.
   */
  bool supportsMethodFamily(const std::string& methodFamily) const final {
    return methodFamily == "GFN0";
  }
  /**
   * @brief Run the calculation.
   * @param dummy
   * @return const Scine::Utils::Results& The results.
   */
  virtual const Scine::Utils::Results& calculate(std::string dummy) final;

 private:
  static std::mutex _mtx;
};

} /* namespace Xtb */
} /* namespace Scine */

#endif /* XTB_GFN0WRAPPER_H_ */
