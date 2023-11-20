/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Internal Includes */
#include "XtbModule.h"
#include <Xtb/Wrapper/GFN0Wrapper.h>
#include <Xtb/Wrapper/GFN1Wrapper.h>
#include <Xtb/Wrapper/GFN2Wrapper.h>
#include <Xtb/Wrapper/GFNFFWrapper.h>
/* External Includes */
#include <Core/DerivedModule.h>
#include <Core/Exceptions.h>
#include <Utils/Settings.h>

namespace Scine {
namespace Xtb {

using InterfaceModelMap =
    boost::mpl::map<boost::mpl::pair<Core::Calculator, boost::mpl::vector<Xtb::GFN0Wrapper, Xtb::GFN1Wrapper, Xtb::GFN2Wrapper, Xtb::GFNFFWrapper>>>;

std::string XtbModule::name() const noexcept {
  return "Xtb";
}

boost::any XtbModule::get(const std::string& interface, const std::string& model) const {
  boost::any resolved = Core::DerivedModule::resolve<InterfaceModelMap>(interface, model);

  // Throw an exception if we could not match an interface or model
  if (resolved.empty()) {
    throw Core::ClassNotImplementedError();
  }

  return resolved;
}

bool XtbModule::has(const std::string& interface, const std::string& model) const noexcept {
  return Core::DerivedModule::has<InterfaceModelMap>(interface, model);
}

std::vector<std::string> XtbModule::announceInterfaces() const noexcept {
  return Core::DerivedModule::announceInterfaces<InterfaceModelMap>();
}

std::vector<std::string> XtbModule::announceModels(const std::string& interface) const noexcept {
  return Core::DerivedModule::announceModels<InterfaceModelMap>(interface);
}

std::shared_ptr<Core::Module> XtbModule::make() {
  return std::make_shared<XtbModule>();
}

std::vector<std::shared_ptr<Core::Module>> moduleFactory() {
  return {XtbModule::make()};
}

} /* namespace Xtb */
} /* namespace Scine */
