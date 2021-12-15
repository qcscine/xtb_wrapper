/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef XTB_XTBSETTINGS_H_
#define XTB_XTBSETTINGS_H_

/* External Includes */
#include <Utils/Settings.h>

namespace Scine {
namespace Xtb {

/**
 * @class
 * @brief The SCINE State for Xtb.
 */
class XtbSettings : public Scine::Utils::Settings {
 public:
  /**
   * @brief Construct a new XtbSettings object.
   */
  XtbSettings();
};

} /* namespace Xtb */
} /* namespace Scine */

#endif /* XTB_XTBSETTINGS_H_ */
