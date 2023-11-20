/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef XTB_XTBSTATE_H_
#define XTB_XTBSTATE_H_

/* External Includes */
#include <Core/BaseClasses/StateHandableObject.h>
#include <Utils/Geometry/AtomCollection.h>

namespace Scine {
namespace Xtb {

/**
 * @class
 * @brief The SCINE State for Xtb.
 */
class XtbState : public Scine::Core::State {
 public:
  /**
   * @brief Construct a new XtbState.
   * @param s The system represented by the atoms.
   */
  XtbState(const Scine::Utils::AtomCollection& s) : system(s){};
  ~XtbState() = default;
  Scine::Utils::AtomCollection system;
};

} /* namespace Xtb */
} /* namespace Scine */

#endif /* XTB_XTBSTATE_H_ */
