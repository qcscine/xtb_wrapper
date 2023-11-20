__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import scine_utilities as utils
from distutils import ccompiler

manager = utils.core.ModuleManager.get_instance()
if not manager.module_loaded('Xtb'):
    shlib_suffix = ccompiler.new_compiler().shared_lib_extension
    module_filename = "xtb.module" + shlib_suffix
    # Look within the python module directory (module is here in the case of
    # python packages) and the lib folder the site packages are in
    current_path = os.path.dirname(os.path.realpath(__file__))
    lib_path = os.path.dirname(os.path.dirname(os.path.dirname(current_path)))
    test_paths = [current_path, lib_path]

    def exists_and_could_load(path):
        full_path = os.path.join(path, module_filename)
        if os.path.exists(full_path):
            try:
                manager.load(full_path)
            except RuntimeError as err:
                print("Could not load {}: {}".format(full_path, err))
                return False

            return True

        return False

    if not any(map(exists_and_could_load, test_paths)):
        raise ImportError('{} could not be located.'.format(module_filename))
