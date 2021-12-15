__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
See LICENSE.txt for details.
"""

import setuptools

# Read README.rst for the long description
with open("README.rst", "r") as fh:
    long_description = fh.read()


class EmptyListWithLength(list):
    """ Makes the wheel a binary distribution and platlib compliant. """

    def __len__(self):
        return 1


# Define the setup
setuptools.setup(
    name="scine_xtb_wrapper",
    version="@Xtb_VERSION@",
    author="ETH Zurich, Laboratory for Physical Chemistry, Reiher Group",
    author_email="scine@phys.chem.ethz.ch",
    description="A wrapper for xtb",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://www.scine.ethz.ch",
    packages=["scine_xtb_wrapper"],
    package_data={"scine_xtb_wrapper": ['*.txt' @xtb_PY_DEPS@]},
    install_requires=["scine_utilities"],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: C++",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    zip_safe=False,
    test_suite='pytest',
    tests_require=['pytest'],
    ext_modules=EmptyListWithLength()
)
