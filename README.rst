SCINE - XTB Wrapper -- A Wrapper for Tight Binding
==================================================

Introduction
------------

SCINE XTB Wrapper is a wrapper around `xtb <https://github.com/grimme-lab/xtb>`_, it
exports:

- GFN0-xTB
- GFN1-xTB (formerly GFN-xTB)
- GFN2-xTB
- GFN-FF

into the SCINE tool chain.
Each method is represented by its own ``Calculator`` and the entire wrapper
constitutes a SCINE module that can be loaded dynamically at runtime.
For more information on these concepts see the ``Scine::Core``
`repository <https://github.com/qcscine/core>`_.

License and Copyright Information
---------------------------------

The SCINE XTB wrapper is distributed under the BSD 3-clause "New" or "Revised"
License. For more license and copyright information, see the file ``LICENSE.txt``
in the repository.

Note: this license does not cover the original `xtb` source code.
For the copyright information of the `xtb` code please follow the linked git
submodule to the developers' repository.

Installation and Usage
----------------------

The wrapper can be built and installed using the following commands::

    export INSTALL_PATH=<desired path>
    git submodule init
    git submodule update
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DSCINE_BUILD_PYTHON_BINDINGS=ON -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH ..
    make -j<number of cores to use>
    make install

This will generate and install both the main `xtb` code and the wrapper in the
form of the file ``xtb.module.so`` that can be used in SCINE.

In order to make XTB available to SCINE the following two environment variables
need to be set::

    export XTBPATH=$INSTALL_PATH/share/xtb
    export SCINE_MODULE_PATH=$SCINE_MODULE_PATH:$INSTALL_PATH/lib

Afterwards, SCINE programs such as `ReaDuct <https://github.com/qcscine/readuct>`_
will pick up XTB's existence and it will be possible to request the implemented
methods.

The SCINE XTB wrapper is also available via Python.
The underlying SCINE module can be loaded and its implemented calculators
accessed using the standard ``scine_utilities`` Python bindings.
A minimal workflow could look like this::

    import scine_utilities as su
    import scine_xtb_wrapper
    
    # Generate Structure
    structure = su.AtomCollection()
    structure.elements = [su.ElementType.H, su.ElementType.H]
    structure.positions = [[-0.7, 0, 0], [0.7, 0, 0]]
    
    # Get calculator
    calculator = su.core.get_calculator('GFN2', 'xtb')

    # Configure calculator
    calculator.structure = structure
    calculator.set_required_properties([su.Property.Energy, su.Property.Gradients])
    
    # Calculate
    results = calculator.calculate()
    print(results.energy)
    print(results.gradients)

How to Cite
-----------

When publishing results obtained with the SCINE XTB wrapper, please cite the corresponding
release as archived on `Zenodo <https://zenodo.org/record/5782861>`_ (DOI
10.5281/zenodo.5782861; please use the DOI of the respective release).

Furthermore, when publishing results obtained with any SCINE module, please cite the following paper:

T. Weymuth, J. P. Unsleber, P. L. Türtscher, M. Steiner, J.-G. Sobez, C. H. Müller, M. Mörchen,
V. Klasovita, S. A. Grimmel, M. Eckhoff, K.-S. Csizi, F. Bosia, M. Bensberg, M. Reiher,
"SCINE—Software for chemical interaction networks", *J. Chem. Phys.*, **2024**, *160*, 222501
(DOI `10.1063/5.0206974 <https://doi.org/10.1063/5.0206974>`_).

This wrapper should also not be mistaken for the actual `xtb` code it wraps.
For the latter code and its citations, we refer you to the original `xtb`
repository. There you will find the references of the actual methods used.
They are listed in the
`README.md <https://github.com/grimme-lab/xtb/blob/master/README.md>`_.

Support and Contact
-------------------

In case you should encounter problems or bugs with the wrapper, please write a
short message to scine@phys.chem.ethz.ch.
