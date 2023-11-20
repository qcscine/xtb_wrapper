from dev.conan.base import ScineConan


class ScineXtbConan(ScineConan):
    name = "scine_xtb_wrapper"
    version = "3.0.0"
    url = "https://github.com/qcscine/xtb_wrapper"
    description = """A wrapper around xtb (https://github.com/grimme-lab/xtb),
it exports GFN0-xTB, GFN1-xTB (formerly GFN-xTB), GFN2-xTB, and GFN-FF into the SCINE tool
chain."""
    options = {
        "shared": [True, False],
        "python": [True, False],
        "microarch": ["detect", "none"]
    }
    default_options = {
        "shared": True,
        "python": False,
        "microarch": "none"
    }
    exports = "dev/conan/*.py"
    exports_sources = [
        "dev/cmake/*", "src/*", "CMakeLists.txt", "README.rst", "LICENSE.txt",
        "dev/conan/hook.cmake", "dev/conan/glue/*"]
    build_requires = "xtb/6.5.1"
    requires = [
        "scine_utilities/9.0.0",
        "boost/[>=1.71.0]",
        "lapack/[>=3.7.1]"
    ]
    cmake_name = "Xtb"
