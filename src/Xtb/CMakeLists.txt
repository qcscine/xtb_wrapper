# Load lists of source files
include(Files.cmake)

# Find/Import dependecies
include(ImportUtilsOS)
import_utils_os()
include(ImportCore)
import_core()

add_library(Xtb SHARED ${XTB_MODULE_FILES})
set_target_properties(Xtb PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(Xtb PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

if(WIN32)
  set_target_properties(Xtb PROPERTIES
    PREFIX ""
    IMPORT_PREFIX ""
    WINDOWS_EXPORT_ALL_SYMBOLS ON
  )
endif()

set_target_properties(Xtb PROPERTIES
  OUTPUT_NAME xtb.module
  PREFIX ""
)

target_link_libraries(Xtb
  PRIVATE
    Scine::UtilsOS
    xtb::xtb
    gfortran
  PUBLIC
    Scine::CoreHeaders
)
set_target_properties(Xtb PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")

# Add namespaced aliases
add_library(Scine::Xtb ALIAS Xtb)
add_library(Scine::XtbModule ALIAS Xtb)

# Install
install(
  TARGETS Xtb
  EXPORT scineXtbTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

# Python Bindings
if(SCINE_BUILD_PYTHON_BINDINGS)
  # Alternative strategy for finding the python executable and library to
  # prefer virtualenvs over system paths.
  if(NOT DEFINED Python_FIND_STRATEGY)
    set(Python_FIND_STRATEGY LOCATION)
  endif()
  if(NOT DEFINED Python_FIND_REGISTRY)
    set(Python_FIND_REGISTRY LAST)
  endif()
  if(NOT DEFINED Python_FIND_VIRTUALENV)
    set(Python_FIND_VIRTUALENV STANDARD)
  endif()
  find_package(Python COMPONENTS Interpreter Development)

  # Copy the module into the python package
  add_custom_command(TARGET Xtb POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:Xtb> ${CMAKE_CURRENT_BINARY_DIR}/scine_xtb_wrapper
    COMMENT "Copying 'xtb.module.so' to 'scine_xtb_wrapper'"
  )

  # Python_ set by FindPython and PYTHON_ set by user via -DPYTHON_EXECUTABLE or set by FindPythonInterp, which is deprecated
  # Prefer PYTHON over Python
  if(Python_EXECUTABLE AND NOT PYTHON_EXECUTABLE)
    set(PYTHON_EXECUTABLE ${Python_EXECUTABLE})
  endif()

  if(NOT PYTHON_EXECUTABLE)
    message(FATAL_ERROR # executable variable is empty, give proper error here instead at install
            "Unable to find Python interpreter, required for Python Bindings. Please install Python or specify the PYTHON_EXECUTABLE CMake variable.")
  endif()

  install(CODE
  "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --prefix=${CMAKE_INSTALL_PREFIX} --upgrade --no-deps ${CMAKE_CURRENT_BINARY_DIR}
                   RESULT_VARIABLE retcode)
   if(NOT \$\{retcode\} EQUAL 0)
     message(FATAL_ERROR \"Fatal error when installing Python module using PIP.\")
   endif()"
  )

  # Copy the setup file and README.rst into the binary dir so you can install
  # using pip or make a distribution with the compiled binary
  file(
    COPY ${CMAKE_CURRENT_SOURCE_DIR}/Python/__init__.py
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scine_xtb_wrapper
  )
  add_custom_command(TARGET Xtb POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:Xtb> ${CMAKE_CURRENT_BINARY_DIR}/scine_xtb_wrapper
  )
  set(_module_name "xtb.module${CMAKE_SHARED_LIBRARY_SUFFIX}")
  include(TargetLibName)
  target_lib_type(Scine::UtilsOS _utils_libtype)
  if (_utils_libtype STREQUAL "SHARED_LIBRARY")
    if(APPLE)
      set_target_properties(Xtb PROPERTIES
        BUILD_WITH_INSTALL_RPATH ON
        INSTALL_RPATH "@loader_path;@loader_path/../lib"
      )
    elseif(UNIX)
      set_target_properties(Xtb PROPERTIES
        BUILD_WITH_INSTALL_RPATH ON
        INSTALL_RPATH "\$ORIGIN;\$ORIGIN/../lib"
      )
    endif()
    add_custom_command(TARGET Xtb POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:Scine::UtilsOS> ${CMAKE_CURRENT_BINARY_DIR}/scine_xtb_wrapper
      COMMENT "Copying dependent shared libraries into python package directory"
    )
    target_lib_filename(Scine::UtilsOS _utils_name)
    set(xtb_PY_DEPS ", \"${_utils_name}\", \"${_module_name}\"")
    unset(_utils_name)
  else()
    set(xtb_PY_DEPS ", \"${_module_name}\"")
  endif()
  unset(_utils_libtype)

  # Add setup, readme and licenses
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Python/setup.py
    ${CMAKE_CURRENT_BINARY_DIR}/setup.py
    @ONLY
  )
  unset(xtb_PY_DEPS)
  file(
    COPY ${PROJECT_SOURCE_DIR}/README.rst
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )
  file(COPY ${PROJECT_SOURCE_DIR}/LICENSE.txt ${EXTERNAL_LICENSE_FILES}
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scine_xtb_wrapper
  )
endif()

# Component-level cmake package files
scine_install_component_cmake_files(
  COMPONENT ScineXtb
  EXPORT_NAME scineXtbTargets
)
