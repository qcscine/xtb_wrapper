cmake_minimum_required(VERSION 3.9)

# Setup the Scine Xtb Project
project(Xtb
  VERSION 1.0.1
  DESCRIPTION "SCINE XTB Wrapper."
)

# Missing modeled dependencies of xtb
enable_language(Fortran)

# Set the module path for universal cmake files inclusion
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/dev/cmake)

# Find xtb
find_package(xtb QUIET)
if(NOT TARGET lib-xtb-static)
  include(DownloadProject)
  download_project(
    PROJ xtb
    GIT_REPOSITORY https://github.com/grimme-lab/xtb.git
    GIT_TAG v6.4.1
    QUIET
    UPDATE_DISCONNECTED 1
  )
  add_subdirectory(${xtb_SOURCE_DIR} ${xtb_BINARY_DIR})
  if(TARGET lib-xtb-static)
    message(STATUS "xtb was not found in your PATH, so it was downloaded.")
  else()
    string(CONCAT error_msg
      "xtb was not found in your PATH and could not be established through "
      "a download. Try specifying xtb_DIR or altering CMAKE_PREFIX_PATH to "
      "point to a candidate xtb installation base directory."
    )
    message(FATAL_ERROR ${error_msg})
  endif()
endif()

# Component setup
include(ComponentSetup)
scine_setup_component()

# Enable testing
if(SCINE_BUILD_TESTS)
  enable_testing()
endif()

# Subdirectories
add_subdirectory(src)

# Documentation
scine_component_documentation(UtilsOSDocumentation)
