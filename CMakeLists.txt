cmake_minimum_required(VERSION 3.9)

# Setup the Scine Xtb Project
project(Xtb
  VERSION 2.0.0
  DESCRIPTION "SCINE XTB Wrapper."
)

# Missing modeled dependencies of xtb
enable_language(Fortran)

# Set the module path for universal cmake files inclusion
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/dev/cmake)

# Find xtb
find_package(xtb QUIET)
if(NOT TARGET lib-xtb-static)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(XTB QUIET xtb)
  if(XTB_FOUND)
    add_library(xtb::xtb INTERFACE IMPORTED)
    target_link_libraries(
      xtb::xtb
      INTERFACE
      ${XTB_LINK_LIBRARIES}
    )
    target_include_directories(
      xtb::xtb
      INTERFACE
      ${XTB_INCLUDE_DIRS}
    )
  else()
    include(DownloadProject)
    download_project(
      PROJ xtb
      GIT_REPOSITORY https://github.com/grimme-lab/xtb.git
      GIT_TAG v6.5.1
      QUIET
      UPDATE_DISCONNECTED 1
    )
    add_subdirectory(${xtb_SOURCE_DIR} ${xtb_BINARY_DIR})
    if(TARGET lib-xtb-static)
      message(STATUS "xtb was not found in your PATH, so it was downloaded.")
    else()
      string(CONCAT error_msg
        "xtb was not found in your PATH and could not be downloaded."
        "Try specifying xtb_DIR or altering CMAKE_PREFIX_PATH to "
        "point to a candidate xtb installation base directory."
      )
      message(FATAL_ERROR ${error_msg})
    endif()
  endif()
endif()

if(${xtb_WITH_OpenMP})
  find_dependency("OpenMP" REQUIRED)
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
