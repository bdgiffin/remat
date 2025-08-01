#------------------------------------------------------------------------------
# REMAT: Bit-reversible physics library
#------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.14)

project(REMAT LANGUAGES C CXX)

#------------------------------------------------------------------------------
# Setup BLT
#------------------------------------------------------------------------------
# Set BLT_SOURCE_DIR to default location, if not set by user 
if(DEFINED BLT_SOURCE_DIR)
    # Support having a shared BLT outside of the repository if given a BLT_SOURCE_DIR
    if(NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR "Given BLT_SOURCE_DIR does not contain SetupBLT.cmake")
    endif()
else()
    set(BLT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/cmake/blt" CACHE PATH "")
    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR
            "The BLT is not present. "
            "Either run the following two commands in your git repository: \n"
            "    git submodule init\n"
            "    git submodule update\n"
            "Or add -DBLT_SOURCE_DIR=/path/to/blt to your CMake command." )
    endif()
endif()

# Tune BLT to our needs
if (NOT BLT_CXX_STD)
    set(BLT_CXX_STD "c++17" CACHE STRING "")
endif()

set(ENABLE_ALL_WARNINGS ON CACHE BOOL "")
#set(ENABLE_WARNINGS_AS_ERRORS ON CACHE BOOL "")

include(${BLT_SOURCE_DIR}/SetupBLT.cmake)

#------------------------------------------------------------------------------
# For Emscripten builds, enable linking for SHARED .wasm library targets
#------------------------------------------------------------------------------

# Fix from: https://github.com/emscripten-core/emscripten/issues/15276#issuecomment-1039349267
if (EMSCRIPTEN)
   set_property(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)
   set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "-s SIDE_MODULE=1")
   set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "-s SIDE_MODULE=1")
   set(CMAKE_STRIP FALSE)  # used by default in pybind11 on .so modules
endif ()

#------------------------------------------------------------------------------
# Setup Macros and dependencies
#------------------------------------------------------------------------------

include(${PROJECT_SOURCE_DIR}/cmake/tpl/SetupTPL.cmake)

#------------------------------------------------------------------------------
# Build src
#------------------------------------------------------------------------------

add_subdirectory(src)

#------------------------------------------------------------------------------
# Configure examples
#------------------------------------------------------------------------------

add_subdirectory(examples)

#------------------------------------------------------------------------------
# Configure unit tests
#------------------------------------------------------------------------------

add_subdirectory(tests)

#------------------------------------------------------------------------------
# Configure default installation
#------------------------------------------------------------------------------

install(TARGETS REMAT
	LIBRARY DESTINATION lib)
install(TARGETS REMAT
	DESTINATION package)
install(FILES ${REMAT_python}
	DESTINATION package)
install(FILES ${REMAT_headers}
	DESTINATION include)
message(STATUS "REMAT_python:  " ${REMAT_python})
message(STATUS "REMAT_headers: " ${REMAT_headers})

# Sets the installation prefix (where to install files)
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/install)
