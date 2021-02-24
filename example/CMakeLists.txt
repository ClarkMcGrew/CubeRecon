cmake_minimum_required(VERSION 3.0)

# Set the project name and version.  The project name will be used to
# import this into other CMAKE files.
project(CubeExample VERSION 1.0.0)

message("Example using the output of the CubeRecon Reconstruction -- ${VERSION}")

# Define the options that can be set in the cache, or on the cmake
# command line.
set(CMAKE_BUILD_TYPE Debug)

# Check to see if this is MACOS
if(APPLE)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif (APPLE)

# Make sure that ROOT is available.  ROOT is absolutely required.
find_package(ROOT REQUIRED
  COMPONENTS Geom EG Physics Matrix MathCore Tree RIO)
if(ROOT_FOUND)
  include(${ROOT_USE_FILE})
endif(ROOT_FOUND)

# CubeRecon needs the edep-sim I/O package.
find_package(CubeRecon)
if(CubeRecon_FOUND)
  message("CubeRecon found")
endif(CubeRecon_FOUND)

# Build the basic simulation
add_subdirectory(app)

