message("Cube Recon Tools")

# Set the default verbosity of the compilation...
set(CMAKE_VERBOSE_MAKEFILE TRUE)

# Configure the dependencies
find_package(ROOT REQUIRED
  COMPONENTS Geom Physics Matrix MathCore Tree RIO)
if(ROOT_FOUND)
  include(${ROOT_USE_FILE})
endif(ROOT_FOUND)

# Define the source and include files that should be used for the io
# part of CubeRecon.
set(source
  ToolPrimaryId.cxx ToolG4Hits.cxx ToolMainTrajectory.cxx ToolTrueDirection.cxx
  ToolCubeTruth.cxx ToolContained.cxx ToolRecon.cxx ToolInitialize.cxx
  )

set(includes
  ToolPrimaryId.hxx ToolG4Hits.hxx ToolMainTrajectory.hxx ToolTrueDirection.hxx
  ToolCubeTruth.hxx ToolContained.hxx ToolRecon.hxx ToolInitialize.hxx
  )

# Make sure the current directories are available for the later
# compilation.  This is required to make ROOT_GENERATE_DICTIONARY
# work.
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Build the library.
add_library(cuberecon_tools SHARED ${source} )

target_include_directories(cuberecon_tools PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:include/CubeRecon>")

target_link_libraries(cuberecon_tools
  PUBLIC cuberecon_io ${ROOT_LIBRARIES})

# Install the io library for CubeRecon
install(TARGETS cuberecon_tools
  EXPORT CubeReconTargets
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include/CubeRecon )

# Install the header files.
install(FILES ${includes} DESTINATION include/CubeRecon)
