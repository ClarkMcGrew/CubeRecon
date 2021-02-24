message("Cube Recon Library")

# Set the default verbosity of the compilation...
set(CMAKE_VERBOSE_MAKEFILE TRUE)

# Configure the dependencies
find_package(ROOT REQUIRED
  COMPONENTS Geom EG Physics Matrix MathCore Tree RIO)
if(ROOT_FOUND)
  include(${ROOT_USE_FILE})
endif(ROOT_FOUND)

# Define the source and include files that should be used for the library
# part of CubeRecon.
set(source
  CubeERepSim.cxx
  CubeHitUtilities.cxx CubeClusterManagement.cxx
  CubeCreateTrack.cxx  CubeMakeUsed.cxx
  CubeTrackFit.cxx CubePCATrackFit.cxx CubeStochTrackFit.cxx
  CubeTimeSlice.cxx CubeMakeHits3D.cxx CubeHits3D.cxx CubeShareCharge.cxx
  CubeRecon.cxx CubeCleanHits.cxx CubeClusterHits.cxx
  CubeTreeRecon.cxx CubeSpanningTree.cxx
  CubeFindKinks.cxx CubeGrowClusters.cxx CubeGrowTracks.cxx CubeMergeXTalk.cxx
  CubeBuildPairwiseVertices.cxx CubeVertexFit.cxx
  )

set(includes
  CubeERepSim.hxx
  CubeHitUtilities.hxx  CubeClusterManagement.hxx
  CubeSafeLine.hxx CubeCreateTrack.hxx CubeMakeUsed.hxx
  CubeTrackFit.hxx CubePCATrackFit.hxx CubeStochTrackFit.hxx
  CubeTimeSlice.hxx CubeMakeHits3D.hxx CubeHits3D.hxx CubeShareCharge.hxx
  CubeRecon.hxx CubeCleanHits.hxx CubeTreeRecon.hxx
  CubeClusterHits.hxx CubeSpanningTree.hxx
  CubeFindKinks.hxx CubeGrowClusters.hxx CubeGrowTracks.hxx  CubeMergeXTalk.hxx
  CubeBuildPairwiseVertices.hxx CubeVertexFit.hxx
  )

# Make sure the current directories are available for the later
# compilation.  This is required to make ROOT_GENERATE_DICTIONARY
# work.
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Build the library.
add_library(cuberecon SHARED ${source})

target_include_directories(cuberecon PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:include/CubeRecon>")

target_link_libraries(cuberecon PUBLIC cuberecon_io ${ROOT_LIBRARIES})

# Install the library for CubeRecon
install(TARGETS cuberecon
  EXPORT CubeReconTargets
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include/CubeRecon )

# Install the header files.
install(FILES ${includes} DESTINATION include/CubeRecon)
