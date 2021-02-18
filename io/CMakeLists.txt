message("Cube Recon I/O")

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
  CubeInfo.cxx CubeHandle.cxx CubeEvent.cxx
  CubeHit.cxx CubeHitSelection.cxx
  CubeCorrValues.cxx CubeReconState.cxx CubeVertexState.cxx
  CubeClusterState.cxx CubeShowerState.cxx CubeTrackState.cxx
  CubeReconObject.cxx CubeReconNode.cxx
  CubeReconVertex.cxx CubeReconCluster.cxx
  CubeReconShower.cxx CubeReconTrack.cxx
  CubeAlgorithm.cxx CubeAlgorithmResult.cxx
  CubeG4Hit.cxx CubeG4Trajectory.cxx
  TUnitsTable.cxx
  )

set(includes
  CubeInfo.hxx CubeLog.hxx CubeHandle.hxx CubeEvent.hxx
  CubeHit.hxx CubeHitSelection.hxx
  CubeCorrValues.hxx CubeReconState.hxx CubeVertexState.hxx
  CubeClusterState.hxx CubeShowerState.hxx CubeTrackState.hxx
  CubeReconObject.hxx CubeReconNode.hxx
  CubeReconVertex.hxx CubeReconCluster.hxx
  CubeReconShower.hxx CubeReconTrack.hxx
  CubeAlgorithm.hxx CubeAlgorithmResult.hxx
  CubeG4Hit.hxx CubeG4Trajectory.hxx
  TUnitsTable.hxx CubeRecurseGeometry.hxx
  )

# Make sure the current directories are available for the later
# compilation.  This is required to make ROOT_GENERATE_DICTIONARY
# work.
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Build the dictionary for the i/o classes.
ROOT_GENERATE_DICTIONARY(G__cuberecon_io
  CubeHandle.hxx CubeEvent.hxx
  CubeHit.hxx CubeHitSelection.hxx
  CubeCorrValues.hxx CubeReconState.hxx CubeVertexState.hxx
  CubeClusterState.hxx CubeShowerState.hxx CubeTrackState.hxx
  CubeReconObject.hxx CubeReconNode.hxx
  CubeReconVertex.hxx CubeReconCluster.hxx
  CubeReconShower.hxx CubeReconTrack.hxx
  CubeAlgorithm.hxx CubeAlgorithmResult.hxx
  CubeG4Hit.hxx CubeG4Trajectory.hxx
  OPTIONS -inlineInputHeader
  LINKDEF cuberecon_io_LinkDef.h)

# Build the library.
add_library(cuberecon_io SHARED ${source} G__cuberecon_io.cxx)

target_include_directories(cuberecon_io PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:include/CubeRecon>")

target_link_libraries(cuberecon_io PUBLIC ${ROOT_LIBRARIES})

# Install the io library for CubeRecon
install(TARGETS cuberecon_io
  EXPORT CubeReconTargets
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include/CubeRecon )

# Install the header files.
install(FILES ${includes} DESTINATION include/CubeRecon)

# If this is ROOT6 or later, then install the rootmap and pcm files.
if(${ROOT_VERSION} VERSION_GREATER 6)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libcuberecon_io.rootmap
    DESTINATION lib)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libcuberecon_io_rdict.pcm
    DESTINATION lib)
endif(${ROOT_VERSION} VERSION_GREATER 6)
