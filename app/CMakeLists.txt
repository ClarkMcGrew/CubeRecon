# Build all of the applications here.

# Configure the dependencies
find_package(ROOT REQUIRED
  COMPONENTS Geom EG Physics Matrix MathCore Tree RIO
  Eve TreePlayer Ged RGL Gui Graf Graf3d)
if(ROOT_FOUND)
  include(${ROOT_USE_FILE})
endif(ROOT_FOUND)

# Add a program to translate an ERepSim event into simple hits.
add_executable(cubeERepTranslate cubeERepTranslate.cxx)
target_link_libraries(cubeERepTranslate LINK_PUBLIC cuberecon)
install(TARGETS cubeERepTranslate RUNTIME DESTINATION bin)

# Add a program to generate voxel hits.
add_executable(cubeMakeHits3D cubeMakeHits3D.cxx)
target_link_libraries(cubeMakeHits3D LINK_PUBLIC cuberecon)
install(TARGETS cubeMakeHits3D RUNTIME DESTINATION bin)

# Add a program to run the reconstruction.
add_executable(cubeRecon cubeRecon.cxx)
target_link_libraries(cubeRecon LINK_PUBLIC cuberecon)
install(TARGETS cubeRecon RUNTIME DESTINATION bin)

# Add a program to run the reconstruction.
add_executable(cubeStrip cubeStrip.cxx)
target_link_libraries(cubeStrip LINK_PUBLIC cuberecon)
install(TARGETS cubeStrip RUNTIME DESTINATION bin)

# Add a test program
add_executable(testTrackTiming.exe testTrackTiming.cxx)
target_link_libraries(testTrackTiming.exe LINK_PUBLIC
  cuberecon_io cuberecon_tools)
install(TARGETS testTrackTiming.exe RUNTIME DESTINATION bin)

# Add a test program
add_executable(testHitTiming.exe testHitTiming.cxx)
target_link_libraries(testHitTiming.exe LINK_PUBLIC
  cuberecon_io cuberecon_tools)
install(TARGETS testHitTiming.exe RUNTIME DESTINATION bin)

# Add a test program
add_executable(testNodeCharge.exe testNodeCharge.cxx)
target_link_libraries(testNodeCharge.exe LINK_PUBLIC
  cuberecon_io cuberecon_tools)
install(TARGETS testNodeCharge.exe RUNTIME DESTINATION bin)

# Add a test program
add_executable(testTrackEff.exe testTrackEff.cxx)
target_link_libraries(testTrackEff.exe LINK_PUBLIC
  cuberecon_io cuberecon_tools)
install(TARGETS testTrackEff.exe RUNTIME DESTINATION bin)

# Add a test program
add_executable(testCC0PiTimes.exe testCC0PiTimes.cxx)
target_link_libraries(testCC0PiTimes.exe LINK_PUBLIC
  cuberecon_io cuberecon_tools ${ROOT_LIBRARIES} )
install(TARGETS testCC0PiTimes.exe RUNTIME DESTINATION bin)

# Add a test program
add_executable(testVertexResolution.exe testVertexResolution.cxx)
target_link_libraries(testVertexResolution.exe LINK_PUBLIC
  cuberecon_io cuberecon_tools)
install(TARGETS testVertexResolution.exe RUNTIME DESTINATION bin)

# Add a test program
add_executable(testRawHits.exe testRawHits.cxx)
target_link_libraries(testRawHits.exe LINK_PUBLIC
  cuberecon_io cuberecon_tools)
install(TARGETS testRawHits.exe RUNTIME DESTINATION bin)

# Add a test program
add_executable(testECalEff.exe testECalEff.cxx)
target_link_libraries(testECalEff.exe LINK_PUBLIC
  cuberecon_io cuberecon_tools)
install(TARGETS testECalEff.exe RUNTIME DESTINATION bin)

# Add a test program
add_executable(testCCInc.exe testCCInc.cxx)
target_link_libraries(testCCInc.exe LINK_PUBLIC
  cuberecon_io cuberecon_tools)
install(TARGETS testCCInc.exe RUNTIME DESTINATION bin)
