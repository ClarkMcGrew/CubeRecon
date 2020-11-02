# Build all of the applications here.

# Add a program to test writing the objets.
add_executable(cubeWriteTest cubeWriteTest.cxx)
target_link_libraries(cubeWriteTest LINK_PUBLIC cuberecon)
install(TARGETS cubeWriteTest RUNTIME DESTINATION bin)

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
