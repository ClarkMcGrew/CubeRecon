cmake_minimum_required(VERSION 3.0)

set(CMAKE_BUILD_TYPE Debug)

# Set the project name and version
project(cube-disp)
set(VERSION_MAJOR "1")
set(VERSION_MINOR "0")
set(VERSION_PATCH "0")
set(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

message("Cube Reconstruction Display -- ${VERSION}")

# Set the default verbosity of the compilation...
set(CMAKE_VERBOSE_MAKEFILE TRUE)

# Configure the dependencies
find_package(ROOT REQUIRED
  COMPONENTS Geom Physics Matrix MathCore Tree RIO
  Eve TreePlayer Ged RGL Gui Graf Graf3d)
if(ROOT_FOUND)
  include(${ROOT_USE_FILE})
endif(ROOT_FOUND)

if(ROOT_Eve_LIBRARY)
  message("Building the cube-disp event display")
else(ROOT_Eve_LIBRARY)
  message("NOT BUILDING THE EVENT DISPLAY: eve not found")
  return()
endif(ROOT_Eve_LIBRARY)

# Define the source and include files that should be used for the display.
set(source
  TEventDisplay.cxx TEventChangeManager.cxx TGUIManager.cxx
  TShowHits.cxx THitChangeHandler.cxx
  TG4HitChangeHandler.cxx TTrajectoryChangeHandler.cxx
  TFitChangeHandler.cxx TFindResultsHandler.cxx
  TMatrixElement.cxx TReconClusterElement.cxx TReconTrackElement.cxx
  TReconShowerElement.cxx
  )

set(includes
  TEventDisplay.hxx TEventChangeManager.hxx TGUIManager.hxx
  TRecurseGeoManager.hxx
  TShowHits.hxx THitChangeHandler.hxx
  TG4HitChangeHandler.hxx TFitChangeHandler.hxx TFindResultsHandler.hxx
  TMatrixElement.hxx  TReconClusterElement.hxx TReconTrackElement.hxx
  TReconShowerElement.hxx
  )

# Make sure the current directories are available for the root
# dictionary generation.
include_directories(${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../io
  ${CMAKE_CURRENT_SOURCE_DIR}/../src)

# Build the dictionary for the i/o classes.
ROOT_GENERATE_DICTIONARY(G__cubedisp
  TEventChangeManager.hxx
  OPTIONS -inlineInputHeader
  LINKDEF cubedisp_LinkDef.h)

# Build the library.
add_library(cubedisp SHARED
  ${source}
  G__cubedisp.cxx)
target_include_directories(cubedisp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(cubedisp PUBLIC cuberecon_io ${ROOT_LIBRARIES})
install(TARGETS cubedisp LIBRARY DESTINATION lib)

# Build the display
add_executable(cube-disp eventDisplay.cxx)
target_link_libraries(cube-disp LINK_PUBLIC cubedisp)
install(TARGETS cube-disp RUNTIME DESTINATION bin)

# If this is ROOT6 or later, then install the rootmap and pcm files.
if(${ROOT_VERSION} VERSION_GREATER 6)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libcubedisp.rootmap
    DESTINATION lib)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libcubedisp_rdict.pcm
    DESTINATION lib)
endif(${ROOT_VERSION} VERSION_GREATER 6)
