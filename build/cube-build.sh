#!/bin/bash
if [ ${#CUBE_ROOT} == 0 ]; then
    echo CUBE is not setup yet.
    echo You can also run cmake and make by hand...
    exit 1
fi


BUILD_LOCATION=${CUBE_ROOT}/${CUBE_TARGET}
if [ ! -d ${BUILD_LOCATION} ]; then
    mkdir -p ${BUILD_LOCATION}
fi

if [ ! -d ${BUILD_LOCATION} ]; then
    echo Unable to access build location at ${BUILD_LOCATION}
    exit 1
fi

cd ${BUILD_LOCATION}

if [ ${#1} != 0 -a "x${1}" == "xforce" ]; then
    shift
    echo Reconfigure build.
    if [ -f  CMakeCache.txt ]; then
	rm CMakeCache.txt
    fi
    if [ -d CMakeFiles ]; then
	rm -rf CMakeFiles
    fi
fi

if [ ! -f CMakeCache.txt ]; then
    cmake -DCMAKE_INSTALL_PREFIX=${CUBE_ROOT}/${CUBE_TARGET} ${CUBE_ROOT}
fi

VERBOSE=true make $* || exit 1
make install || exit 1

echo Use "cube-build force" to force a clean build.
