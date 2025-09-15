#!/bin/bash

SOURCE_PATH=../pvfmm-af175f382320e6de8c73afa8133b9ec32a433145

# You can invoke this shell script with additional command-line
# arguments.  They will be passed directly to CMake.
#
EXTRA_ARGS=$@

#
# Each invocation of CMake caches the values of build options in a
# CMakeCache.txt file.  If you run CMake again without deleting the
# CMakeCache.txt file, CMake won't notice any build options that have
# changed, because it found their original values in the cache file.
# Deleting the CMakeCache.txt file before invoking CMake will insure
# that CMake learns about any build options you may have changed.
# Experience will teach you when you may omit this step.
#
rm -f CMakeCache.txt

#
# Enable all primary stable Trilinos packages.
#
cmake \
  -D CMAKE_INSTALL_PREFIX:FILEPATH="$SFTPATH" \
  -D CMAKE_INSTALL_LIBDIR=lib \
  -D CMAKE_BUILD_TYPE:STRING="Release" \
  -D CMAKE_CXX_COMPILER:STRING="mpicxx" \
  -D CMAKE_CXX_FLAGS:STRING="$CXXFLAGS" \
  -D PVFMM_EXTENDED_BC:BOOL=ON \
  -D MKL_INCLUDE_DIR:FILEPATH="$MKL_INCLUDE_DIRS" \
  -D MKL_FFTW_INCLUDE_DIR:FILEPATH="$MKL_INCLUDE_DIRS/fftw" \
  -D MKL_SDL_LIBRARY:STRING="$MKL_LIB_DIRS/libmkl_rt.so" \
  $EXTRA_ARGS \
  $SOURCE_PATH
