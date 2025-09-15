#! /bin/bash

cmake \
  -D CMAKE_CXX_COMPILER=mpicxx \
  -D CMAKE_C_COMPILER=mpicc \
  -D CMAKE_CXX_FLAGS="-O3 -march=broadwell" \
  -D CMAKE_BUILD_TYPE=Release \
  -D SFTPATH="${HOME}/ceph/envs/stb_intel_13_rocky8/" \
  -D pvfmm_DIR="${HOME}/ceph/envs/stb_intel_13_rocky8/share/pvfmm/" \
../
