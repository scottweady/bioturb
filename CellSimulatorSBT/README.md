# CellSimulator
Cells simulated as sphero-cylinders with growth, collision, and simple hydrodynamics

## 0. Get the code and its submodule
If you get this code as a git repository, just type
```bash
git checkout master
git pull
git submodule init && git submodule update
```
This will automatically get the submodule `SimToolbox` at the correct commit.

If you get this code as a zip file, go to `https://github.com/wenyan4work/SimToolbox` and download the master branch of `SimToolbox` as a zip file. 
Then, unzip it to the folder `SimToolbox`. 
There should be a file at `SimToolbox/CMakeLists.txt`.

## 1. Setup your toolchain properly
Read the document `SimToolbox/ToolChainSetup.md` to setup your toolchain correctly either on your Mac or on your Linux workstation. It is critical to setup thing exactly as documented otherwise you may runinto weird problems.

## 2. Compile this code
This code now uses the cmake build system to make it easier to compile things. 

First, make a directory for the cmake configuration and compilation files, under the base folder of this directory:
```bash
mkdir build
```

Second, go to that folder and run the script if using gcc on Linux or using clang on Mac:
```bash
cd build
bash ../gcc-do-cmake-release.sh
```

If using Intel compilers, run the script for Intel instead:
```bash
cd build
bash ../intel-do-cmake-release.sh
```

After the script completes, you may or may not have to set the correct paths for dependent libraries manually.
Under the folder `build`, type:
```bash
ccmake ../
```

You will see an interactive cmake interface, showing the paths for those dependence libraries:
```bash
 CMAKE_BUILD_TYPE                 Release                                      
 CMAKE_INSTALL_PREFIX             /usr/local                                   
 Eigen3_DIR                       /home/wyan/local/share/eigen3/cmake          
 PVFMM_INCLUDE_DIR                /home/wyan/local/include/pvfmm               
 PVFMM_LIBRARY                    /home/wyan/local/lib/pvfmm/libpvfmm.a        
 TRNG_INCLUDE_DIR                 /home/wyan/local/include                     
 TRNG_LIBRARY                     /home/wyan/local/lib/libtrng4.a              
 Trilinos_DIR                     /home/wyan/local/lib/cmake/Trilinos          
 yaml-cpp_DIR                     /home/wyan/local/lib/cmake/yaml-cpp       
```
If you installed Eigen3, TRNG and yaml with homebrew, their paths will be under `/usr/local`. If you installed them manually as instructed in `SimToolbox/CMakeLists.txt`, you will find all paths under `/your/home/folder/local`.
If any of these paths show `Not Found` or incorrect, manually type in the correct path. 

Then follow the directions on the lower half of the screen:
```bash
CMAKE_BUILD_TYPE: Choose the type of build, options are: None(CMAKE_CXX_FLAGS or
Press [enter] to edit option Press [d] to delete an entry   CMake Version 3.10.2
Press [c] to configure       Press [g] to generate and exit
Press [h] for help           Press [q] to quit without generating
Press [t] to toggle advanced mode (Currently Off)
``` 

First press `c` to configure. If no error happens, press `g` to generate and exit. 

Then you will exit the cmake interface and come back to the bash terminal line. Still under the directory `build`, simply type
```bash
make -j4
```
If you computer have more cores, type `make -jN` where `N` is the number of your cores to compile things in parallel.

After the compilation, you should see the last line of output as this:
```bash
[100%] Built target CellSimulator.X
```

This is the executable file. Copy this to anywhere you want to run the program, together with the runConfig.txt (required, you must have one) and posInitial.dat (optional, if you have one).