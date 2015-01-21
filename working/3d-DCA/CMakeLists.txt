# Minimum version of CMake required. Don't touch.
cmake_minimum_required(VERSION 2.8)

# Set the name of your project
project(go)

# Include macros from the SBEL utils library


##########################################################
# If using Makefiles, be sure to set CMAKE_BUILD_TYPE to # 
# select different build configurations.                 #
# Possible values:                                       #
# - RELEASE - fully optimized for speed                  #
# - DEBUG - unoptimized, with debugging symbols          #
# - RELWITHDEBUGINFO - optimized, with debugging symbols #
# - MINSIZEREL - optimized for size                      #
# - [not set] - unoptimized                              #
##########################################################


find_package(CUDA)
#cuda_include_directories(/usr/local/cuda/include/cusp)
cuda_add_executable(go cudaDCAdriver.cu Initialize.cu hybridDCA.cu RK45.cu DCAfuncts.cu Disassemble.cpp solveBCs.cu Assemble.cpp)

