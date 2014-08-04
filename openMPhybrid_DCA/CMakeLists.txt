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
cuda_add_executable(go hybridDCAdriver3.cu funct_bin/Assemble.cu funct_bin/Disassemble.cu funct_bin/findCutoff.cu funct_bin/DCAfuncts.cu funct_bin/Functs.cu funct_bin/hybridDCA.cu funct_bin/Initialize.cu funct_bin/npy.cu funct_bin/RK45.cu funct_bin/solveBCs.cu funct_bin/updatezetas.cu h_code/Assemble_setup.cu h_code/Disassemble_setup.cu h_code/Init_setup.cu)

