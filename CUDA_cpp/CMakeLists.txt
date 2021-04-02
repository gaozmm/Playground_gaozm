cmake_minimum_required(VERSION 3.17)
project(CUDA_cpp CUDA)

set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_FLAGS "-Wall")
#set(C:/Users/Gao_z/CLionProjects/Playground_gaozm/CUDA_cpp main.cu)
add_executable(CUDA_cpp main.cu)

set_target_properties(
        CUDA_cpp
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)