cmake_minimum_required(VERSION 3.13)
project(1D_Diffusion)

#   -- Source Code Path --
set(DIFFUSION_CODE_PATH ${CMAKE_SOURCE_DIR}/code)



#   == SOURCE FILES ==
file(GLOB_RECURSE prms ${DIFFUSION_CODE_PATH}/prms/Parameters.h)
file(GLOB_RECURSE clss ${DIFFUSION_CODE_PATH}/clss/*.*)
file(GLOB_RECURSE hdrs ${DIFFUSION_CODE_PATH}/hdrs/*.h*)
file(GLOB_RECURSE main ${DIFFUSION_CODE_PATH}/main.cpp)



#   == COMPILER FLAGS ==
#   -- Set Standard --
set(CMAKE_CXX_STANDARD 17)



#   == EXECUTABLE CREATION ==
set(SOURCE_FILES ${prms} ${clss} ${hdrs} ${main})
add_executable(1D_Diffusion ${SOURCE_FILES})



#   == INCLUDE HEADER FILES ==
target_include_directories(1D_Diffusion PUBLIC ${DIFFUSION_CODE_PATH}/prms)
target_include_directories(1D_Diffusion PUBLIC ${DIFFUSION_CODE_PATH}/clss)
target_include_directories(1D_Diffusion PUBLIC ${DIFFUSION_CODE_PATH}/hdrs)