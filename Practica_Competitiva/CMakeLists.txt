cmake_minimum_required(VERSION 3.0)
project(Practica_Competitiva)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(Practica_Competitiva main.cpp)

if(OpenMP_CXX_FOUND)
    target_link_libraries(Practica_Competitiva PUBLIC OpenMP::OpenMP_CXX)
endif()
