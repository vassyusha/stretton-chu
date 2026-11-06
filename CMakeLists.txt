cmake_minimum_required(VERSION 3.10)
project(stretton_chu LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_library(stretton_chu_lib src/Solver.cpp)
target_include_directories(stretton_chu_lib PUBLIC src)
target_compile_options(stretton_chu_lib PUBLIC -Wall -Wextra)

add_executable(stretton_chu src/main.cpp)
target_link_libraries(stretton_chu PRIVATE stretton_chu_lib)
