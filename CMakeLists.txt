cmake_minimum_required(VERSION 3.0.0)
project(project_smoke VERSION 0.1.0)

include(CTest)
enable_testing()

add_executable(project_smoke main.cpp)

target_compile_options(project_smoke PUBLIC -O3)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
