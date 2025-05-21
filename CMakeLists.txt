cmake_minimum_required(VERSION 3.15...3.30)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

# Create executable file, not library
#Python_add_library(TMscore MODULE TMscore_py/TMscore.cpp WITH_SOABI)
set(tmscore_src "TMscore_py/TMscore.cpp")
add_executable(TMscore ${tmscore_src})

install(TARGETS TMscore DESTINATION ${SKBUILD_PROJECT_NAME})