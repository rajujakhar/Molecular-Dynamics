

cmake_minimum_required(VERSION 2.8.9)

# Defining the project name
project (SiWiR2_MD)

set(CMAKE_CXX_FLAGS "-Wall -Winline -Wshadow -pedantic")

add_definitions(-std=c++11)

# Bringing in the include directories
include_directories(include)


#Adding the sources using the set command
set(SOURCES test/main.cpp src/ParameterReader.cpp src/ProblemFormulate.cpp src/Particle.cpp)

#Setting the executable file
add_executable(mdsim ${SOURCES})



