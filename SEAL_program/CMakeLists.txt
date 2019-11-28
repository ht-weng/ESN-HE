cmake_minimum_required(VERSION 3.12)

project(MAIN VERSION 3.4.4 LANGUAGES CXX)


add_executable(main main.cpp)

# Executable will be in ../bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${MAIN_SOURCE_DIR}/)

find_package(SEAL 3.4.4 EXACT REQUIRED)

# Link Microsoft SEAL
target_link_libraries(main SEAL::seal)