cmake_minimum_required(VERSION 3.14)
project(maf_counter_count LANGUAGES CXX)

# Check if using GNU and if the version is less than 13.2.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "13.2")
  message(WARNING "Your current GNU compiler version (${CMAKE_CXX_COMPILER_VERSION}) is less than 13.2 which is required for proper std::filesystem support.")
  find_program(CLANGPP clang++)
  if(CLANGPP)
    message(STATUS "Found clang++ at ${CLANGPP}. Consider using it with -DCMAKE_CXX_COMPILER=${CLANGPP}")
  else()
    message(FATAL_ERROR "No suitable compiler found. Please upgrade your GNU compiler to 13.2 or later, or install clang++ and re-run CMake with -DCMAKE_CXX_COMPILER=clang++.")
  endif()
endif()

# Set C++ standard and enforce C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Force -O3 optimizations in Release mode
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Binaries will be placed in build/bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Add include directories:
# - The project's main include directory
# - Boost.Multiprecision submodule include directory
# - Boost.Config submodule include directory (required by Multiprecision)
include_directories(
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/include/boost_multiprecision/include
  ${CMAKE_SOURCE_DIR}/include/boost_config/include
)

# Find required system libraries
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

# ----------------------------------------------------------------------------
# 1) maf_counter_count
# ----------------------------------------------------------------------------
add_executable(maf_counter_count src/maf_counter_count.cpp)
target_link_libraries(maf_counter_count PRIVATE Threads::Threads ZLIB::ZLIB)

# ----------------------------------------------------------------------------
# 2) maf_counter_dump
# ----------------------------------------------------------------------------
add_executable(maf_counter_dump src/maf_counter_dump.cpp)
target_link_libraries(maf_counter_dump PRIVATE Threads::Threads)

# ----------------------------------------------------------------------------
# 3) maf_counter_tools
# ----------------------------------------------------------------------------
add_executable(maf_counter_tools src/maf_counter_tools.cpp)
# If maf_counter_tools needs zlib as well, link it (uncomment next line if needed):
target_link_libraries(maf_counter_tools PRIVATE Threads::Threads ZLIB::ZLIB)

# ----------------------------------------------------------------------------
# 4) maf_counter_proteomes_count
# ----------------------------------------------------------------------------
add_executable(maf_counter_proteomes_count src/maf_counter_proteomes_count.cpp)
# Assuming it requires zlib as well (like maf_counter_count)
target_link_libraries(maf_counter_proteomes_count PRIVATE Threads::Threads ZLIB::ZLIB)

# ----------------------------------------------------------------------------
# 5) maf_counter_proteomes_dump
# ----------------------------------------------------------------------------
add_executable(maf_counter_proteomes_dump src/maf_counter_proteomes_dump.cpp)
# Assuming it does NOT use zlib (like maf_counter_dump)
target_link_libraries(maf_counter_proteomes_dump PRIVATE Threads::Threads)
