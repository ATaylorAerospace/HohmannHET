# Author: A Taylor | Purpose: CMake build configuration | Ref: Vallado/Curtis
cmake_minimum_required(VERSION 3.20)

project(hohmann_het
    VERSION 1.0.0
    DESCRIPTION "Low-thrust orbital transfer library: Keplerian dynamics, HET propulsion, and mission optimization."
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---------------------------------------------------------------------------
# Header-only interface library
# ---------------------------------------------------------------------------
add_library(hohmann_het INTERFACE)
target_include_directories(hohmann_het INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(hohmann_het INTERFACE cxx_std_20)

# Strict warnings on supported compilers
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(hohmann_het INTERFACE -Wall -Wextra -Wpedantic)
elseif(MSVC)
    target_compile_options(hohmann_het INTERFACE /W4)
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
option(HOHMANN_HET_BUILD_TESTS "Build GoogleTest suite" ON)
if(HOHMANN_HET_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# ---------------------------------------------------------------------------
# Install rules
# ---------------------------------------------------------------------------
install(DIRECTORY include/ DESTINATION include)
install(TARGETS hohmann_het EXPORT hohmann_het_targets)
install(EXPORT hohmann_het_targets
    FILE hohmann_het-config.cmake
    DESTINATION lib/cmake/hohmann_het
)
