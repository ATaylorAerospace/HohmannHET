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

# Strict warnings for this project's own translation units (the test suite).
# These are deliberately NOT attached to the INTERFACE library: doing so would
# force -Wall/-Wextra/-Wpedantic (or /W4) onto every downstream target that
# links hohmann_het via add_subdirectory(), overriding the consumer's own
# warning policy.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(HOHMANN_HET_WARNING_FLAGS -Wall -Wextra -Wpedantic)
elseif(MSVC)
    set(HOHMANN_HET_WARNING_FLAGS /W4)
else()
    set(HOHMANN_HET_WARNING_FLAGS "")
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
# Build the GoogleTest suite by default only when this is the top-level
# project. When consumed via add_subdirectory() the default is OFF so that a
# downstream build does not fetch/compile GoogleTest and three test binaries
# it never asked for. Consumers can still opt in with
# -DHOHMANN_HET_BUILD_TESTS=ON.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(HOHMANN_HET_BUILD_TESTS_DEFAULT ON)
else()
    set(HOHMANN_HET_BUILD_TESTS_DEFAULT OFF)
endif()
option(HOHMANN_HET_BUILD_TESTS "Build GoogleTest suite" ${HOHMANN_HET_BUILD_TESTS_DEFAULT})
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
