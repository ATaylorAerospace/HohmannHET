# Author: A Taylor | Purpose: CMake test configuration | Ref: Vallado/Curtis
include(FetchContent)

FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
# Prevent GTest from overriding compiler/linker settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# ---------------------------------------------------------------------------
# Helper macro: add a gtest executable and register it with CTest
# ---------------------------------------------------------------------------
macro(add_hohmann_test target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} PRIVATE hohmann_het GTest::gtest_main GTest::gmock)
    target_compile_features(${target} PRIVATE cxx_std_20)
    add_test(NAME ${target} COMMAND ${target})
endmacro()

add_hohmann_test(test_dynamics     test_dynamics.cpp)
add_hohmann_test(test_propulsion   test_propulsion.cpp)
add_hohmann_test(test_optimization test_optimization.cpp)
