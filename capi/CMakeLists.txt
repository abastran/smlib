# SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
# SPDX-License-Identifier: Apache-2.0
#
# SMLib + SM_API + the C API as ONE static library, `smlib`, with include/smlib.h as its only
# public header. Nothing from USD, Python, Packman or premake is needed:
#
#   cmake -S capi -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build --parallel
#   ctest --test-dir build
#
# See capi/README.md for each platform and for cross builds.

cmake_minimum_required(VERSION 3.20)
project(smlib_capi C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# Rust links the MSVC runtime as a DLL (/MD); so must we.
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")

option(SMLIB_CAPI_TESTS "Build the C API smoke test" ON)
option(SMLIB_CAPI_FUZZ "Build the libFuzzer target (clang only)" OFF)
set(SMLIB_RAND_SEED "20240101" CACHE STRING "Seed for SMLib's random sampling (empty: the clock, as upstream)")

set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/../source)

if(SMLIB_CAPI_FUZZ)
  # SMLib itself is instrumented too, so the fuzzer sees what happens inside it.
  add_compile_options(-fsanitize=fuzzer-no-link,address,undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address,undefined)
endif()

file(GLOB SMLIB_SOURCES CONFIGURE_DEPENDS ${SRC}/SMLib/src/*.cpp)
# Precompiled-header stubs and the DLL entry point have no place in a static library.
list(FILTER SMLIB_SOURCES EXCLUDE REGEX "/(StdAfx|SmDllMain)\\.cpp$")
file(GLOB SMAPI_SOURCES CONFIGURE_DEPENDS ${SRC}/SM_API/src/*.cpp)
list(FILTER SMAPI_SOURCES EXCLUDE REGEX "/StdAfx\\.cpp$")

set(SMLIB_DEFINES NDEBUG NOMINMAX SM_STATIC)
if(NOT SMLIB_RAND_SEED STREQUAL "")
  list(APPEND SMLIB_DEFINES SM_RAND_SEED=${SMLIB_RAND_SEED}u)
endif()

if(MSVC)
  set(SMLIB_WARNINGS /w /bigobj /EHsc)
else()
  set(SMLIB_WARNINGS -w)
endif()

# SMLib/inc and SM_API/inc each have a StdAfx.h; each part must find its own first.
add_library(smlib_core OBJECT ${SMLIB_SOURCES})
target_include_directories(smlib_core PRIVATE ${SRC}/SMLib/inc)
target_compile_definitions(smlib_core PRIVATE ${SMLIB_DEFINES})
target_compile_options(smlib_core PRIVATE ${SMLIB_WARNINGS})

add_library(smlib_api OBJECT ${SMAPI_SOURCES})
target_include_directories(smlib_api PRIVATE ${SRC}/SM_API/inc ${SRC}/SMLib/inc)
target_compile_definitions(smlib_api PRIVATE ${SMLIB_DEFINES})
target_compile_options(smlib_api PRIVATE ${SMLIB_WARNINGS})

# The shims: every capi/src/*.cpp, one file per area, all sharing src/smlib_shim.h.
file(GLOB CAPI_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_library(smlib STATIC ${CAPI_SOURCES} $<TARGET_OBJECTS:smlib_core> $<TARGET_OBJECTS:smlib_api>)
target_include_directories(smlib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_include_directories(smlib SYSTEM PRIVATE ${SRC}/SM_API/inc ${SRC}/SMLib/inc)
target_compile_definitions(smlib PRIVATE ${SMLIB_DEFINES})
if(MSVC)
  target_compile_options(smlib PRIVATE /W3 /EHsc)
else()
  target_compile_options(smlib PRIVATE -Wall -Wextra)
endif()
# A C program linking the archive needs the C++ runtime; tell CMake consumers so.
set_target_properties(smlib PROPERTIES LINKER_LANGUAGE CXX)

install(TARGETS smlib ARCHIVE DESTINATION lib)
install(FILES include/smlib.h DESTINATION include)

if(SMLIB_CAPI_TESTS)
  enable_testing()
  add_executable(smlib_smoke tests/smoke.c)
  target_link_libraries(smlib_smoke PRIVATE smlib)
  set_target_properties(smlib_smoke PROPERTIES LINKER_LANGUAGE CXX)
  if(NOT MSVC)
    target_link_libraries(smlib_smoke PRIVATE m)
  endif()
  add_test(NAME smlib_smoke COMMAND smlib_smoke)
  # The operations from plain C (tests/ops_kernel.cpp, built below, covers them with SM_API).
  add_executable(smlib_ops_test tests/ops.c)
  target_link_libraries(smlib_ops_test PRIVATE smlib)
  set_target_properties(smlib_ops_test PROPERTIES LINKER_LANGUAGE CXX)
  if(NOT MSVC)
    target_link_libraries(smlib_ops_test PRIVATE m)
  endif()
  add_test(NAME smlib_ops COMMAND smlib_ops_test)
  set(SMOKE_DUMP_DIR ${CMAKE_CURRENT_BINARY_DIR}/smoke-dumps)
  file(MAKE_DIRECTORY ${SMOKE_DUMP_DIR})
  set_tests_properties(smlib_smoke PROPERTIES ENVIRONMENT "SMLIB_SMOKE_DUMP_DIR=${SMOKE_DUMP_DIR}")
  # History through booleans of primitives; runs SM_API itself, so it sees SMLib's headers.
  add_executable(smlib_history_test tests/history.cpp)
  target_include_directories(smlib_history_test SYSTEM PRIVATE src ${SRC}/SM_API/inc ${SRC}/SMLib/inc)
  target_compile_definitions(smlib_history_test PRIVATE ${SMLIB_DEFINES})
  target_link_libraries(smlib_history_test PRIVATE smlib)
  add_test(NAME smlib_history COMMAND smlib_history_test)
  # One C++ test per other shim file (tests/<area>.cpp), free to build its inputs with SM_API.
  file(GLOB CAPI_TESTS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
  list(FILTER CAPI_TESTS EXCLUDE REGEX "/(history|flat_roundtrip)\\.cpp$")
  foreach(test ${CAPI_TESTS})
    get_filename_component(name ${test} NAME_WE)
    add_executable(smlib_test_${name} ${test})
    target_link_libraries(smlib_test_${name} PRIVATE smlib)
    target_include_directories(smlib_test_${name} SYSTEM PRIVATE ${SRC}/SM_API/inc ${SRC}/SMLib/inc)
    target_compile_definitions(smlib_test_${name} PRIVATE ${SMLIB_DEFINES} SMLIB_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
    add_test(NAME smlib_test_${name} COMMAND smlib_test_${name})
  endforeach()

  # Flat arrays: SmBrep -> arrays -> SmBrep round trips, and hostile arrays. The test builds its
  # solids with SM_API directly and hands them to the C API as native files.
  add_executable(smlib_flat_roundtrip tests/flat_roundtrip.cpp)
  target_link_libraries(smlib_flat_roundtrip PRIVATE smlib)
  target_include_directories(smlib_flat_roundtrip SYSTEM PRIVATE ${SRC}/SM_API/inc ${SRC}/SMLib/inc)
  target_compile_definitions(smlib_flat_roundtrip PRIVATE ${SMLIB_DEFINES})
  add_test(NAME smlib_flat_roundtrip COMMAND smlib_flat_roundtrip ${CMAKE_CURRENT_BINARY_DIR}/flat-roundtrip)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/flat-roundtrip)
endif()

if(SMLIB_CAPI_FUZZ)
  add_executable(smlib_fuzz fuzz/fuzz_capi.cpp)
  target_link_libraries(smlib_fuzz PRIVATE smlib)
  target_compile_options(smlib_fuzz PRIVATE -fsanitize=fuzzer,address,undefined)
  target_link_options(smlib_fuzz PRIVATE -fsanitize=fuzzer,address,undefined)
  # The history shims, through a real boolean (run through SM_API, so it sees SMLib's headers).
  add_executable(smlib_fuzz_history fuzz/fuzz_history.cpp)
  target_include_directories(smlib_fuzz_history SYSTEM PRIVATE src ${SRC}/SM_API/inc ${SRC}/SMLib/inc)
  target_compile_definitions(smlib_fuzz_history PRIVATE ${SMLIB_DEFINES})
  target_link_libraries(smlib_fuzz_history PRIVATE smlib)
  target_compile_options(smlib_fuzz_history PRIVATE -fsanitize=fuzzer,address,undefined)
  target_link_options(smlib_fuzz_history PRIVATE -fsanitize=fuzzer,address,undefined)
  # One more target per shim file that needs bodies the C API can't make: fuzz/fuzz_<area>.cpp.
  file(GLOB CAPI_FUZZ CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/fuzz/fuzz_*.cpp)
  list(FILTER CAPI_FUZZ EXCLUDE REGEX "/fuzz_(capi|history)\\.cpp$")
  foreach(target ${CAPI_FUZZ})
    get_filename_component(name ${target} NAME_WE)
    string(REPLACE "fuzz_" "smlib_fuzz_" name ${name})
    add_executable(${name} ${target})
    target_link_libraries(${name} PRIVATE smlib)
    target_include_directories(${name} SYSTEM PRIVATE ${SRC}/SM_API/inc ${SRC}/SMLib/inc)
    target_compile_definitions(${name} PRIVATE ${SMLIB_DEFINES})
    target_compile_options(${name} PRIVATE -fsanitize=fuzzer,address,undefined)
    target_link_options(${name} PRIVATE -fsanitize=fuzzer,address,undefined)
  endforeach()
endif()
