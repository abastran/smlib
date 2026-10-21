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

add_library(smlib STATIC src/smlib_c.cpp $<TARGET_OBJECTS:smlib_core> $<TARGET_OBJECTS:smlib_api>)
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
  set(SMOKE_DUMP_DIR ${CMAKE_CURRENT_BINARY_DIR}/smoke-dumps)
  file(MAKE_DIRECTORY ${SMOKE_DUMP_DIR})
  set_tests_properties(smlib_smoke PROPERTIES ENVIRONMENT "SMLIB_SMOKE_DUMP_DIR=${SMOKE_DUMP_DIR}")
endif()

if(SMLIB_CAPI_FUZZ)
  add_executable(smlib_fuzz fuzz/fuzz_capi.cpp)
  target_link_libraries(smlib_fuzz PRIVATE smlib)
  target_compile_options(smlib_fuzz PRIVATE -fsanitize=fuzzer,address,undefined)
  target_link_options(smlib_fuzz PRIVATE -fsanitize=fuzzer,address,undefined)
endif()
