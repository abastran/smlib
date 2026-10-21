# SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
# SPDX-License-Identifier: Apache-2.0
#
# Cross-build for Windows (MSVC ABI) from Linux or macOS with clang-cl and the Windows SDK and CRT
# that `cargo xwin` downloads (the same ones the Rust side links against):
#
#   cmake -S capi -B build-win -G Ninja -DCMAKE_BUILD_TYPE=Release \
#     -DCMAKE_TOOLCHAIN_FILE=$PWD/capi/cmake/clang-cl-xwin.cmake -DSMLIB_XWIN_ARCH=x86_64
#
# SMLIB_XWIN_ARCH: x86_64 (default), aarch64 or i686 (xwin downloads i686 only when asked:
# `xwin --arch x86`). XWIN_DIR: the xwin cache (default ~/.cache/cargo-xwin/xwin).

set(CMAKE_SYSTEM_NAME Windows)
set(SMLIB_XWIN_ARCH "x86_64" CACHE STRING "x86_64, aarch64 or i686")
if(NOT DEFINED XWIN_DIR)
  if(DEFINED ENV{XWIN_DIR})
    set(XWIN_DIR "$ENV{XWIN_DIR}")
  else()
    set(XWIN_DIR "$ENV{HOME}/.cache/cargo-xwin/xwin")
  endif()
endif()

if(SMLIB_XWIN_ARCH STREQUAL "aarch64")
  set(CMAKE_SYSTEM_PROCESSOR ARM64)
  set(_lib_arch aarch64)
elseif(SMLIB_XWIN_ARCH STREQUAL "i686")
  set(CMAKE_SYSTEM_PROCESSOR X86)
  set(_lib_arch x86)
else()
  set(CMAKE_SYSTEM_PROCESSOR AMD64)
  set(_lib_arch x86_64)
endif()
set(_target ${SMLIB_XWIN_ARCH}-pc-windows-msvc)

# The LLVM tools; distributions often keep the unversioned llvm-lib/-rc/-mt names out of PATH.
find_program(_clang_cl NAMES clang-cl REQUIRED NO_CMAKE_FIND_ROOT_PATH)
get_filename_component(_llvm_bin "${_clang_cl}" REALPATH)
get_filename_component(_llvm_bin "${_llvm_bin}" DIRECTORY)
find_program(_lld_link NAMES lld-link HINTS ${_llvm_bin} REQUIRED NO_CMAKE_FIND_ROOT_PATH)
find_program(_llvm_lib NAMES llvm-lib HINTS ${_llvm_bin} REQUIRED NO_CMAKE_FIND_ROOT_PATH)
find_program(_llvm_rc NAMES llvm-rc HINTS ${_llvm_bin} NO_CMAKE_FIND_ROOT_PATH)
find_program(_llvm_mt NAMES llvm-mt HINTS ${_llvm_bin} NO_CMAKE_FIND_ROOT_PATH)
set(CMAKE_C_COMPILER ${_clang_cl})
set(CMAKE_CXX_COMPILER ${_clang_cl})
set(CMAKE_LINKER ${_lld_link})
set(CMAKE_AR ${_llvm_lib})
if(_llvm_rc)
  set(CMAKE_RC_COMPILER ${_llvm_rc})
endif()
if(_llvm_mt)
  set(CMAKE_MT ${_llvm_mt})
endif()
set(CMAKE_C_COMPILER_TARGET ${_target})
set(CMAKE_CXX_COMPILER_TARGET ${_target})

set(_inc "/imsvc${XWIN_DIR}/crt/include /imsvc${XWIN_DIR}/sdk/include/ucrt /imsvc${XWIN_DIR}/sdk/include/um /imsvc${XWIN_DIR}/sdk/include/shared")
set(CMAKE_C_FLAGS_INIT "${_inc}")
set(CMAKE_CXX_FLAGS_INIT "${_inc}")
set(_libs "/libpath:${XWIN_DIR}/crt/lib/${_lib_arch} /libpath:${XWIN_DIR}/sdk/lib/um/${_lib_arch} /libpath:${XWIN_DIR}/sdk/lib/ucrt/${_lib_arch}")
set(CMAKE_EXE_LINKER_FLAGS_INIT "${_libs}")
set(CMAKE_SHARED_LINKER_FLAGS_INIT "${_libs}")
set(CMAKE_MODULE_LINKER_FLAGS_INIT "${_libs}")
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
