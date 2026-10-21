# SMLib C API (SolidCraft)

One static library, `smlib`, holding SMLib, SM_API and a C boundary in front of them, with
`include/smlib.h` as its only public header. It's built with CMake alone. Nothing from USD,
Python, Packman or premake is involved, and OpenCASCADE is never part of it.

```sh
cmake -S capi -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build          # the C smoke test
```

A clean build compiles about 235 files, which takes about a minute at `-j8`. The archive holds
everything, so a C or Rust program links `libsmlib.a` (`smlib.lib`) plus the platform's C++
runtime: `libstdc++` on Linux, `libc++` on macOS, and the MSVC runtime as a DLL (`/MD`, as Rust
links it) on Windows.

## The boundary

`include/smlib.h` sets the rules for every function:
- **Only C data crosses:** scalars, caller-owned arrays with their lengths, and opaque handles.
- **Status codes, never exceptions.** Every call returns an `int32_t`:
  - `SMLIB_OK` (0);
  - a positive SMLib status such as `SM_ERR` (1001) or `SM_ERR_INVALID_INPUT` (1007);
  - a negative `SMLIB_ERR_*` from the shim.
- **Exceptions are caught.** Every entry point catches every C++ exception and returns it as
  `SMLIB_ERR_EXCEPTION`.
- **Inputs are never changed.** Where SM_API edits or deletes its operands, the shim works on a
  copy.
- **Handles belong to the caller,** who frees them with `smlib_brep_free`.
- **Calls are serialised.** SMLib keeps process-global state, so one lock serialises every call.

So far the API covers what's needed to prove the build: init, copy, free, box, volume, counts
and native file read/write. SolidCraft's kernel adapter adds the rest, one shim at a time, each
with a case in the fuzz target.

## Diagnostics

- **Log:** `smlib_set_log` receives one line per call, written after it returns. Example:
  `smlib_make_box(origin=[0,0,0], size=[10,20,30]) -> 0 (0.21 ms)`.
- **Replay dumps:** with `smlib_set_dump_dir` set, every failing call leaves its line in
  `<n>-<function>.txt` and each B-rep it was given as a native `<n>-<function>-<input>.smb`. To
  replay a call, read those back with `smlib_brep_read` and call the function again.
- **Crash reports:** `smlib_current_call` copies the call in flight into a caller buffer. It takes
  no lock and doesn't allocate, so a crash handler can record it.

## Fuzzing

`fuzz/fuzz_capi.cpp` is a libFuzzer target that calls every shim with arbitrary arguments. SMLib
itself is built with ASan and UBSan for it:

```sh
cmake -S capi -B build-fuzz -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DSMLIB_CAPI_FUZZ=ON -DSMLIB_CAPI_TESTS=OFF
cmake --build build-fuzz --target smlib_fuzz
build-fuzz/smlib_fuzz -max_total_time=120 corpus/
```

`fuzz/fuzz_history.cpp` (`smlib_fuzz_history`) does the same for the history shims, through a
real boolean of two tagged boxes.

## Targets

| Target | How | Checked by |
|---|---|---|
| Linux x86_64, aarch64 | native gcc/clang | CI (`.github/workflows/capi.yml`) |
| macOS universal | `-DCMAKE_OSX_ARCHITECTURES="arm64;x86_64"` | CI |
| Windows x64, ARM64 (MSVC) | Visual Studio generator, `-A x64` / `-A ARM64` | CI |
| Windows from Linux/macOS (MSVC ABI) | clang-cl + the `cargo xwin` SDK: `-DCMAKE_TOOLCHAIN_FILE=$PWD/capi/cmake/clang-cl-xwin.cmake [-DSMLIB_XWIN_ARCH=aarch64]` | by hand |
| FreeBSD, RISC-V | not yet in CI; SMLib has no architecture-specific code | — |

On Windows the archive is about 1 GB, although a linked program is only a few MB. Under
`_WIN32`, upstream's headers instantiate every `SmTArray<T>` explicitly in every file
(`SM_TARRAY_TEMPLATE_PREDECLARATION`), so each object is about 7 MB. A clean clang-cl build takes
about 10 minutes at `-j8`. Turning those into `extern template` declarations would fix it; that's
not done yet.

## Changes from upstream

This fork changes upstream SMLib (NVIDIA usd-brep, Apache-2.0) as follows:
- **`SM_RAND_SEED`** (`SmConfig.h`, `SmFace.cpp`, `SmShell.cpp`): upstream seeds `rand()` from
  the clock before its random sampling, so results could differ from run to run. This build
  fixes the seed (`-DSMLIB_RAND_SEED=`; leave it empty for upstream's behaviour).
- **`SM_STATIC`** (`SmTypes.h`, `SmApiTypes.h`, `nurbs.h`): built into a static library, nothing
  is marked `__declspec(dllexport)`/`dllimport` on Windows.
