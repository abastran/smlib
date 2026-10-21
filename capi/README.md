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

So far the API covers what's needed to prove the build (init, copy, free, box, volume, counts,
manifold check and native file read/write), SMLib's own validity check (`smlib_brep_check`),
exact mass properties (`smlib_brep_mass_properties`, `smlib_brep_area`), history and the flat
arrays below. SolidCraft's kernel adapter adds the rest, one shim at a time, each with a case in
a fuzz target.

## Flat B-rep arrays

`smlib_brep_export` turns a B-rep into a `smlib_flat`: plain `int64`/`double` arrays that a
caller reads and writes by id (`smlib_flat_get_i64`/`_f64`, `smlib_flat_set_i64`/`_f64`).
`smlib_brep_import` builds a B-rep from them. SolidCraft maps these arrays onto its own B-rep.
The arrays mirror one member of USD's BrepArray schema (`source/BREP_USD_DATA/inc/UsdBrepArrayData.h`)
with the same layout and units: regions, shells, faceuses, faces, loops, edgeuses with radial
lists, edges, wire edges and vertices; analytic surfaces (plane, cylinder, cone, sphere, torus)
and NURBS; and analytic angles and face ranges in radians. Tokens become small integers and
3-vectors are flattened. `include/smlib.h` documents every array: its meaning, element type,
units and length rule. The conversion is `src/smlib_flat.cpp`, a USD-free port of
`source/BREP_SM_USD/src/SmuConvert.cpp`, and it keeps that converter's behaviour:
- UV trim curves are carried only for faces on NURBS surfaces. On analytic faces, SMLib's trims
  live in its internal NURBS form of the surface rather than the analytic parameterisation, so
  they are not exported and SMLib rebuilds them from the 3D edges on import.
- Edge curves are always exported as NURBS: lines, circles and ellipses as their exact NURBS form,
  and any other curve approximated to the B-rep's tolerance. Import also accepts line, circle
  and ellipse curves.
- A cone with a negative semi-angle is rebuilt in the equivalent positive-angle form on import.
  Its second export therefore differs from its first, and from then on the arrays are stable.
- Attributes, materials and colours are not carried.

Import validates everything first (counts, indices, lengths, kinds, finite values, NURBS
orders and knots, closed loops, radial cycles), so malformed arrays return
`SM_ERR_INVALID_INPUT` without building anything. A failing import's dump holds its arrays as
text, which `smlib_flat_read` reads back for replay. `tests/flat_roundtrip.cpp` round-trips a box,
cylinder, cone, sphere, torus, two filleted boxes and a boolean, then feeds the import hostile
and randomly mutated arrays.

## Diagnostics

- **Log:** `smlib_set_log` receives one line per call, written after it returns. Example:
  `smlib_make_box(origin=[0,0,0], size=[10,20,30]) -> 0 (0.21 ms)`.
- **Replay dumps:** with `smlib_set_dump_dir` set, every failing call leaves its line in
  `<n>-<function>.txt` and each B-rep it was given as a native `<n>-<function>-<input>.smb`. To
  replay a call, read those back with `smlib_brep_read` and call the function again.
- **Crash reports:** `smlib_current_call` copies the call in flight into a caller buffer. It takes
  no lock and doesn't allocate, so a crash handler can record it.

## Fuzzing

`fuzz/fuzz_capi.cpp` is a libFuzzer target that calls every shim with arbitrary arguments, and
each `fuzz/fuzz_<area>.cpp` is one more (`smlib_fuzz_<area>`) for shims that need bodies the C
API can't make yet, broken ones included. SMLib itself is built with ASan and UBSan for them:

```sh
cmake -S capi -B build-fuzz -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DSMLIB_CAPI_FUZZ=ON -DSMLIB_CAPI_TESTS=OFF
cmake --build build-fuzz --target smlib_fuzz
build-fuzz/smlib_fuzz -max_total_time=120 corpus/
```

`fuzz/fuzz_history.cpp` (`smlib_fuzz_history`) does the same for the history shims, through a
real boolean of two tagged boxes.

`SMLIB_FUZZ_SMB=<file.smb>` adds a second starting body for the flat-array case. Point it at a
body with NURBS faces, such as `blend.smb`, which the round-trip test writes to its scratch
directory; that way the NURBS and UV-trim paths get fuzzed too, not only a box.

## Targets

| Target | How | Checked by |
|---|---|---|
| Linux x86_64, aarch64 | native gcc/clang | CI (`.github/workflows/capi.yml`) |
| macOS universal | `-DCMAKE_OSX_ARCHITECTURES="arm64;x86_64"` | CI |
| Windows x64, ARM64 (MSVC) | Visual Studio generator, `-A x64` / `-A ARM64` | CI |
| Windows from Linux/macOS (MSVC ABI) | clang-cl + the `cargo xwin` SDK: `-DCMAKE_TOOLCHAIN_FILE=$PWD/capi/cmake/clang-cl-xwin.cmake [-DSMLIB_XWIN_ARCH=aarch64]` | by hand |
| FreeBSD, RISC-V | not yet in CI; SMLib has no architecture-specific code | — |

On Windows the archive is about 150 MB and a clean clang-cl build takes a few minutes, now that
a static build no longer instantiates every `SmTArray<T>` in every file (see below).

## Changes from upstream

This fork changes upstream SMLib (NVIDIA usd-brep, Apache-2.0) as follows:
- **`SM_RAND_SEED`** (`SmConfig.h`, `SmFace.cpp`, `SmShell.cpp`): upstream seeds `rand()` from
  the clock before its random sampling, so results could differ from run to run. This build
  fixes the seed (`-DSMLIB_RAND_SEED=`; leave it empty for upstream's behaviour).
- **`SM_STATIC`** (`SmTypes.h`, `SmApiTypes.h`, `nurbs.h`): built into a static library, nothing
  is marked `__declspec(dllexport)`/`dllimport` on Windows.
- **Seam-curve leak** (`SmFace.cpp`, `SmFace::HasSeamProblem`): it copied its seam
  classification curve and then disowned the original, leaking a few KB on every call for each
  periodic face (cylinders, cones, spheres, tori). AssertValid runs it on every such face. The
  original is now handed over instead of copied.
- **`SM_TARRAY_TEMPLATE_PREDECLARATION`** (`SmTArray.h`): under `SM_STATIC` it is empty on Windows
  too, as on every other platform, instead of instantiating each listed `SmTArray<T>` in every
  file. The Windows archive drops from about 1 GB to 146 MB; linked programs are byte-for-byte
  the same size. It also lets MinGW compile `SmHashTable.h`.
- **Apex leak** (`SmPrimitiveCreation.cpp`, `CreateCone`): a cone with a zero radius at one end
  built that end's circle and, having no cap to give it to, leaked it.
