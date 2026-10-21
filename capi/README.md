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
history, the flat arrays below and the operations. SolidCraft's kernel adapter adds the rest,
one shim at a time, each with a case in a fuzz target.

## Operations

`src/smlib_ops.cpp`: `smlib_boolean` (union, subtract, intersect), `smlib_fillet`,
`smlib_chamfer`, `smlib_shell`, `smlib_offset_faces` and `smlib_thicken`. Each copies its inputs and works on the copies, so tagged
operands (`smlib_history_tag`) give a result whose history `smlib_history_read` reads. An empty
boolean result is `SMLIB_ERR_EMPTY`; several solids come back as one body
(`smlib_brep_solid_count`). A solid input that gives anything but a closed manifold solid is
`SMLIB_ERR_NOT_SOLID`.

Edges and faces are chosen by points, never by SMLib's internal order: each point takes the
nearest edge (face), which must lie within the caller's tolerance and be nearer than any other
by more than the body's modelling tolerance (`SMLIB_ERR_PICK_NONE`, `SMLIB_ERR_PICK_AMBIGUOUS`,
with the point's index). Fillets and chamfers take each picked edge's whole tangent chain
(edges meeting within 0.5°). A chamfer's distance is the setback along each face at the edge's
midpoint; it is cut by SMLib's constant-radius solver with a linear section, whose ball radius
is the setback times tan(a/2), a being the angle between the faces there. A shell's thickness
is negative inward, positive outward. An inward wall that leaves no room inside the body is
`SMLIB_ERR_NO_ROOM`, refused before SMLib's offset sees it (which corrupts memory there,
tangency-labs/solidcore#380): the shim looks for a point inside the solid farther than the wall
(and a margin) from every face but the open ones, certified as the centre of a ball inside the
body that it reached from a face sample along the inward normal, then from ball to ball. Found
none, it refuses; so a wall at or past the body's inradius never reaches SMLib, whatever the
body's shape. `smlib_offset_faces` moves the picked faces along their outward normals onto
their surfaces' offsets (`SmOffsetGeometryCreation::FaceuseOffset`) and lets SMLib's local
operation (`SmBrep::LocalOperation`) extend or trim the faces around them; moved in past all of
the body behind them is `SMLIB_ERR_NO_ROOM` by the same test. `smlib_thicken` makes a solid of a
sheet with `SmApiShellBrepFull` (ruled walls along the open edges), a symmetric one after
moving the sheet back by half with `SmApiOffsetBrepFull`.

`tests/ops.c` checks them from C, `tests/ops_kernel.cpp` with SM_API-made cylinders and SMLib's
regression cases, both against analytic volumes and `smlib_brep_check` level 2.

exact mass properties (`smlib_brep_mass_properties`, `smlib_brep_area`), history, the flat
arrays and surface/surface intersection below. SolidCraft's kernel adapter adds the rest, one
shim at a time, each with a case in a fuzz target.

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
- Edge curves are exported as lines, circles and ellipses when SMLib holds them as `SmLine`,
  `SmCircle` and `SmEllipse`. B-splines stay NURBS, and any other curve kind is approximated
  as NURBS to the B-rep's tolerance. (This differs from SmuConvert, which exports every edge
  curve as NURBS; its import already read the analytic kinds.) Analytic edges are written in the
  arrays and units that import reads: the range is in radians for circles and ellipses. Before
  writing an analytic edge, export rebuilds the curve exactly as import will and checks it
  pointwise; if the check fails, the edge stays NURBS. An imported analytic edge has the
  parameterisation of its kind rather than the original's. On a NURBS face, the edge's UV trims
  are re-mapped to match, or dropped when the two parameterisations are not affinely related,
  in which case SMLib rebuilds them.
- SMLib's own primitives (`SmApiCreateCylinder`, `SmApiCreateCone`, `SmApiCreateTorus`) build
  their circles as plain B-splines, so those export as NURBS. So do booleans' intersection curves.
- A cone with a negative semi-angle is rebuilt in the equivalent positive-angle form on import.
  Its second export therefore differs from its first, and from then on the arrays are stable.
- Attributes, materials and colours are not carried.

Import validates everything first (counts, indices, lengths, kinds, finite values, NURBS
orders and knots, closed loops, radial cycles), so malformed arrays return
`SM_ERR_INVALID_INPUT` without building anything. A failing import's dump holds its arrays as
text, which `smlib_flat_read` reads back for replay. `tests/flat_roundtrip.cpp` round-trips a box,
cylinder, cone, sphere, torus, two filleted boxes, a boolean and an inclined cut, checking each
one's exported curve kinds, then feeds the import hostile
and randomly mutated arrays.

## Surface/surface intersection

`src/smlib_ssi.cpp` intersects two surfaces described in plain C, without a B-rep:

- `smlib_surface_plane`, `_cylinder`, `_cone`, `_sphere`, `_torus` and `_nurbs` take exactly the
  fields of an SSI corpus (`ssi-corpus`) case file (`cases/data/<id>.json`) and build the surface the way the corpus runner
  (`runner/ssi_runner.cpp`) does: the same SMLib constructors, parameterisations and domains.
  Angles are in degrees, as in the case files. Every argument is validated first (finite,
  bounded, non-parallel axes, positive radii, consistent counts, non-decreasing knot vectors with
  a non-empty domain, positive weights), and anything else returns `SMLIB_ERR_ARG`.
- `smlib_ssi(a, b, method, approx_tol, samples_per_curve, &result)` runs
  `SmSurface::GlobalSurfaceIntersect` (0), `SmAdvSurfaceIntersector` (1) or
  `SmSurfaceIntersector` (2), as the runner's `gsi`, `adv` and `basic` modes do, with both
  surfaces' edges seeding curves and `approx_tol <= 0` meaning SMLib's default.
- Each curve is sampled as the runner samples it: `samples_per_curve` parameters spaced evenly
  over the 3D curve's natural interval, each UV curve likewise over its own interval, and both
  surfaces evaluated at the UV samples. The result gives per curve its SMLib type (crossing,
  tangent, touching and so on), whether the method returned UV curves (GlobalSurfaceIntersect's
  analytic special cases return none), and the samples copied into caller buffers. TOUCHING curves
  are listed again as isolated points.
- A failing call's dump holds both surfaces as a corpus case file, `<n>-smlib_ssi-case.txt`, that
  the corpus runner replays.

Through this API, all 119 corpus cases give bit-for-bit the runner's output (curve counts,
types and every sample), in all three modes. One difference: where SMLib cannot evaluate a
sample, the API writes NaN, while the runner drops the sample. `tests/ssi.cpp` checks corpus-like
cases against their closed forms and rejects hostile inputs; `fuzz/fuzz_ssi.cpp`
(`smlib_fuzz_ssi`, run with `-timeout=10`) fuzzes the constructors and all three methods.

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

`fuzz/fuzz_ops.cpp` (`smlib_fuzz_ops`) chains the operations on boxes, half the values pulled
into ranges SMLib can build so the kernel runs too, the rest arbitrary.

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
- **ProtoEdge left in its ProtoVertex** (`SmProtoTopology.cpp`, `SmProtoEdge::~SmProtoEdge`): a
  ProtoVertex lists a ProtoEdge once per end, so one bounding both ends of an edge (closed, or
  merged by `CombineProtoVertex`) lists it twice, but the destructor removed it once. The
  ProtoVertex's destructor then wrote through the freed ProtoEdge (heap use-after-free in
  `~SmProtoTopologyManager`, ASan, on box − cylinder and on shells). Every entry is removed now.
- **Implanted curves deleted** (`SmProtoTopology.cpp`, `SmEdgeDefinition::Destruct` and
  `Imprint`): `MakeEdgeInFace` can fail after it has made a new Edge of the curves it was given;
  the EdgeDefinition then deleted them anyway, and the Edge used them after the free. Curves with
  an owner (which `MakeEdgeInFace` sets) are left to their Brep now.
- **Failed imprint ignored** (`SmTopologyIntersector.cpp`, `SmMerge.cpp`):
  `IntersectInsertRelate` overwrote the status of `ImprintIntersections` with the next step's,
  and `SmMerge::PiecewiseMerge` ignored its first `IntersectInsertRelate`, so a boolean went on
  with half-imprinted, inconsistent Breps (use-after-free reads in later steps). Both stop with
  the error now.
- **Vertex used after a squeeze** (`SmProtoTopology.cpp`,
  `SmProtoTopologyManager::GlueCoincidentBrepVertices`): for coincident vertices joined by
  several edges it went on through the edges after squeezing or gluing had already deleted the
  second vertex, gluing and squeezing with the freed vertex. It stops once that vertex is gone.
- **Matrix type read uninitialised** (`NL_Math.cpp`, `N_InitRealMatrix`): it left the matrix type
  unset, and `N_RealMatrixIsNULL` reads it (UBSan: invalid enum load, from fillets that join a
  tangent chain). It is set now.
- **Removed ProtoEdge left in ProtoFaces** (`SmProtoTopology.cpp`,
  `SmProtoTopologyManager::RemoveProtoEdge`, `SplitIntersectingProtoEdges`): a ProtoEdge
  combined into another and deleted stayed in the coincident-edge lists of ProtoFaces (and the
  overlap lists of other ProtoEdges), which `GlueCoincidentBrepVertices` then read. Removing a
  ProtoEdge drops those references; a combined one is replaced by the one it went into.
