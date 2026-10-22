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
exact mass properties (`smlib_brep_mass_properties`, `smlib_brep_area`), the exact box
(`smlib_brep_bounds`), history and the flat
arrays below. SolidCraft's kernel adapter adds the rest, one shim at a time, each with a case in
a fuzz target.
history, the flat arrays below and the operations. SolidCraft's kernel adapter adds the rest,
one shim at a time, each with a case in a fuzz target.

## Operations

`src/smlib_ops.cpp`: `smlib_boolean` (union, subtract, intersect), `smlib_fillet`,
`smlib_chamfer`, `smlib_blend`, `smlib_fillet_variable`, `smlib_shell`, `smlib_offset_faces` and `smlib_thicken`. Each copies its inputs and works on the copies, so tagged
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

`smlib_blend` gives each picked edge (and its tangent chain) its own size and one of four
blends: a round (circular section, constant radius, as `smlib_fillet`), a curvature-continuous
round (SMLib's G2 blend-curve section between the contact lines of the round of that radius), a
chord round (circular section at a constant width from contact to contact, SMLib's
constant-distance solver) and the chamfer above. `smlib_fillet_variable` rounds one edge with a
radius that follows a clamped B-spline law (degree 1 to 3) over the edge's parameter, 0 at the
end nearest a given point; it becomes an `SmBSplineFilletLaw` on SMLib's variable-radius solver,
whose law curve carries the radius as x and, as SMLib's own laws do, the parameter itself as y
(control points at the Greville abscissae), so a constant law can be extended past the edge's ends.
Where the law has a corner the round turns and is split into faces there; a law whose radius is
zero inside the edge (the round would pinch to a point) is `SMLIB_ERR_UNSUPPORTED`.

`smlib_merge_faces` removes a body's topological edges and vertices
(`SmBrep::RemoveTopologicalEdgesAndVertices`): faces on one surface that meet along an edge
become one, and two edges on one curve meeting at a vertex nothing else uses become one; the
edges under the given points stay, so faces split on purpose stay apart. SMLib's booleans do the
same at their end (`SmMerge::BooleanPostProcess`); this is it for any other body. History tags
follow the merge (a merged face carries both faces' tags).

`tests/ops.c` checks them from C, `tests/ops_kernel.cpp` with SM_API-made cylinders and SMLib's
regression cases, both against analytic volumes and `smlib_brep_check` level 2.

exact mass properties (`smlib_brep_mass_properties`, `smlib_brep_area`), the exact box
(`smlib_brep_bounds`), history, the flat
arrays and surface/surface intersection below. SolidCraft's kernel adapter adds the rest, one
shim at a time, each with a case in a fuzz target.

## Sweeps

`src/smlib_sweeps.cpp`: `smlib_extrude` (straight through `SmApiCreateLinearSweep`, drafted
through `SmApiCreateTaperExtrude`), `smlib_revolve` (`SmApiCreateRotationalSweep`), `smlib_sweep`
(`SmApiCreateCurveSweep`, the path's curves joined into one B-spline, the profile carried as it
lies relative to the path; `smlib_sweep_placed` also says where it placed the profile at its
first, middle and last sections; round a closed path the walls are sewn to each other where the ends
meet, with no caps; along a straight path `SmApiCreateLinearSweep` and, turning with arcs about
one axis, `SmApiCreateRotationalSweep`, whose walls are exact
surfaces where the curve sweep fits B-splines), `smlib_loft`
(`SmApiCreateSkinPrimitive`: a section's loops, outline then holes, are one profile, each lofted to
the loop in the same place in the next section), `smlib_loft_to_point` (a pyramid: the planar
base and a triangle from each of its lines to the apex, stitched and oriented as
`SmApiCreatePyramid` builds its own), `smlib_screw_sweep` and `smlib_sweep_stations`.
`smlib_project_curves` projects a curve list onto one picked face along a direction
(`SmBrep::CreateParallelProjectionCurves` restricted to that face, as `SmApiProjectCurve` runs it
on all of them) and hands back the projected curves as the B-splines SMLib made
(`smlib_curves`), with their intervals. `smlib_face_section` gives one face's section by a plane
(`SmBrep::CreatePlanarSectionCurves`, the curves on that face), also as `smlib_curves`.
Profiles, paths and sections cross as curve lists (lines, circular arcs, NURBS; `include/smlib.h`
gives the layout), checked before SMLib sees them and built into SMLib curves inside the call.
With `cap` the result must be a closed manifold solid (`SMLIB_ERR_NOT_SOLID` otherwise).

Two of SM_API's behaviours shape the shims. Its draft extrusion takes the draft's sign for the
side (positive in, negative out), not the 90 to 180 degrees `SmApiBrep.h` documents for an
outward draft, and asserts a single chain of curves, so a drafted profile must be one loop. Its
curve sweep carries the profile on its own frame along the path, which along a helix turns the
profile away from a screw motion (3 turns of a thread's tooth came out 1 to 2% off in volume);
`smlib_screw_sweep` therefore places the profile by the screw motion itself at stations along
the turns and skins the walls through them, cubic, the stations close enough for the caller's
tolerance. SM_API's curve sweep has no twist and its scale law is uniform in the path's parameter,
and nothing in SM_API follows guide rails round a closed profile; `smlib_sweep_stations` takes the
caller's stations instead (affine maps: the profile twisted, scaled or stretched and carried
along the path), places the profile exactly by each (its control points mapped), skins the walls
through the sections (`SmApiCreateSkinPrimitive`, cubic) and caps them, or, round a closed path
where the last section lies on the first, glues the walls' ends to their starts
(`SmApiStitchIntoSolid`).

`tests/sweeps.cpp` checks each against closed-form volumes and SMLib's check, and feeds hostile
curve lists; `fuzz/fuzz_sweeps.cpp` (`smlib_fuzz_sweeps`) fuzzes all seven.

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
- Tolerances (an addition to SmuConvert): besides the B-rep's one `TOLERANCE`, export writes
  each edge's and vertex's own (`EDGE_TOLERANCE`, `VERTEX_TOLERANCE`, as twice its zone
  tolerance), and import gives each edge and vertex its own when they are present. So one loose
  edge (a fillet's approximated rail) no longer makes the whole body loose on the way back.
  Building an edge again may raise its tolerance (SMLib measures what it needs), never lower it.
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
- Shell order (an addition to SmuConvert): SMLib takes a region's first shell as the one that
  encloses it, and a body whose first shell is an inner one fails its check (test 12, "inner
  shell is not geometrically inside its containing outer shell"; tangency-labs/solidcore#388).
  Export writes the enclosing shell first. Import takes a region's shells in any order and puts
  first the one that SMLib's containment test (`SmShell::IsInnerShellOf`, the ray casts of
  test 12) finds enclosing all the others, trying the largest bounding box first; a region
  where none can be shown to enclose the others keeps its order.

Import validates everything first (counts, indices, lengths, kinds, finite values, NURBS
orders and knots, closed loops, radial cycles), so malformed arrays return
`SM_ERR_INVALID_INPUT` without building anything. A failing import's dump holds its arrays as
text, which `smlib_flat_read` reads back for replay. `tests/flat_roundtrip.cpp` round-trips a box,
cylinder, cone, sphere, torus, two filleted boxes, a boolean and an inclined cut, checking each
one's exported curve kinds, imports an island inside a void with its regions and shells in
every order, then feeds the import hostile
and randomly mutated arrays.

## Surface/surface intersection

`src/smlib_ssi.cpp` intersects two surfaces described in plain C, without a B-rep:

- `smlib_surface_plane`, `_cylinder`, `_cone`, `_sphere`, `_torus` and `_nurbs` take exactly the
  fields of an SSI corpus (`ssi-corpus`) case file (`cases/data/<id>.json`) and build the surface the way the corpus runner
  (`runner/ssi_runner.cpp`) does: the same SMLib constructors, parameterisations and domains.
  Angles are in degrees, as in the case files. Every argument is validated first (finite,
  bounded, non-parallel axes, positive radii, consistent counts, non-decreasing knot vectors with
  a non-empty domain, positive weights), and anything else returns `SMLIB_ERR_ARG`.
- `smlib_ssi(a, b, method, approx_tol, &result)` runs `SmSurface::GlobalSurfaceIntersect` (0),
  `SmAdvSurfaceIntersector` (1) or `SmSurfaceIntersector` (2), as the runner's `gsi`, `adv` and
  `basic` modes do, with both surfaces' edges seeding curves and `approx_tol <= 0` meaning
  SMLib's default.
- Each curve is kept exactly as SMLib holds it: `smlib_ssi_curve_nurbs` copies its 3D curve and
  its UV curves on `a` and `b` as (rational) B-splines (degree, every knot, cartesian control
  points, weights). SMLib's lines, circles and ellipses are B-spline curves, and its intersection
  curves B-splines it fitted within `approx_tol`, so nothing is sampled or approximated here; a
  curve SMLib does not hold as a B-spline is `SMLIB_ERR_UNSUPPORTED`, a UV curve the method did
  not return `SMLIB_ERR_EMPTY`. The result also gives per curve its SMLib type (crossing,
  tangent, touching and so on), SMLib's deviation estimate, the 3D curve's natural interval and
  whether the method returned UV curves (GlobalSurfaceIntersect's analytic special cases return
  none). TOUCHING curves are listed again as isolated points, at the start of their curves.
- A failing call's dump holds both surfaces as a corpus case file, `<n>-smlib_ssi-case.txt`, that
  the corpus runner replays.

A caller that compares with the corpus runner's output evaluates the curves itself where the
runner samples them (evenly over each curve's natural interval). `tests/ssi.cpp` checks
corpus-like cases against their closed forms, evaluating the exact curves, and rejects hostile
inputs; `fuzz/fuzz_ssi.cpp` (`smlib_fuzz_ssi`, run with `-timeout=10`) fuzzes the constructors and
all three methods.

## Planar regions

`smlib_planar_faces` (`src/smlib_sweeps.cpp`) makes a sheet of planar faces from a curve list of
closed loops in one plane, in any order (`SmApiCreatePlanarFaces`, which nests them: outer loops,
their holes, islands in the holes). `smlib_boolean2d` (`src/smlib_ops.cpp`) combines two such
sheets on one plane as regions of it (`SmApiBoolean2d`): SMLib's topology intersector splits both
inputs' edges where they cross, and the result's edges are copies of the inputs' own curves
(lines, circles and NURBS, none refitted), so the regions read back exactly.

## Imprint

`smlib_imprint(body, tool, &out)` cuts each face of `body` that `tool` meets along the
intersection (SMLib's boolean with `SM_BO_IMPRINT`, the tool kept), leaving one body whose faces
are the pieces. Trim and Split Face build on it and choose which pieces to keep; a tool that
meets nothing is `SMLIB_ERR_NO_CUT`.

## Split

`smlib_split_faces(body, faces, n, tools, n_tools, &out, made)` cuts the chosen faces where
each tool surface (an `smlib_surface`, over its whole domain) meets them, as SMLib's plane cutter
cuts a subset of faces (`SmPlaneCutter::MergeIntersection`): each face's surface is intersected
with the tool over the face's UV domain (`CreatePlanarSectionCurves` for a plane,
`GlobalSurfaceIntersect` otherwise) and the curves go to `SmBrep::MergeCurvesOnSurface` with the
tool as the intersecting surface, which trims them to the face and splits it. Unlike the plane
cutter, a failed merge fails the call. `made[i]` counts the faces face `faces[i]` became besides
itself; a call that splits nothing is `SMLIB_ERR_NO_CUT`.

`smlib_split_body(body, tool, &back, &front)` slices a manifold solid with a one-face sheet on
the tool surface (`SM_BO_SLICE`), once as it is and once turned over (`SwapUV`; the slice reads
the side from the surface's normal, not from the face's use), giving the part behind the surface
and the part in front of it. A tool that misses leaves both parts the whole body.

`smlib_face_points(body, out, n_faces)` gives a point inside each face
(`SmFace::CalculateAnInternalPoint`), so a trim can tell which side of its tool each piece lies on.

## Extend

`smlib_extend_sheet(sheet, point, tol, distance, mode, &out)` continues a one-face sheet past the
side its edge nearest `point` runs along. The face must be its surface's whole parameter box;
anything else is `SMLIB_ERR_UNSUPPORTED`. The surface is SMLib's extended surface
(`SmBSplineSurface::CreateExtendedSurface`: `SMLIB_EXTEND_NATURAL` continues its end spans,
`_TANGENT` and `_CURVATURE` add a G1 or G2 run-off). SMLib's analytic surfaces are B-splines too
(subclasses of `SmBSplineSurface`), so they take the same path; `CopyAndAddAnalytics` gives the
result back its analytic kind where it is one.

## Untrim

`smlib_untrim_face(body, point, tol, holes, n_holes, outer, in_place, &out)` restores the face
nearest `point` towards its surface. Holes picked by points (or all, `n_holes` < 0) are deleted
with `SmFace::DeleteInnerLoop` on a copy where their edges are lamina. With `outer` the face is
rebuilt over its surface's parameter box (`SmBrep::CreateFaceFromSurface`: a free-form surface's
own box, an analytic one's the tight box of its outer loop's UV trims; one all the way round
a periodic direction closes on the surface's seam) and each hole kept is
rebuilt in it from its own edges (`MakeVertexLoop`, `MakeEdgeInFace`, the disc it splits off
deleted, as `CreateFacesWithCurves` does). The face alone comes back, or with `in_place` the body
with the face replaced; the new face and edges carry the attributes (history tags) of the ones
they stand for.

## Sheets

`src/smlib_sheets.cpp`: `smlib_stitch` sews a body's faces along the vertices and edges that
coincide within the caller's tolerance (`SmBrep::StitchFaces`, once, at that tolerance; not
`StitchAndOrient`, which caps the tolerance and stitches again at a larger one), then orients them
with `SmBrep::OrientTrimmedSurfaces`: outward where every edge is glued, else only so neighbours
agree. It reports whether the result is a closed manifold solid. `smlib_offset_sheet` moves every
face of a sheet onto its surface's offset along its normal (`SmApiOffsetBrepFull`, as
`smlib_thicken` moves a sheet back by half). Tags travel through both.

## Diagnostics

- **Log:** `smlib_set_log` receives one line per call, written after it returns. Example:
  `smlib_make_box(origin=[0,0,0], size=[10,20,30]) -> 0 (0.21 ms)`.
  It also receives each warning, message and error SMLib reports itself through
  `smos_ErrorMessage` (`WARN`, `SE`, explicit reports), tagged with the call in flight: the
  shims install SMLib's error callback (`smos_SetErrorCallback`, kept per thread by
  `SmThreadLocalStorage`) on the calling thread at the start of every call.
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
- **A curve sweep that keeps the profile where it lies** (`SmPrimitiveCreation.h`,
  `SmSweepOptions::bProfileAtPathStart`; `SmPrimitiveCreation::CreateCurveSweep`): upstream
  sweeps a profile with its first point (or its bounding-box centre) moved onto the path, so a
  profile beside the path loses its offset round a bend. With the new option, off by default,
  the path's start is the profile's reference point: the profile is carried as it lies relative
  to the path. `smlib_sweep` uses it.
- **Knot checks before clamping** (`NL_Iges.cpp`, `N_IgesClampSrf`): it assumed finite,
  nondecreasing knots and bounds inside them. A NaN knot left a span index unset, and the knot
  copy then wrote past the new knot vector, a heap overflow. A surface of revolution given a
  sweep domain from 360 degrees made such a knot. Such knots and bounds are now refused
  (`NL_YES`), so the caller's operation fails instead of corrupting memory.
- **An arc moved past its turn** (`SmEllipse.cpp`, `SmEllipse::AdjustSTEPInterval`): it found the
  new interval's NURBS parameters through the polar converter of the current arc, which clamps an
  angle outside that arc to its ends. A full circle given [360, 540] degrees (a torus face whose
  minor range starts a turn on) had both ends clamped to one parameter, so the circle had no
  parametric length, and the torus rebuilt on it had NaN knots (the overflow above, then a failed
  import). An angle outside a full-turn arc is now taken whole turns back into it, the NURBS
  parameter moving by the circle's NURBS length per turn; outside a partial arc, the arc's NURBS
  interval continues in proportion to the angle.
- **Caps of a curve sweep whose path bends** (`SmPrimitiveCreation.cpp`,
  `SmPrimitiveCreation::CreateCurveSweep`): it took the start cap's curves from the first side
  faces and the end cap's from the last, one per profile curve. A swept surface is split into a
  face per G1 piece of the path (a line, then an arc), so with a profile of two loops (an
  annulus) the start cap got a face from the junction and the end cap missed one, and the sweep
  failed (1001). Each side face is now marked by whether its own surface starts or ends where
  the swept surface does along the path, and the caps take those.
- **Intersection edges that run past their vertex** (`SmProtoTopology.cpp`,
  `SmEdgeDefinition::Imprint`): an intersection curve's ProtoVertex can come from another
  intersection, and where the curve crosses a face's edge at a grazing angle (a thread groove's
  flank meeting a shaft's rim at the end face) that vertex lay on the curve a few tolerances inside
  its end. The edge ran past its vertex and lay along the rim within tolerance up to its end
  (SmFace test #34, "partially coincident edges"), with a vertex gap that only a larger tolerance
  covered. Each end of the edge's interval now moves to the curve's point nearest its ProtoVertex
  when that point is inside the interval, within tolerance and nearer than the end, and the UV
  curves are trimmed with it. Where the vertex instead lies past the curve's end (the curve stops
  short of it, as a thread flank's curve on a nut's bore does 3.1e-5 before the bore's rim), the
  curve and its UV curves are continued to the curve's point nearest the vertex, each by the exact
  polynomial of its end span (`sm_CreateExactlyExtendedCurve`), when that point is within tolerance
  of the vertex and each UV curve's new end lies on its face there. The point each end goes to is
  the ProtoVertex's position only when the ProtoVertex is one Vertex (every vertex definition in
  it imprinted onto the same Vertex in each Brep); otherwise it is the Vertex the end is
  imprinted onto, where that is one point in both Breps. A ProtoVertex collects every vertex
  definition within the capture distance (0.01 mm), which may resolve into several Vertices, and
  its position is then their mean: a pleated loft joined to a slab, whose ruled edge crosses the
  slab's window top 0.0042 mm from the window's corner, had its curve along the window top,
  already ending at its own Vertex, continued 0.0021 mm past it to that mean, which left a gap
  of that size at the next edge (SmLoop, "Edge/Edge Gap").
- **A moved face that leaves its tangent neighbour** (`SmOffsetGeometryCreation.cpp`,
  `ExtendTrimEdges`, `EdgeuseOffset`): a local operation (`SmBrep::LocalOperation`) rebuilds a
  tangent edge by lifting it onto one face's new surface and dropping it onto the other's. When
  only one of the two faces moved (a boss wall beside its base round, one piece of a split
  face), the drop crossed the gap, the failure was ignored, and the body came back with an edge
  as loose as the gap. A drop farther than ten times the edge's tolerance now fails the local
  operation (offset and shell, which close such gaps later, are unchanged).
- **An inward shell with no cavity** (`SmOffsetExecutive.cpp`, `ShellBrep`): an inward wall at
  or past the body's inradius leaves an empty offset (or one that is no closed solid), and the
  difference then handed the body back unhollowed as a success. Such a shell now fails
  (`SM_ERR_INVALID_INPUT` for an empty offset, `SM_ERR` for an open one).
- **Offset edges on a closed intersection curve** (`SmOffsetGeometryCreation.cpp`,
  `TrimExtendedEdge`): two trim points on a closed curve bound two arcs, and the edge always took
  the interval between them. The rim of a hemisphere, two arcs of one circle, became the same
  long arc twice and the cavity no solid. The edge now takes the arc that runs where its base
  edge does when that arc is an interval of the curve (a trim point on the seam), and reports no
  curve when it would cross the seam.
- **Pointer maps in insertion order** (`SmMapPtrToPtr.h`, `SmMapPtrToPtr`, `SmMapPtrToPtrs`; through
  them `SmRelation`): they were `std::unordered_map`s and `std::unordered_set`s keyed by
  pointers, and handed out their keys and values (`GetAllKeys`, `GetSeconds`, `GetSecond`) in
  hash order, which is address order. The offset and shell walk them to process edges and
  faces, so the same call could give a different answer wherever the heap put its objects (a
  box with a cylinder on its side shelled to 1001 or 1006 by turns in one process, and a host
  whose threads allocate meanwhile moves the heap with its load). They now hand everything out
  in the order it was inserted, so every run does the same work.
- **A validity check linear in the body** (`SmLoopuse.cpp`, `SmTopology.cpp`/`.h`, `SmEdgeuse.cpp`,
  `SmFace.cpp`, `SmBrep.cpp`, new `SmBoxSearch.h`): `SmBrep::ValidatePointers`, `AssertValid` and
  `ContainsCoincidentTopology` took time quadratic in a body's topology (over ten minutes on a
  50,000-face part). A loop's vertices and edges are de-duplicated with a hash set past 32 of
  them; `SmFace::PointClassify` no longer lists the face's vertices to test the first one at the
  point; within `AssertValid` and `ContainsCoincidentTopology` (`SmTopologyListWalkScope`) each
  topology list, each loop's edgeuses and each face's vertices are gathered once instead of once
  per member; and the pairwise tests (coincident vertices, edges, faces, vertex/edge, vertex/face
  and edge/face; a face's coincident vertices and edges) visit only the pairs that a box search,
  or a shared vertex or edge, leaves possible, in the order the all-pairs loops did. Every test
  still runs on every pair that could fail it, so the reports are the same.
- **Trimming curves of a degenerate face** (`SmBrep.cpp`, `CreateFaceInRegionFromSurface`): it
  made a face's boundary curves and then let `MakeFaceWithCurves` refuse a surface degenerate to
  a curve or a point, which returns before it takes the curves, so they leaked (a box at
  `z = 1e200`). The surface is refused before the curves are made.
- **Linear sections on arcs revolved** (`SmFilletSolver.cpp`, `CalcTorusFilletGeom`): a round
  on a circular arc between surfaces of revolution about the arc's axis was made exactly, as the
  revolved section (a torus), but a linear section (a chamfer) was always traced as the
  intersection of the offset surfaces. Where the offset of one face collapses (a chamfer as wide
  as the corner round it runs round, whose offset is the round's axis), that intersection is a
  single point and the trace failed (SM_ERR). A linear section is now revolved too, a cone.
  The same trace failed on a hole's rims under `SmApiChamferFillet` (constant distance), which
  now gives a closed solid there.
- **Rail points of a section centered on the axis** (`SmFilletSolver.cpp`,
  `CalcTorusFilletGeom`): the revolved section was placed by the rail points of the
  intersection points at its ends, and oriented by the edge point nearest the section's
  center. A round as wide as the corner round it runs round (R5 rim beside R5 vertical
  corners) is centered on that round's axis, where the offset faces collapse: the rail points
  then lie anywhere on the circle the ball touches, and every edge point is as near, so the
  torus could turn the long way round, and the rollover test failed (SM_ERR). With the center on
  the axis, the rail points are now turned into the half-planes through the edge's ends, and
  the edge's start orients the sweep.
- **Rounds whose center curve is a point** (`SmFilletSolver.cpp`, new `CalcSphereFilletGeom`):
  the round of an edge between two cylinders as wide as the round, whose axes meet (the seam
  where two R2 bottom rounds of a box meet, in the tangent chain of an R2 vertical edge), has
  every ball centered where the axes meet: the offset cylinders are their axes, and their
  intersection, which SMLib traced, is that one point (SM_ERR). Such a round is now made exactly,
  as the patch of that ball between the circles it touches the cylinders along.
- **Radial order at a tangent edge** (`SmEdge.cpp`, `SmEdge::AssertValid`, test 9): it ordered
  the faces round an edge by their binormals alone. Faces tangent along the edge (a cylinder or a
  torus a boolean leaves touching a plane, four faces round the tangency) have the same binormal,
  so the check read their right order as out of sequence and refused a valid body. Where two
  faces round an edge leave it within 0.5° of each other, it now orders every face by the
  direction to its point 20 edge tolerances into the face, square to the edge.
- **Analytic copies keep the parameterization** (`SmBSplineSurface.cpp`,
  `SmBSplineSurface::CopyAndAddAnalytics`): splitting a face copies its surface for the new face,
  as an analytic one where the shape is one, and the new face keeps the old UV domain and trim
  curves. The analytic types build their NURBS form in their own parameterization (a circle as
  SMLib's two-span arcs) and match only the domain, so a NURBS face parameterized otherwise (a
  rational half cylinder in four spans) came out of a boolean with trims millimetres off its edges
  (edge tolerances of 0.78 mm, ValidatePointers failing). An analytic copy that does not give the
  same points at the same parameters is replaced by an exact copy.
- **A filleted vertex with no corner** (`SmFilletExecutive.cpp`, `CreateFilletCorners`,
  `CollectTangentEdgeSequence`; `SmFilletCorner.cpp`, `SmFilletCorner::Create`): the vertex/corner
  map holds each filleted vertex's array of edgeuses until its corner replaces it. Where no corner
  type handles a vertex (four filleted edges among five, at the line where a hole's wall touches
  a face, its rim's tangent chain running on along that face), the array stayed in the map, and
  the smooth-sequence walk called it as a corner through its virtual table (a crash). Such a
  vertex now leaves the map, its array and the corner's edgeuse list are freed, and the fillet
  fails (SM_ERR); a walk that finds no corner at a vertex fails too.
- **A round ending past a tangent side edge** (`SmFilletCorner.cpp`,
  `SmFilletNx1Corner::MakeCornerTopology`; `SmFilletGeom.cpp`, `SmFilletEdge::CalcEdgeOnExtendedEdge`):
  at a vertex where the filleted edge's side face runs on tangentially across a side edge (a flat
  top onto a rounded end) and the other side edge has the opposite convexity (a slot's floor
  meeting its wall), the corner took the tangent edge for the side edge to extend and found no end
  face between it and the other side edge, so it used the round, which the fillet never meets
  (1001). The corner now skips tangent side edges round the vertex: when the edge reached past
  them has the filleted edge's convexity, it is the side edge to extend and the end face is the
  face between it and the other side edge (the slot wall). That extended edge runs in the side
  face, so its curve is where the side face meets the end face's extension rather than the
  extension of the edge's own two faces.
- **A radius law that cannot be evaluated** (`SmFilletStandardSolver.cpp`,
  `SmVariableRadiusFS::SetupOffsetValues`; `SmFilletExecutive.cpp`, `SmTemporaryRadiusChange::SetUp`;
  `SmFilletSolver.cpp`, `StepSolve`): both ignored the status of the law's evaluation, which fails
  where the law must be extended past the edge's ends and its curve can't be (a constant law
  curve of zero speed), and set the offset surfaces to an uninitialised radius. They now return
  the failure (`SetUp` changing nothing), and the stepping checks `SetupOffsetValues`.
- **Where a curve sweep places the profile** (`NL_SrfAdv.cpp`, `N_CreateSweepScale`;
  `SmBSplineSurfaceNLib.cpp`, `SmBSplineSurface::CreateSweepSurface`; `SmPrimitiveCreation.cpp`,
  `SmPrimitiveCreation::CreateCurveSweep`; `SmApiBrep.cpp`, `SmApiCreateCurveSweep`): the sweep
  carries the section on local systems it computes from point to point along the path (each
  z the last one taken square to the next tangent), so along a path in space nobody else can say
  where the profile ends up, and a caller sweeping a path with corners piece by piece could not
  start the next piece where the last one left the profile. Each now takes an optional output
  (default none): the systems of the first, middle and last sections it places, and, from
  `CreateCurveSweep` up, the rigid motions carrying the profile as given onto them.
  `smlib_sweep_placed` hands them on.
- **A planar path's sweep frame** (`SmPrimitiveCreation.cpp`, `SmPrimitiveCreation::CreateCurveSweep`):
  the sweep carries the profile on a frame that starts as the one `MakeUnitOrthoVectors` builds on
  the path's start tangent. For a start tangent along +y its z lies in a planar path's plane, and a
  wall swept along a half circle came out with its flat top spread 2.9e-4 mm about its plane (the
  tolerance is 1e-5), so a boolean with another wall on that plane failed (1006). A planar path now
  starts the frame with its z on the plane's normal, and the top is flat to 1e-13.
- **Coincident faces merged across a fold** (`SmSurface.cpp`, `SmSurface::CoverCoincidentSurface`):
  after a boolean, the faces either side of an edge between coincident surfaces are merged onto one
  surface, a B-spline continued past its natural domain to cover the other (a planar one up to a
  hundred times its size). That continuation need not be one sheet: a swept planar sector's
  straight lines run on through its centre and back over it. Points of the other face dropped onto
  the wrong sheet, and the merged face's UV trims lay far from its edges (two swept arc walls
  sharing their top and bottom planes, joined: UV gaps and the wrong volume). A B-spline now covers
  only where its continuation does not fold over the domain the merge uses: the normal may not
  vanish or turn back between neighbouring points of a 21 x 21 grid over that domain; otherwise
  the faces stay as they are, with the edge between them.
- **Inward offsets across rounds thinner than the wall** (`SmOffsetExecutive.cpp`/`.h`,
  `DoExtendedOffset`, new `GetEdgeuseAcrossCollapsedFace`; `SmOffsetGeometryCreation.cpp`,
  `TrimExtendedEdge`): a convex round whose radius is at most the inset distance has no offset
  (its offset surface is NULL), and the offset, so the shell, failed (1001): a box with r1 or
  r0.5 edges shelled 2 inward (tangency-labs/solidcore#397). Such a face now gets no offset face;
  the faces on either side of it are mated across it (their offsets intersected, as at a sharp
  convex edge, the trim points found on the offsets of the faces beyond collapsed neighbours), and
  the offset edge along its end, which the mating shrinks to a corner, is squeezed out. A collapsed
  face that cannot be bridged so (not four-sided, or a collapsed face across it) is refused with
  `SM_ERR_DEGENERATE_SURFACE`.
- **Triangles held to the chord tolerance** (`SmTess.cpp`/`.h`, new
  `SmTess::RefineTrianglesToChordHeight`, called from `AddQuadBoundaries`): the subdivision tree
  bounds each cell's chord height, but the triangles made in it were never measured. The cell strip
  beside a sphere's seam, filled by ear clipping without interior vertices from a reflex corner,
  became a fan of slivers 13 times the tolerance off the surface (r 12 at 0.05 and 0.01 mm), and
  the two triangles of a cell on a sphere or torus lay up to 1.26 times it off
  (tangency-labs/solidcore#390). Each triangle is now measured at its UV centroid and inner edge
  midpoints, along the surface normal, and one beyond the tolerance has its longest inner edge
  bisected (in both triangles on it, so the mesh stays closed) until none is; boundary edges and
  edges shorter than twice the tolerance are left alone. `tests/tess_kernel.cpp` checks spheres,
  tori and cylinders at 0.1 to 0.005 mm.
- **Variable rounds whose radius has a corner** (`SmFilletExecutive.h`/`.cpp`,
  `SmBSplineFilletLaw::GetCornerParameters`, `SmFilletExecutive::DoFillet`, `DoFilleting`,
  `NoteRadiusCornerSurfaces`, `SplitFilletFacesAtRadiusCorners`; `SmFilletSolver.h`,
  `SmFilletStandardSolver.h`/`.cpp`, `SmVariableRadiusFS::GetRadiusCornerParameters`,
  `IsAtRadiusCorner`, `GetRadiusParameter`, `SmCircularCrossSectionFSG::CreateSurface`;
  `SmFilletSolver.cpp`, `CalcFilletGeom`; `SmFilletIntersector.cpp`, `ComputeNextPoint`,
  `FlushCurve`): a radius law whose slope jumps inside the edge (straight 4 to 2 and back to 4)
  turns the round's center curve and rails there. The trace's tangent there was a central
  difference across the corner, the curves were fitted through it with that one averaged
  tangent, and one smooth surface was skinned across the turn: the round missed the law past
  the corner (35.1 mm³ removed for 40.06) and failed the check (control net doubling back).
  Now the corner is a through point of the trace, a step that would cross it ends on it, the
  spans on either side take the tangents of their own side, the surface is skinned piece by
  piece between corners and the pieces joined there, and once the round is in the body its
  face is split along the turn into smooth faces (each surface cut back to its face, the new
  edge the exact isoparameter curve).
- **Face subdivision counted its V knots as U** (`SmBrep.cpp`,
  `SmBrep::SubdivideFaceAtDiscontinuities`): the V pass took the number of U knots to split at
  and read past the V knots (a null read for a surface with a U discontinuity only).
- **The clean-up after filleting kept to the fillet** (`SmMerge.h`/`.cpp`, `SmMerge::BooleanPostProcess`;
  `SmFilletExecutive.cpp`, `DoFilleting`): it removed every topological edge in the body,
  between the target's own faces too. On a sphere in six patches (each on its own sphere
  surface, as an importer hands over a cube projected onto it) cut flat, an R1 round on the rim
  merged the five patches into one face whose surface no patch's parameterization could carry:
  loops of the wrong orientation, rails up to half a millimetre off their surfaces. It now
  removes no edge the target had and no edge between two of its faces, as booleans clean up
  only the edges of their intersection.
- **A round that runs over a hole** (`SmFilletExecutive.h`/`.cpp`, `DoFilleting`, new
  `RailsLeaveTheirFaces`): topology insertion, the executive's default, puts each round's rails
  into the faces they lie on and trims those faces back to them, and intersects the round with
  nothing else. An R4 round on a plate's edge whose contact line runs across a through hole
  must be trimmed by the hole's wall; it wasn't, and the body came back open. When a rail runs
  off the face it lies on (a point of it outside that face by both the 3D and the UV
  classification), the fillet faces are now merged into the body one at a time instead
  (`m_bDoPiecewiseMerge`), which intersects them with the rest of it.
- **Offset surfaces extended the same way whatever was evaluated first** (`SmOffsetSurface.cpp`,
  `sm_ExtendSurface`, new `sm_ExtendToward`, `SmOffsetSurface::EvaluatePoint`,
  `EvaluatePlainOrSimple`): evaluated beyond its base surface's domain, an offset surface
  extends the base (G1R), and the extension's shape depends on its length and on how its other
  sides were extended. It was made to reach one domain size past the first point evaluated out
  there, and extended again from itself for the next, so two copies of one offset surface
  evaluated in another order disagreed beyond the base: an R4 round on the rim where a half
  cylinder meets an extruded B-spline ramp found its end points on the fillet solver's copy and
  traced from them on another, 0.03 to 0.07 mm off, and failed (1006, then 1001). Now every
  extension is made from the base, all of its sides that can grow at once, each by the base
  domain's size times the least power of two that reaches the point, and never less than the
  extension it replaces: the same points give the same extension.
- **Precise face properties integrate each smooth piece on its own** (`SmFace.cpp`,
  `SmFace::ComputePreciseProperties`, `SmODEAreaIFEO::Evaluate`): area, volume and moments are an
  integral along each UV trim curve of an integral along the iso-u line from the face's lowest v,
  each was one adaptive Runge-Kutta run over its whole interval. The integrand's derivatives
  jump at the trim curve's interior knots, where the trim curve crosses a knot line u = k of the
  surface, and, along the iso-u line, at the surface's v knots (a rational circle's double knots, in
  either). A step across one is outside the method's error model, so its error estimate does not
  bound its error, and the requested accuracy did not hold (tangency-labs/solidcore#384, formerly
  KNOWN_ISSUES.md "faces bounded by circles"): at 1e-4 a sphere's area came out 1.6e-4 off and a
  torus's volume 2.5e-4, at 1e-3 a cylinder's volume 4.4e-4, at 1e-1 a cylinder's centroid 0.09 mm.
  Each piece between those breaks is now integrated on its own from a half-length first step, and
  the pieces summed (the crossings of u = k found per knot span by a sign change of u(t) - k between
  4 (degree + 1) points, then bisection). Against closed forms a cylinder, a sphere, a torus and a
  block with a round hole have area and volume within 2e-5 at 1e-4 and 5e-9 at 1e-8. Fixing only one
  of the two directions made tight accuracies worse, the remaining one's error no longer cancelled
  by the other's. 1.0.1 already split the trim curve at its own knots, with its retry from ever finer
  first steps; this splits at the surface's knot lines too, in both directions, and keeps the retry
  on each piece.
- **Faces an extended offset already mated are not intersected again** (`SmOffsetExecutive.cpp`,
  `SmOffsetExecutive::PiecewiseMerge`): the offset's self-intersection merge skips the surface
  pairs known to share an edge already (an offset face and its edge's cap face, faces across a
  tangent edge). With extended convex edges an edge gets no cap face: the offsets of its two
  faces are mated along it instead, and the merge intersected them again. On a cone frustum
  whose side is two half-cone faces, shelled open at the top, that plane/offset-cone curve ran
  along the inner rim the mating had made and ended 0.0036 mm short of the inner seam, so the
  rim was split there into a tiny edge and three spurious vertices, and the shell's difference
  then left partially coincident edges (SmFace test #34). Those mated pairs now count as
  sharing their edge too.
- **Which side a piecewise-merged fillet removes** (`SmFilletExecutive.cpp`, `DoPiecewiseMerge`;
  `SmAttribute.h`, `SM_AI_FILLET_EDGECURVE`): after merging the fillet faces into the body, the
  merge decides which region each fillet face cuts off by the side away from the cross
  section's centre of curvature, and skipped faces without curvature across them. A chamfer's
  faces are flat, so its wedge kept its material, the chamfer faces were removed again and the
  run returned success with the body unchanged (a 9 mm chamfer whose rails cross openings).
  The side is now the one the filleted edge lies on (a copy of the edge's curve, named from the
  fillet surface by an attribute the merge's copies carry), the same side on a round; and a
  merge in which no fillet face separates two regions fails (`SM_FILERR_BAD_MERGE`).
- **A split ProtoEdge's upper EdgeDefinition** (`SmProtoTopology.cpp`, `SmEdgeDefinition::Split`):
  the upper child is a copy of the EdgeDefinition, and the copy constructor gives it the
  original's owner, the lower ProtoEdge. `SmProtoEdge::AddEdgeDefinition` then removes it from
  that owner, which never listed the copy; once the removal's status was passed on, every split of
  a ProtoEdge that carries an EdgeDefinition failed (an inward shell of the ladder part failed in
  its offset's self-intersection merge). The copy now starts without an owner.
- **A rail that runs along its face's rim** (`SmFilletSolver.cpp`, `TestRollOver`): the rollover
  test classifies each rail against its face and treats every interval boundary as the rail
  crossing a side edge. A blend whose setback reaches exactly to a hole's rim (a 5 mm chamfer or
  round on a cube edge with openings 5 mm in) touches the rim at its corners and runs along its
  edge between them, so the classification came back on the face, on the edge, on the face;
  `sm_FindRolloverSideEdge` found an interval ending at a vertex, not an edge crossing, and the
  blend failed ("Unable to Process Rollover"). Intervals on the face or on its boundary are now
  one interval: only the parts of a rail outside its face roll over.
- **A full-turn revolution nests the shells it leaves inside one another**
  (`SmPrimitiveCreation.cpp`, `SmPrimitiveCreation::CreateRotationalSweep`, the full-turn
  analytic path): it builds a face per profile curve and lets `StitchAndOrient` orient them,
  which puts every closed shell in the infinite region as a solid of its own. A profile whose
  boundary meets the axis more than twice (a half annulus, both arcs ending on the axis) or
  holds a hole revolves into shells inside one another; the inner one came out as a second
  ball, the hollow ball's volume the sum of both (294604.5 for 241560.7), and SMLib's check let
  it through. Each shell found inside another's region (by a ray from one of its vertices) is
  now moved into the region of the shell directly around it, turned round where it bounds a
  void, and the regions' void flags are set from the nesting. A partial turn was right already.
- **A draft extrusion's corner pivots round its own joint** (`SmPrimitiveCreation.cpp`,
  `sm_FindPivotPoint`, used by `SmPrimitiveCreation::CreateTaperExtrude`): a corner the offset goes
  round is filled by a ruled face from its offset curve to the joint of the two profile curves it
  lies between, and the joint was the first pair of coincident ends found. A loop of two curves
  (a half disc's diameter and arc, a lens of two arcs) has two joints between the same two
  curves, so both corners took the same pivot and an outward draft came out open
  (`SMLIB_ERR_NOT_SOLID` from the shim). Of the joints, the one nearest the corner's offset curve is
  now its pivot: a half disc r20 drafted 3° outward over 10 is a valid solid of 6555.914 mm³
  (6555.908 by quadrature; the same loop with its diameter split in two gave that already).
- **A loop closed through its vertex** (`SmTrimmingTools.cpp`, `SmTrimmingTools::CheckLoop`,
  SmLoop test 0 "Edge/Edge Gap"): where two edges of a loop meet, it measured the gap between
  their lifted UV ends against the two edges' tolerances. SMLib's tolerance model (`SmTol.h`,
  3a) connects an edge end to its vertex within their XSectTol3d, so two ends in one vertex's
  neighbourhood may lie farther apart than the edges' tolerances, as its own comment noted. A
  union of two cylinders had an intersection curve end 1.3e-5 mm from its vertex (edge
  tolerance 5.8e-6, vertex 2.6e-5) and failed the check. Each end is now measured against the
  vertex point with its edge/vertex XSectTol3d.
- **A failed evaluation taken for a value** (`SmLocalSolve1d.cpp`, `SolveOscillating`, `SolveIt`;
  `SmLocalSolveNd.cpp`, `SolveIt`): an evaluator that fails leaves no function value (the curve
  property evaluator zeroes its outputs first). The 1D oscillation search skipped a failed guess
  and then read the last guess's output as a sign change, so bisection returned F = 0 as a
  minimum: the minimum radius of curvature of a straight UV trim on a cylinder came out 0. A 1D
  solve without an interval went on from a failed first evaluation. The N-D solver summed the
  equations of a failed clamped evaluation, could take a failed scaled one, and kept the first
  point's Jacobian and found flag when it took the scaled point. A sign change now counts only
  at a guess that was evaluated, a failed start without an interval returns its status, and the
  N-D step takes its equations, Jacobian and found flag from the one evaluation it uses (or
  stops when none succeeded).
- **A radius of curvature never read** (`SmBSplineCurve.cpp`, `SmBSplineCurve::AssertValid`,
  tests 13 and 14): the level-2 check ran the minimum-radius-of-curvature analysis and then
  compared the radius it had not yet read (still `SM_BIG_DOUBLE`) with the tolerance, so a curve
  with a vanishing radius of curvature was never reported and the analysis was wasted. The
  smallest radius among the analysis' solutions is now compared.
- **A torus's minor interval never checked** (`SmTorus.cpp`, `SmTorus::AssertValid`, tests 3 and 4;
  `SmSurfOfRevolution.cpp`, `SmSurfOfRevolution::AssertValid`, test 13): the test that a torus's
  V (minor) interval is at most a turn read the U interval, so a torus with a minor interval past
  a turn passed it. The U and V tests, and the surface of revolution's U test, also reported
  under the ids of other tests (2 and 3, "Same Parameterization" / "Parameterization"). The V test
  now reads the V interval, and each test reports under its own id.
- **Two curves that change sides cross; they don't coincide** (`SmCurve.cpp`,
  `SmCurve::SimpleCoincidenceChecker`, which `SmCCIGlobalSolver` runs on each pair of spans): a
  pair of spans whose end points lie within tolerance of the other curve, the curves staying within
  tolerance between, was a coincidence range between those ends. Two equal circles a little apart
  (r 5, centres 1e-4 mm apart, the walls of a prism on a ring) graze within tolerance over a
  stretch; their spans end at the circles' tops, so the range ran from one top to the other, an
  interval set by the subdivision, while the faces' own classification found the one crossing
  between them. The union then shared a 5e-5 mm edge among four faces (solidcore#480). A range is
  now a coincidence only where the separation keeps one side; where it points to opposite sides at
  the range's ends, the answer is the crossing a curve/curve intersection finds there. Curves
  that coincide (equal circles, the separation numeric zero) still give their range.
- **A boolean keeps the second operand's UV trims** (`SmMerge.cpp`, `SmMerge::MakeFaceuse`, new
  `sm_CopyUVTrimCurvesOfOtherFace`): the faces a boolean takes over from its second operand are
  rebuilt in the first by `SmBrep::MakeFace`, which was called without the operand's trims
  (`SmBrep.h`: `cpOptUVCurves`), so every such edgeuse came out with none; SMLib then made the ones
  it needed by projecting the edge (`SmEdgeuse::GetOrCreateUVTrimCurve`), and the rest stayed
  missing. On a thread's groove cut from a shaft the two long helical edges between its root and
  flanks were left without trims (tangency-labs/solidcore#474). Where the new edge runs on the
  copy of the operand edge's curve over the same interval and the face on the copy of its surface,
  the operand's trim is exact as it stands and is now copied over.
- **Curves coincide only where their separation vanishes** (`SmCurve.cpp`,
  `SmCurve::SimpleCoincidenceChecker`): after a crossing is a point (above), a pair of spans
  that stays within tolerance with the separation on one side but not zero (a tangency, two
  curves running close) is the point where they come nearest, not a range between the spans'
  ends. Contact is what the geometry gives exactly: crossings and tangencies are points, a
  coincidence range only where the curves are the same.
- **Faces that bound nothing within tolerance are removed after a boolean** (`SmMerge.cpp`,
  `sm_RemoveFacePairsBoundingNothing`, at the end of `SmMerge::ManifoldBoolean`): the
  contacts are exact (a tangency is an edge), and exact contacts within tolerance of each
  other can leave a solid thinner than tolerance. A box cut by a cylinder tangent to its face
  0.0003 from its corner edge, the edge 3e-8 off the cylinder, kept a 3e-8 sliver between a
  strip of the face and a strip of the cylinder, its caps loops of no area (solidcore#429). Two
  adjacent Faces that each lie on the other within tolerance (five points along every edge and
  an internal point within tolerance of the other's surface, inside its box) and face opposite
  ways are removed; the Edges left on one Face are stitched, and an Edge one Loop then runs
  along both ways is deleted. A representability step; it changes no contact.
- **Area counts the faces between solid and void only** (`SmBrep.cpp`,
  `SmBrep::ComputePreciseProperties`): a face whose two sides lie in different regions added its
  area, and its volume, for every faceuse it integrated, whichever of the two behaviours selected
  it. A face between two solid regions therefore counted its area once per side, although
  `SmMassPropertyBehavior` (`SmTopoTypes.h`) has the default area "only for those faces which lie
  between a solid region and a void region"; three solid regions of a merge of two overlapping
  10 mm cubes measured 1200 mm² instead of 800. Under area behaviour 2 the void side's faceuse
  also added its volume. Each faceuse now adds its volume and mass where it bounds a non void
  region, and its area and area moments where the area behaviour selects it.
- **A spindle torus's self-intersection checked without a report** (`SmTorus.cpp`,
  `SmTorus::AssertValid`, new test 5): for a torus whose major radius is below its minor one the
  check put the apex angle at 90 + acos(R/r) instead of 90 + asin(R/r), compared it with the U
  (major) interval instead of the V (minor) one, and on a hit cleared the result without a
  report. Every spindle torus with a whole turn of U failed silently, even its outer sheet alone.
  An apex may now bound the V domain but not lie inside it (as a sphere's poles, SmSphere tests 4
  and 5), and test 5 reports how far inside it lies.
- **The solvers' last dropped statuses** (`SmGlobalSolver.cpp`, `BranchMayContainAnswers`;
  `SmLocalSolve1d.cpp`, `SolveIt`; all three solver files now build with `-Werror=unused-value`):
  the branch culling for perspective projections went on with directions it could not unitize,
  and the surface-projected test with offset-box samples it could not evaluate, both giving a
  bound too small to cull by; a branch that cannot be bounded is now kept. The bisection at the
  end of a 1D solve writes nothing when it fails, so its status is discarded with that reason.
