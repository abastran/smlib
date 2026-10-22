/* SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors */
/* SPDX-License-Identifier: Apache-2.0 */

/*
 * SMLib's C API for SolidCraft: a plain C boundary in front of SM_API.
 *
 * Rules every function here follows:
 *  - Only C data crosses: scalars, caller-owned flat arrays with their lengths, and opaque
 *    handles. No C++ type, exception or SMLib pointer escapes.
 *  - Every function returns an int32_t status: SMLIB_OK, a positive SMLib status (SM_ERR_*,
 *    for example 1001 = SM_ERR, 1007 = SM_ERR_INVALID_INPUT), or a negative SMLIB_ERR_*.
 *  - Results are written through out-pointers only on SMLIB_OK.
 *  - Inputs are never changed. Operations that SM_API performs in place work on a copy.
 *  - A handle returned by this API is owned by the caller and freed with smlib_brep_free.
 *  - SMLib keeps process-global state, so calls are serialised inside the library; any thread
 *    may call, one call runs at a time.
 */
#ifndef SMLIB_C_H
#define SMLIB_C_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Bumped whenever a declaration in this header changes. */
#define SMLIB_ABI_VERSION 3

#define SMLIB_OK 0
/* A C++ exception escaped SMLib (out of memory, most likely). */
#define SMLIB_ERR_EXCEPTION (-1000)
/* An output buffer was too small; the needed length was written to *out_len. */
#define SMLIB_ERR_BUFFER (-1001)
/* An argument was rejected before SMLib saw it (NULL, not finite, out of range). */
#define SMLIB_ERR_ARG (-1002)
/* SMLib reported success but returned no body, or an empty one. */
#define SMLIB_ERR_EMPTY (-1003)

/* An SMLib B-rep, owned by the caller. */
typedef struct smlib_brep smlib_brep;

/* ---- Diagnostics ------------------------------------------------------------------------- */

/* Receives one line per API call once it returns, such as
 *   `smlib_make_box(origin=[0,0,0], size=[10,20,30]) -> 0 (0.21 ms)`.
 * `level` is 0 for a call that succeeded and 1 for one that failed. The line is valid only during
 * the callback. The callback runs with the library lock held: it must not call back into this
 * API. */
typedef void (*smlib_log_fn)(void *user, int32_t level, const char *line);

/* Sets (or, with NULL, clears) the log callback. */
void smlib_set_log(smlib_log_fn log, void *user);

/* A directory (UTF-8) where every failing call leaves what it needs to be replayed: a
 * `<n>-<function>.txt` with the call's line, a native SMLib `<n>-<function>-<input>.smb` for
 * each B-rep it was given, and `<n>-<function>-<input>.txt` (smlib_flat_read reads it) for each
 * set of flat arrays. NULL or "" turns this off (the default). */
int32_t smlib_set_dump_dir(const char *dir);

/* The call running right now, for crash reports: copies its `function(arguments)` (truncated to
 * `cap - 1` bytes, NUL-terminated) into `buf` and returns its length, 0 when no call is running.
 * Takes no lock and allocates nothing, so a crash handler may call it. */
int64_t smlib_current_call(char *buf, int64_t cap);

/* ---- Native files (for replaying dumps) --------------------------------------------------- */

/* Reads a native SMLib .smb file (binary unless `ascii`). */
int32_t smlib_brep_read(const char *path, int32_t ascii, smlib_brep **out);

/* Writes `brep` to a native SMLib .smb file (binary unless `ascii`). */
int32_t smlib_brep_write(const smlib_brep *brep, const char *path, int32_t ascii);

/* ---- B-reps ------------------------------------------------------------------------------- */

/* SMLIB_ABI_VERSION of the library actually linked. */
int32_t smlib_abi_version(void);

/* Creates SMLib's context. Optional (every call does it on first use); call it once at start-up
 * to keep the cost out of the first operation. */
int32_t smlib_init(void);

/* Frees a B-rep; NULL is ignored. */
void smlib_brep_free(smlib_brep *brep);

/* A deep copy of `brep`. */
int32_t smlib_brep_copy(const smlib_brep *brep, smlib_brep **out);

/* `brep` moved by the affine `matrix`: rows 0-2 the images of the x, y and z axes, row 3 the
 * translation (column 3 is ignored), as SmBrep::Transform takes it. Reflections mirror the body;
 * a non-uniform scale or a shear turns its geometry into NURBS; tolerances scale with it. */
int32_t smlib_brep_transform(const smlib_brep *brep, const double matrix[4][4], smlib_brep **out);

/* An axis-aligned box from its lower corner `origin` with edge lengths `size`, in mm. */
int32_t smlib_make_box(const double origin[3], const double size[3], smlib_brep **out);

/* The enclosed volume at `relative_accuracy` (clamped by SMLib to [1e-8, 1e-1]). */
int32_t smlib_brep_volume(const smlib_brep *brep, double relative_accuracy, double *out);

/* Face, edge and vertex counts, in that order. */
int32_t smlib_brep_counts(const smlib_brep *brep, int64_t out[3]);

/* ---- Primitives (src/smlib_primitives.cpp) ----------------------------------------------- */

/* A cylinder about +Z: its base circle's centre `origin`, radius and height (mm, > 0). */
int32_t smlib_make_cylinder(const double origin[3], double radius, double height, smlib_brep **out);

/* A sphere from its centre and radius (mm, > 0). */
int32_t smlib_make_sphere(const double center[3], double radius, smlib_brep **out);

/* A cone or frustum about +Z: its base centre `origin`, base radius (> 0), top radius (>= 0,
 * 0 for a pointed cone, different from the base radius) and height (> 0), in mm. */
int32_t smlib_make_cone(const double origin[3], double base_radius, double top_radius, double height, smlib_brep **out);

/* A ring torus about +Z from its centre, with major radius > minor radius > 0 (mm). */
int32_t smlib_make_torus(const double center[3], double major_radius, double minor_radius, smlib_brep **out);

/* ---- History (src/smlib_history.cpp) ------------------------------------------------------ */

/* Face and edge history out of an operation: tag each operand, run the operation on the tagged
 * copies, read the result's history. SMLib carries the tags itself through every split, merge,
 * copy and coincidence the operation makes, so any operation works, and the history is what the
 * operation did, never a guess from geometry.
 *
 * Indices are positions in a body's face and edge order (the order SMLib lists them in, which
 * smlib_brep_copy keeps): operand indices in the tagged copy's, result indices in the result's. */

#define SMLIB_HISTORY_VERSION 1

/* How a result face relates to its origins (solidcraft_geom::brep::history::Change). */
#define SMLIB_CHANGE_KEPT 0      /* one origin, still that whole face, boundary unchanged */
#define SMLIB_CHANGE_MODIFIED 1  /* one origin, trimmed differently */
#define SMLIB_CHANGE_SPLIT 2     /* one origin, which other result faces also came from */
#define SMLIB_CHANGE_MERGED 3    /* several origins (coincident or joined faces) */
#define SMLIB_CHANGE_GENERATED 4 /* no origin: new to the operation */

/* A copy of `brep` whose face i and edge i are tagged as face i and edge i of operand `operand`
 * (0 <= operand < 2^30). Tags a tagged body afresh. */
int32_t smlib_history_tag(const smlib_brep *brep, int32_t operand, smlib_brep **out);

/* The history of `result`, as int64 values in `buf` (`cap` of them):
 *   SMLIB_HISTORY_VERSION, face count F, edge count E, then
 *   F faces, in face order:  change, n, then n (operand, face) origins, sorted;
 *   E edges, in edge order:  m, then the m result faces it bounds (ascending),
 *                            n, then n (operand, edge) origins, sorted.
 * An edge with no origin (n = 0) is new: it lies where its faces meet (a boolean's intersection
 * curve). Untagged faces and edges have no origin. Writes the length to *out_len; when it is more
 * than `cap`, writes nothing else and returns SMLIB_ERR_BUFFER (call again with a buffer that
 * big; `buf` may be NULL when `cap` is 0). */
int32_t smlib_history_read(const smlib_brep *result, int64_t *buf, int64_t cap, int64_t *out_len);

/* ---- Validity (smlib_check.cpp) ----------------------------------------------------------- */

/* Runs SMLib's own validity checks on `brep`: SmBrep::ValidatePointers, then, only when the
 * pointers are sound, AssertValid over the whole topology graph (SM_WALK). Both run at test
 * `level`: 0 the fewest tests, 1 adds the geometry checks, 2 every test (the slow ones too).
 *
 * The number of failed tests goes to *out_failures, 0 meaning SMLib finds the body valid. The
 * report, one line per failed test, is written NUL-terminated to `buf` and its length (without
 * the NUL) to *out_len. A failing check is still SMLIB_OK; the status only says the check ran.
 * When `cap` is too small, nothing goes to `buf` or *out_failures: the call returns
 * SMLIB_ERR_BUFFER with the needed length in *out_len (`buf` may be NULL with `cap` 0). */
int32_t smlib_brep_check(const smlib_brep *brep, int32_t level, char *buf, int64_t cap, int64_t *out_len,
                         int64_t *out_failures);

/* ---- Mass properties (smlib_measure.cpp) ------------------------------------------------- */

/* Exact mass properties of a closed solid at unit density, integrated by SMLib on the exact
 * surfaces (no tessellation) at `relative_accuracy`, clamped to [1e-8, 1e-1]. At 1e-8 analytic
 * primitives come out within 2e-7 of their closed forms. Into `out`, in mm units:
 *   [0] area, [1] volume, [2..5) centroid x, y, z,
 *   [5..8)  moments of inertia Ixx, Iyy, Izz about axes through the centroid, parallel to the
 *           world axes: the integrals of (y^2+z^2, z^2+x^2, x^2+y^2) dV,
 *   [8..11) products of inertia Iyz, Izx, Ixy about the centroid: the integrals of (yz, zx, xy) dV
 *           relative to it. The inertia tensor's off-diagonal terms are their negatives.
 * For a density rho, the mass is rho * volume and the inertia scales by rho. A body that is not a closed manifold solid
 * (a sheet, an open shell) is rejected with SM_ERR_INVALID_INPUT; smlib_brep_area measures
 * those. */
int32_t smlib_brep_mass_properties(const smlib_brep *brep, double relative_accuracy, double out[11]);

/* Total area of every face, each counted once, so sheets and open shells too, at
 * `relative_accuracy`, clamped to [1e-8, 1e-1]. */
int32_t smlib_brep_area(const smlib_brep *brep, double relative_accuracy, double *out);

/* The exact axis-aligned box of `brep`: the extrema of its vertices, edges and faces on their
 * exact curves and surfaces (SmBrep::CalculateTightBoundingBox), not grown by tolerances or
 * control hulls. Into `lo` the least x, y, z, into `hi` the greatest, in mm. A body with no
 * geometry (an empty box) is SM_ERR_INVALID_INPUT. */
int32_t smlib_brep_bounds(const smlib_brep *brep, double lo[3], double hi[3]);

/* SMLib's own check that `brep` is one closed, manifold solid: *out = 1 if so, else 0. */
int32_t smlib_brep_is_manifold_solid(const smlib_brep *brep, int32_t *out);

/* ---- Flat B-rep arrays -------------------------------------------------------------------
 *
 * One B-rep as plain arrays, laid out like one member of USD's BrepArray schema (the arrays of
 * source/BREP_USD_DATA/inc/UsdBrepArrayData.h) with the same units, so a caller can map an SMLib
 * B-rep onto its own topology and back without USD. Every array is either int64 ("i64") or
 * double ("f64"); 3-vectors are flattened x,y,z. Indices are 0-based and local to this B-rep.
 * Lengths are in model units (mm for SolidCraft); analytic angles and analytic face ranges
 * are in radians, as in USD.
 *
 * Topology (as in SMLib and USD):
 *  - Region 0 is the infinite region. Regions own consecutive blocks of shells, shells own
 *    consecutive blocks of faceuses, faces own consecutive blocks of loops, loops own
 *    consecutive blocks of edgeuses, all in array order.
 *  - Every region but the infinite one has one shell that encloses it, the region's others
 *    lying inside that one. Export writes it first; import takes a region's shells in any order
 *    and puts first the one SMLib's containment test (ray casts) finds enclosing the others.
 *  - Each face has exactly two faceuses (one per side), so #faceuses = 2 * #faces.
 *  - An edgeuse record is one mated pair of SMLib edgeuses: the use of an edge by one loop
 *    (seen from the face's "same" side) and its mate on the other side.
 *  - Edges here are face edges; wire edges (shells made of one bare edge) are separate.
 *    Vertices here are edge and loop vertices; point shells keep their point separately.
 *
 * Geometry is packed per kind, in order of appearance: the k-th face of type PLANE uses entry k
 * of the FACE_PLANE_* arrays, and so on. NURBS control points are given unweighted (cartesian)
 * with their weights alongside (a weights array may also be empty, meaning all 1); knot vectors
 * are full ("flat"): #knots = #points + order.
 *
 * What is not carried: UV trim curves of faces on analytic surfaces (SMLib's are expressed in
 * its internal NURBS form of the analytic surface, which is not the USD parameterisation), so
 * SMLib rebuilds them from the 3D edges on import. Edge curves keep their SMLib class where it
 * is a line, circle or ellipse (SmLine, SmCircle, SmEllipse); B-splines stay NURBS, and other
 * curve kinds are approximated as NURBS to the B-rep's tolerance; surfaces that are neither
 * analytic nor NURBS are approximated likewise. A line, circle or ellipse edge comes back with
 * the parameterisation its kind defines (a line's by the distance along its direction vector,
 * a circle's or ellipse's by angle in degrees, as SMLib has it), not the original's; on a
 * NURBS face, its UV trims are re-mapped to match, or dropped (and rebuilt by SMLib) when the
 * two parameterisations are not affinely related. Attributes, materials and colours are not
 * carried.
 *
 * Length rules below: n_x is the length of the array named. "per X" means one entry per X. */

/* Opaque, caller-owned set of flat arrays for one B-rep. */
typedef struct smlib_flat smlib_flat;

/* --- int64 arrays --- */
/* Shells per region. Length n_regions (>= 1). */
#define SMLIB_FLAT_REGION_SHELL_COUNT 1
/* Per region: 0 = solid region, 1 = void region. Length n_regions. */
#define SMLIB_FLAT_REGION_TYPE 2
/* Faceuses per shell. Length n_shells = sum(REGION_SHELL_COUNT). */
#define SMLIB_FLAT_SHELL_FACEUSE_COUNT 3
/* Wire edges per shell (0 or 1). Length n_shells. */
#define SMLIB_FLAT_SHELL_WIRE_EDGE_COUNT 4
/* Per shell: 1 = a point shell (no faceuses, no wire edges, one point in
 * SHELL_POINT_POSITION), 0 = otherwise. Length n_shells. */
#define SMLIB_FLAT_SHELL_POINT_TYPE 5
/* Face index of each faceuse. Length n_faceuses = sum(SHELL_FACEUSE_COUNT) = 2 * n_faces. */
#define SMLIB_FLAT_FACEUSE_FACE 6
/* Per faceuse: 1 = the side the surface normal points to ("same"), 0 = the other side.
 * Length n_faceuses. */
#define SMLIB_FLAT_FACEUSE_ORIENTATION 7
/* Loops per face (>= 1; the first is the outer loop). Length n_faces. */
#define SMLIB_FLAT_FACE_LOOP_COUNT 8
/* Per face, its surface kind: one of SMLIB_FLAT_SURFACE_*. Length n_faces. */
#define SMLIB_FLAT_FACE_SURFACE_TYPE 9
/* Per face: 1 = "rectangular" (trimmed only by its UV range), 0 = "general". Length n_faces. */
#define SMLIB_FLAT_FACE_TRIM_TYPE 10
/* Edgeuses per loop; 0 = a vertex loop (a single vertex, see LOOP_VERTEX).
 * Length n_loops = sum(FACE_LOOP_COUNT). */
#define SMLIB_FLAT_LOOP_EDGEUSE_COUNT 11
/* Per loop: the vertex of a vertex loop, -1 for an edgeuse loop. Length n_loops. */
#define SMLIB_FLAT_LOOP_VERTEX 12
/* Edge of each edgeuse. Length n_edgeuses = sum(LOOP_EDGEUSE_COUNT). */
#define SMLIB_FLAT_EDGEUSE_EDGE 13
/* Per edgeuse: 1 = the loop runs along the edge's curve direction ("same"), 0 = against it.
 * Length n_edgeuses. */
#define SMLIB_FLAT_EDGEUSE_ORIENTATION 14
/* Per edgeuse: the next edgeuse of the same edge in a right-hand traversal around the edge
 * (the radial list; it must cycle through all of the edge's edgeuses). Length n_edgeuses. */
#define SMLIB_FLAT_EDGEUSE_NEXT_RADIAL 15
/* Per edgeuse: 1 = "topEntry" (the traversal around the edge meets this pair's top use
 * first), 0 = "bottomEntry". Length n_edgeuses. */
#define SMLIB_FLAT_EDGEUSE_RADIAL_ENTRY 16
/* Per edge, its curve kind: one of SMLIB_FLAT_CURVE_*. Length n_edges. */
#define SMLIB_FLAT_EDGE_CURVE_TYPE 17
/* Start and end vertex of each edge. Length 2 * n_edges. */
#define SMLIB_FLAT_EDGE_VERTICES 18
/* Per wire edge, its curve kind. Length n_wire_edges = sum(SHELL_WIRE_EDGE_COUNT). */
#define SMLIB_FLAT_WIRE_EDGE_CURVE_TYPE 19
/* Start and end vertex of each wire edge. Length 2 * n_wire_edges. */
#define SMLIB_FLAT_WIRE_EDGE_VERTICES 20
/* Per vertex: 1 = it has a point in VERTEX_POSITION (always 1 here; USD's "BrepPointAPI").
 * Length n_vertices; every vertex must be used by an edge, a wire edge or a vertex loop. */
#define SMLIB_FLAT_VERTEX_POINT_TYPE 21
/* Per NURBS edge: control point count (>= order) and order (degree + 1, in [2, 18]).
 * Length n_nurbs_edges each. */
#define SMLIB_FLAT_EDGE_NURBS_VERTEX_COUNT 22
#define SMLIB_FLAT_EDGE_NURBS_ORDER 23
/* The same for NURBS wire edges. */
#define SMLIB_FLAT_WIRE_EDGE_NURBS_VERTEX_COUNT 24
#define SMLIB_FLAT_WIRE_EDGE_NURBS_ORDER 25
/* Per edgeuse: control point count and order of its UV trim curve, both 0 when it has none.
 * Length n_edgeuses, or 0 when no edgeuse has one. Only faces of type NURBS carry trims. */
#define SMLIB_FLAT_EDGEUSE_NURBS_VERTEX_COUNT 26
#define SMLIB_FLAT_EDGEUSE_NURBS_ORDER 27
/* Per NURBS face: control point counts and orders in U and V. Length n_nurbs_faces each. */
#define SMLIB_FLAT_FACE_NURBS_U_VERTEX_COUNT 28
#define SMLIB_FLAT_FACE_NURBS_V_VERTEX_COUNT 29
#define SMLIB_FLAT_FACE_NURBS_U_ORDER 30
#define SMLIB_FLAT_FACE_NURBS_V_ORDER 31
/* 1 = written by SMLib (import trusts it and skips the healer), 0 = from elsewhere (import
 * heals it). Length 1. */
#define SMLIB_FLAT_SOURCE 32

/* --- double arrays --- */
/* The B-rep's intersection tolerance (twice its largest face/edge/vertex tolerance). Length 1. */
#define SMLIB_FLAT_TOLERANCE 101
/* Bounding box: min x,y,z then max x,y,z (informational; import ignores it). Length 6. */
#define SMLIB_FLAT_EXTENT 102
/* Per face: u_min, v_min, u_max, v_max. NURBS: the face's native UV domain. Analytic: the
 * STEP parameterisation, angles in radians: plane (u, v) distances along the placement's
 * x and y axes; sphere (longitude, latitude); cylinder and cone (angle, distance along the
 * axis); torus (major angle, minor angle). Length 4 * n_faces. */
#define SMLIB_FLAT_FACE_RANGE 103
/* Per edge: start and end curve parameter (radians for circle and ellipse curves).
 * Length 2 * n_edges. */
#define SMLIB_FLAT_EDGE_RANGE 104
/* Per wire edge: the same. Length 2 * n_wire_edges. */
#define SMLIB_FLAT_WIRE_EDGE_RANGE 105
/* Point of each point shell, x,y,z. Length 3 * (number of point shells). */
#define SMLIB_FLAT_SHELL_POINT_POSITION 106
/* Point of each vertex, x,y,z. Length 3 * n_vertices. */
#define SMLIB_FLAT_VERTEX_POSITION 107
/* NURBS edge curves: control points x,y,z (3 * sum(count)), knots (sum(count + order)),
 * weights (> 0, sum(count)). */
#define SMLIB_FLAT_EDGE_NURBS_CVS 108
#define SMLIB_FLAT_EDGE_NURBS_KNOTS 109
#define SMLIB_FLAT_EDGE_NURBS_WEIGHTS 110
/* Circle edge curves: center, axis (normal), reference direction (x axis, angle 0), 3 each per
 * circle; radius, 1 per circle. */
#define SMLIB_FLAT_EDGE_CIRCLE_CENTER 111
#define SMLIB_FLAT_EDGE_CIRCLE_AXIS 112
#define SMLIB_FLAT_EDGE_CIRCLE_REF_DIRECTION 113
#define SMLIB_FLAT_EDGE_CIRCLE_RADIUS 114
/* Line edge curves: origin and direction, 3 each per line (point = origin + t * direction). */
#define SMLIB_FLAT_EDGE_LINE_ORIGIN 115
#define SMLIB_FLAT_EDGE_LINE_DIRECTION 116
/* Ellipse edge curves: center, axis, reference direction (3 each); x and y radius (1 each). */
#define SMLIB_FLAT_EDGE_ELLIPSE_CENTER 117
#define SMLIB_FLAT_EDGE_ELLIPSE_AXIS 118
#define SMLIB_FLAT_EDGE_ELLIPSE_REF_DIRECTION 119
#define SMLIB_FLAT_EDGE_ELLIPSE_X_RADIUS 120
#define SMLIB_FLAT_EDGE_ELLIPSE_Y_RADIUS 121
/* The same fourteen arrays for wire edges. */
#define SMLIB_FLAT_WIRE_EDGE_NURBS_CVS 122
#define SMLIB_FLAT_WIRE_EDGE_NURBS_KNOTS 123
#define SMLIB_FLAT_WIRE_EDGE_NURBS_WEIGHTS 124
#define SMLIB_FLAT_WIRE_EDGE_CIRCLE_CENTER 125
#define SMLIB_FLAT_WIRE_EDGE_CIRCLE_AXIS 126
#define SMLIB_FLAT_WIRE_EDGE_CIRCLE_REF_DIRECTION 127
#define SMLIB_FLAT_WIRE_EDGE_CIRCLE_RADIUS 128
#define SMLIB_FLAT_WIRE_EDGE_LINE_ORIGIN 129
#define SMLIB_FLAT_WIRE_EDGE_LINE_DIRECTION 130
#define SMLIB_FLAT_WIRE_EDGE_ELLIPSE_CENTER 131
#define SMLIB_FLAT_WIRE_EDGE_ELLIPSE_AXIS 132
#define SMLIB_FLAT_WIRE_EDGE_ELLIPSE_REF_DIRECTION 133
#define SMLIB_FLAT_WIRE_EDGE_ELLIPSE_X_RADIUS 134
#define SMLIB_FLAT_WIRE_EDGE_ELLIPSE_Y_RADIUS 135
/* UV trim curves, packed in edgeuse order for edgeuses with a non-zero count: control points
 * u,v (2 * sum(count)), knots (sum(count + order)), weights (sum(count)). Each curve runs in
 * the direction of its edgeuse; its knot interval may be any interval and is mapped onto the
 * edge's range on import. */
#define SMLIB_FLAT_EDGEUSE_NURBS_CVS 136
#define SMLIB_FLAT_EDGEUSE_NURBS_KNOTS 137
#define SMLIB_FLAT_EDGEUSE_NURBS_WEIGHTS 138
/* NURBS surfaces: control points x,y,z, U-major (point (i, j) at index i * v_count + j),
 * 3 * sum(u_count * v_count); U knots sum(u_count + u_order); V knots sum(v_count + v_order);
 * weights (> 0) sum(u_count * v_count). */
#define SMLIB_FLAT_FACE_NURBS_CVS 139
#define SMLIB_FLAT_FACE_NURBS_U_KNOTS 140
#define SMLIB_FLAT_FACE_NURBS_V_KNOTS 141
#define SMLIB_FLAT_FACE_NURBS_WEIGHTS 142
/* Analytic surfaces. Each has an origin (center), an axis (z axis: the normal of a plane, the
 * axis of revolution otherwise) and a reference direction (x axis), 3 each per surface; the
 * y axis is axis x reference direction. Then: */
#define SMLIB_FLAT_FACE_SPHERE_CENTER 143
#define SMLIB_FLAT_FACE_SPHERE_AXIS 144
#define SMLIB_FLAT_FACE_SPHERE_REF_DIRECTION 145
#define SMLIB_FLAT_FACE_SPHERE_RADIUS 146 /* 1 per sphere */
#define SMLIB_FLAT_FACE_PLANE_ORIGIN 147
#define SMLIB_FLAT_FACE_PLANE_AXIS 148
#define SMLIB_FLAT_FACE_PLANE_REF_DIRECTION 149
#define SMLIB_FLAT_FACE_CYLINDER_ORIGIN 150
#define SMLIB_FLAT_FACE_CYLINDER_AXIS 151
#define SMLIB_FLAT_FACE_CYLINDER_REF_DIRECTION 152
#define SMLIB_FLAT_FACE_CYLINDER_RADIUS 153 /* 1 per cylinder */
#define SMLIB_FLAT_FACE_CONE_ORIGIN 154
#define SMLIB_FLAT_FACE_CONE_AXIS 155
#define SMLIB_FLAT_FACE_CONE_REF_DIRECTION 156
#define SMLIB_FLAT_FACE_CONE_RADIUS 157 /* 1 per cone: the radius at the origin (v = 0) */
#define SMLIB_FLAT_FACE_CONE_SEMI_ANGLE 158 /* 1 per cone, radians, in (-pi/2, pi/2) */
#define SMLIB_FLAT_FACE_TORUS_ORIGIN 159
#define SMLIB_FLAT_FACE_TORUS_AXIS 160
#define SMLIB_FLAT_FACE_TORUS_REF_DIRECTION 161
#define SMLIB_FLAT_FACE_TORUS_MAJOR_RADIUS 162 /* 1 per torus */
#define SMLIB_FLAT_FACE_TORUS_MINOR_RADIUS 163 /* 1 per torus */
/* Per face edge (as EDGE_CURVE_TYPE) and per vertex (as VERTEX_POINT_TYPE): its own tolerance,
 * in TOLERANCE's terms (twice the edge's or vertex's zone tolerance). Export writes both; import
 * reads them when given (each finite, > 0, at most 1), else every edge and vertex takes
 * TOLERANCE. An approximated edge (a fillet's rail) may need more than the body's tolerance, and
 * only it then carries more. */
#define SMLIB_FLAT_EDGE_TOLERANCE 164
#define SMLIB_FLAT_VERTEX_TOLERANCE 165

/* FACE_SURFACE_TYPE values. */
#define SMLIB_FLAT_SURFACE_NURBS 0
#define SMLIB_FLAT_SURFACE_PLANE 1
#define SMLIB_FLAT_SURFACE_CYLINDER 2
#define SMLIB_FLAT_SURFACE_CONE 3
#define SMLIB_FLAT_SURFACE_SPHERE 4
#define SMLIB_FLAT_SURFACE_TORUS 5

/* EDGE_CURVE_TYPE and WIRE_EDGE_CURVE_TYPE values. Export writes LINE, CIRCLE and ELLIPSE for
 * SMLib's SmLine, SmCircle and SmEllipse curves and NURBS for everything else. */
#define SMLIB_FLAT_CURVE_NURBS 0
#define SMLIB_FLAT_CURVE_LINE 1
#define SMLIB_FLAT_CURVE_CIRCLE 2
#define SMLIB_FLAT_CURVE_ELLIPSE 3

/* An array longer than this many elements is rejected. */
#define SMLIB_FLAT_MAX_LEN ((int64_t)1 << 27)

/* A new, empty set of arrays (NULL when out of memory). */
smlib_flat *smlib_flat_new(void);

/* Frees a set of arrays; NULL is ignored. */
void smlib_flat_free(smlib_flat *flat);

/* The arrays of `brep`, as a new set owned by the caller. */
int32_t smlib_brep_export(const smlib_brep *brep, smlib_flat **out);

/* A new B-rep built from `flat`, which is validated first: any malformed count, index, length,
 * kind or value returns SM_ERR_INVALID_INPUT (1007) and builds nothing. */
int32_t smlib_brep_import(const smlib_flat *flat, smlib_brep **out);

/* Read access to array `id`: *data points at *len elements, valid until the array is set again
 * or `flat` is freed (*data may be NULL when *len is 0). An unknown id, or an id of the other
 * element type, returns SMLIB_ERR_ARG. */
int32_t smlib_flat_get_i64(const smlib_flat *flat, int32_t id, const int64_t **data, int64_t *len);
int32_t smlib_flat_get_f64(const smlib_flat *flat, int32_t id, const double **data, int64_t *len);

/* Replaces array `id` with a copy of `len` elements at `data` (`data` may be NULL when `len` is
 * 0). An unknown id, the other element type, a negative `len` or one over SMLIB_FLAT_MAX_LEN
 * returns SMLIB_ERR_ARG. */
int32_t smlib_flat_set_i64(smlib_flat *flat, int32_t id, const int64_t *data, int64_t len);
int32_t smlib_flat_set_f64(smlib_flat *flat, int32_t id, const double *data, int64_t len);

/* Writes `flat` to a text file (one line per non-empty array: name, type, length, values),
 * the format of smlib_brep_import's failure dumps, and reads one back. */
int32_t smlib_flat_write(const smlib_flat *flat, const char *path);
int32_t smlib_flat_read(const char *path, smlib_flat **out);

/* ---- Operations (src/smlib_ops.cpp) ----------------------------------------------------- */

/* SMLib reported success, but a solid input gave a result that is not a closed manifold solid. */
#define SMLIB_ERR_NOT_SOLID (-1100)
/* A selection point is farther than the tolerance from every edge (face). */
#define SMLIB_ERR_PICK_NONE (-1101)
/* A selection point is as near to two edges (faces) as the body's modelling tolerance can tell. */
#define SMLIB_ERR_PICK_AMBIGUOUS (-1102)
/* An inward offset leaves no room inside the body: no point of it lies farther than the wall
 * from every face that is offset (smlib_shell, smlib_offset_faces). */
#define SMLIB_ERR_NO_ROOM (-1103)
/* SMLib reported success but the operation did nothing: a blend (smlib_fillet, smlib_chamfer,
 * smlib_blend, smlib_fillet_variable, smlib_chamfer_sides) that made no face, a face move or draft (smlib_offset_faces, smlib_draft) that left no face on the
 * surface it was moved onto. The call's log line names it. */
#define SMLIB_ERR_NO_EFFECT (-1107)

/* Selecting edges and faces. Operations that act on some edges or faces take them as points
 * (`points`: x, y, z for each of `n_points`, at most 4096). Each point selects the edge (face)
 * nearest to it; it must lie within `tolerance` (> 0) of it, and every other edge (face) within
 * `tolerance` must be farther by more than the body's modelling tolerance. A point that selects
 * nothing fails the call with SMLIB_ERR_PICK_NONE, one that cannot decide with
 * SMLIB_ERR_PICK_AMBIGUOUS; either way, when `out_failed_point` is not NULL, the index of that
 * point is written to it (the one exception to "outputs only on SMLIB_OK"). Points that select
 * the same edge (face) count once. */

#define SMLIB_BOOLEAN_UNION 0
/* a minus b. */
#define SMLIB_BOOLEAN_SUBTRACT 1
#define SMLIB_BOOLEAN_INTERSECT 2

/* A boolean's result may be closed but not manifold: solids touching each other along an edge
 * (an edge four or more faces share, in pairs) or at a vertex, as cuts that leave pieces meeting
 * along a line give. Without this flag such a result is SMLIB_ERR_NOT_SOLID. */
#define SMLIB_RESULT_NONMANIFOLD 1

/* A new body: `a` and `b` combined by `op` (SMLIB_BOOLEAN_*); both inputs are left unchanged.
 * A result that is the empty set (an intersection of disjoint bodies, a subtraction of a body
 * that encloses `a`) is SMLIB_ERR_EMPTY. A result of several separate solids (a union of
 * disjoint bodies) is one body: smlib_brep_solid_count tells how many. Solid inputs give a
 * closed result: manifold, or (SMLIB_RESULT_NONMANIFOLD in `flags`, which smlib_boolean always
 * passes) closed with solids touching along an edge or at a vertex; anything else, a sheet or
 * an edge bounding a face on one side only, is SMLIB_ERR_NOT_SOLID. Any other bit in `flags` is
 * SMLIB_ERR_ARG. */
int32_t smlib_boolean(const smlib_brep *a, const smlib_brep *b, int32_t op, smlib_brep **out);
int32_t smlib_boolean_with(const smlib_brep *a, const smlib_brep *b, int32_t op, int32_t flags, smlib_brep **out);

/* A new sheet: `a` and `b`, sheets of planar faces on one plane (smlib_planar_faces), combined as
 * regions of that plane by `op` (SMLIB_BOOLEAN_*: union, `a` less `b`, intersection), with
 * SM_API's 2D boolean (SmApiBoolean2d). The result's edges are the inputs' own curves, split where
 * they cross; both inputs are left unchanged. Faces off the plane of the others are
 * SM_ERR_INVALID_INPUT; an input or a result without faces is SMLIB_ERR_EMPTY. */
int32_t smlib_boolean2d(const smlib_brep *a, const smlib_brep *b, int32_t op, smlib_brep **out);

/* The number of separate solid regions of `brep` (voids not counted). */
int32_t smlib_brep_solid_count(const smlib_brep *brep, int64_t *out);

/* A new body: `brep` with the selected edges rounded at constant `radius` (> 0) by a circular
 * cross-section. Each selected edge brings its whole tangent chain: every edge that continues
 * a rounded edge within 0.5 degrees at a shared vertex is rounded too. */
int32_t smlib_fillet(const smlib_brep *brep, const double *points, int64_t n_points, double tolerance, double radius,
                     smlib_brep **out, int64_t *out_failed_point);

/* A new body: `brep` with the selected edges (and their tangent chains, as for smlib_fillet)
 * cut by an equal-distance chamfer, `distance` (> 0) back from the edge along each face, with
 * bevelled corners. The setback is measured at each edge's midpoint; SMLib keeps the chord
 * across the chamfer constant along an edge, so on an edge whose dihedral angle varies the
 * setback varies with it. */
int32_t smlib_chamfer(const smlib_brep *brep, const double *points, int64_t n_points, double tolerance, double distance,
                      smlib_brep **out, int64_t *out_failed_point);

/* An operation SMLib cannot build exactly yet; the call's log line says which. */
#define SMLIB_ERR_UNSUPPORTED (-1105)

/* What smlib_blend makes of each selected edge. */
/* Rounded by a circular cross-section at constant radius (as smlib_fillet). */
#define SMLIB_BLEND_ROUND 0
/* Rounded by a curvature-continuous (G2) blend-curve cross-section between the contact lines
 * the round of that radius has: the quintic whose first and last three control points lie on
 * the contacts' tangent lines, the inner ones a third and two thirds of the way from each
 * contact to where the two tangent lines meet (on planar faces, zero curvature at the contacts). */
#define SMLIB_BLEND_CURVATURE 1
/* Rounded by a circular cross-section whose width from contact line to contact line (the
 * chord) is constant: the size is that width. */
#define SMLIB_BLEND_CHORD 2
/* Cut by an equal-distance chamfer (as smlib_chamfer): the size is the setback. */
#define SMLIB_BLEND_CHAMFER 3

/* A new body: `brep` with the selected edges blended as `kind` (SMLIB_BLEND_*), each at its
 * point's size (`sizes`, one per point, each > 0). Each selected edge brings its whole tangent
 * chain, as for smlib_fillet; an edge reached from several points takes the first point's size.
 * smlib_fillet and smlib_chamfer are this call with one size for every point. */
int32_t smlib_blend(const smlib_brep *brep, int32_t kind, const double *points, const double *sizes, int64_t n_points,
                    double tolerance, smlib_brep **out, int64_t *out_failed_point);

/* A new body: `brep` with the one edge `point` selects (no tangent chain) rounded by a circular
 * cross-section whose radius varies along it: a clamped B-spline of `degree` (1 to 3) over the
 * position along the edge, 0 at its end nearest `start` and 1 at the other, with `n_knots`
 * knots (n_radii + degree + 1, rising from 0 to 1, each end repeated degree + 1 times) and
 * `n_radii` (2 to 100) control radii, each >= 0. The position is the edge curve's parameter,
 * scaled to 0..1, which is the fraction of its length on a line or a circular arc. A point that
 * selects nothing or cannot decide fails as for any selection (no index is written). Where the
 * law has a corner inside the edge (its slope jumps where a knot repeats `degree` times) the
 * round turns, and is split into faces there. A law whose radius is zero inside the edge (the
 * round would pinch to a point) is SMLIB_ERR_UNSUPPORTED. */
int32_t smlib_fillet_variable(const smlib_brep *brep, const double *point, double tolerance, const double *start, int32_t degree,
                              const double *knots, int64_t n_knots, const double *radii, int64_t n_radii, smlib_brep **out);

/* A new body: `brep` with the selected edges and their tangent chains (as for smlib_fillet) cut
 * by an unequal chamfer with bevelled corners: `distance` (> 0) back from the edge along the face
 * whose outward normal at the edge's midpoint leans most along `first` (non-zero), and along the
 * other face `second` (> 0) back, or (`second_is_angle` 1) at the angle `second` (radians, in
 * (0, pi)) the chamfer makes with the first face. Setbacks are measured at each edge's
 * midpoint, as for smlib_chamfer. */
int32_t smlib_chamfer_sides(const smlib_brep *brep, const double *points, int64_t n_points, double tolerance, double distance,
                            double second, int32_t second_is_angle, const double *first, smlib_brep **out, int64_t *out_failed_point);

/* A new body: `brep` hollowed to a wall of `thickness` (non-zero; negative: the wall lies
 * inside the original surface, positive: outside), open at the selected faces. `n_points` may
 * be 0 (and `points` NULL) for a closed hollow body. An inward wall that leaves no room inside
 * the solid (no point of it farther than the wall from every face but the open ones, by a
 * margin of 0.1% of the wall or ten times the body's tolerance) is SMLIB_ERR_NO_ROOM. */
int32_t smlib_shell(const smlib_brep *brep, const double *points, int64_t n_points, double tolerance, double thickness,
                    smlib_brep **out, int64_t *out_failed_point);

/* A new body: the solid `brep` with the selected faces moved `distance` (non-zero; negative:
 * into the body) along their outward normals, each onto its surface's offset, and the faces
 * around them extended or trimmed to meet the moved faces (SMLib's local operation). The
 * topology is kept where the moved faces meet their neighbours as the old ones did. */
int32_t smlib_offset_faces(const smlib_brep *brep, const double *points, int64_t n_points, double tolerance, double distance,
                           smlib_brep **out, int64_t *out_failed_point);

/* The faces of `brep` under the points tilted by `angle` radians (0 < |angle| < 89.94 degrees)
 * about the line where each meets the neutral plane (`neutral_origin`, `neutral_normal`), so
 * that a positive angle leans them in going along `pull`; the faces around them extended or
 * trimmed to meet them (SMLib's local operation, as smlib_offset_faces). Each face must be flat
 * and parallel to `pull`: SMLIB_ERR_UNSUPPORTED otherwise. The hinge is where a line along
 * `pull` from a point inside the face meets the neutral plane. */
int32_t smlib_draft(const smlib_brep *brep, const double *points, int64_t n_points, double tolerance, const double neutral_origin[3],
                    const double neutral_normal[3], const double pull[3], double angle, smlib_brep **out, int64_t *out_failed_point);

/* A new solid: the sheet `brep` (not a closed solid: SMLIB_ERR_ARG) thickened by `thickness`
 * (non-zero) along its surfaces' normals (negative: against them), or by half of it each way
 * when `symmetric` is non-zero. The sheet and its offset are joined by ruled walls along the
 * sheet's open edges. */
int32_t smlib_thicken(const smlib_brep *brep, double thickness, int32_t symmetric, smlib_brep **out);

/* A new body: `brep` with its topological edges and vertices removed (SMLib's
 * SmBrep::RemoveTopologicalEdgesAndVertices): an edge between two faces on one surface (or on
 * surfaces that coincide over the other face) goes, its faces becoming one, and a vertex where
 * only two edges on one curve meet goes, its edges becoming one. The edges under `keep` (each
 * point takes the nearest edge, as the operations above pick them; none when `n_keep` is 0) stay,
 * so faces split along them on purpose stay apart. The shape is unchanged; a solid stays a
 * closed manifold solid (SMLIB_ERR_NOT_SOLID otherwise). */
int32_t smlib_merge_faces(const smlib_brep *brep, const double *keep, int64_t n_keep, double tolerance, smlib_brep **out,
                          int64_t *out_failed_point);

/* ---- Surface/surface intersection (src/smlib_ssi.cpp) ------------------------------------
 *
 * Two surfaces described in plain C, intersected by one of SMLib's SSI entry points, with each
 * intersection curve sampled into caller buffers.
 *
 * The constructors take exactly the fields of an SSI corpus case file (ssi-corpus,
 * cases/data/<id>.json, built as its runner/ssi_runner.cpp builds them), so a case replays
 * through this API unchanged. Angles are in DEGREES, as in the case files and SMLib's analytic
 * surfaces. Each surface's UV domain is SMLib's natural UV domain of the surface built
 * (smlib_surface_domain); for planes that is (u, v) itself, for the other analytic kinds it is
 * the parameterisation of SMLib's internal NURBS form, as the runner uses it.
 *
 * Every input is validated first; anything non-finite, out of range or inconsistent returns
 * SMLIB_ERR_ARG without touching SMLib. Coordinates, lengths and plane ranges are limited to
 * +-SMLIB_SURFACE_MAX_COORD. Axes must be non-zero and not parallel; they need not be unit or
 * orthogonal (SMLib normalises them as the runner's surfaces do). */

/* A surface with its UV domain, owned by the caller (smlib_surface_free). */
typedef struct smlib_surface smlib_surface;

#define SMLIB_SURFACE_MAX_COORD 1e8
#define SMLIB_SURFACE_MAX_DEGREE 25
/* Most control points along one direction, and in all (nu * nv). */
#define SMLIB_SURFACE_MAX_POLES 1000000

/* A plane through `origin` spanned by `xaxis` and `yaxis`: P(u, v) = origin + u*x + v*y with x, y
 * the axes made unit, for u in [u[0], u[1]], v in [v[0], v[1]] (each u[0] < u[1]). */
int32_t smlib_surface_plane(const double origin[3], const double xaxis[3], const double yaxis[3], const double u[2],
                            const double v[2], smlib_surface **out);

/* A cylinder of `radius` about the axis through `origin` along z = unit(xaxis x yaxis), angle
 * `angle_deg` = [a0, a1] from xaxis towards yaxis (a0 < a1, within [-360, 360], a1 - a0 <= 360),
 * for z in `height` = [h0, h1] (h0 < h1) measured from `origin`. */
int32_t smlib_surface_cylinder(const double origin[3], const double xaxis[3], const double yaxis[3], double radius,
                               const double angle_deg[2], const double height[2], smlib_surface **out);

/* A full (360 degree) cone with its apex at `apex`, axis z = unit(xaxis x yaxis), radius
 * h * tan(half_angle_deg) at distance h along z, for h in `height` = [h0, h1], 0 <= h0 < h1
 * (h0 = 0 includes the apex). half_angle_deg in [0.01, 89]. */
int32_t smlib_surface_cone(const double apex[3], const double xaxis[3], const double yaxis[3], double half_angle_deg,
                           const double height[2], smlib_surface **out);

/* A sphere of `radius` at `center`: longitude `lon_deg` (from xaxis towards yaxis; within
 * [-360, 360], span <= 360) and latitude `lat_deg` (within [-90, 90]), each increasing. */
int32_t smlib_surface_sphere(const double center[3], const double xaxis[3], const double yaxis[3], double radius,
                             const double lon_deg[2], const double lat_deg[2], smlib_surface **out);

/* A ring torus at `center` about z = unit(xaxis x yaxis): major radius R (the case files' "R"),
 * minor radius r ("r"), 0 < r < R; `major_deg` the angle round z from xaxis, `minor_deg` the
 * angle round the tube (each within [-360, 360], increasing, span <= 360). */
int32_t smlib_surface_torus(const double center[3], const double xaxis[3], const double yaxis[3], double major_radius,
                            double minor_radius, const double major_deg[2], const double minor_deg[2], smlib_surface **out);

/* A (rational) B-spline surface of degrees 1..SMLIB_SURFACE_MAX_DEGREE with nu x nv control
 * points (nu > degree_u, nv > degree_v). `points` holds nu * nv Cartesian points (x, y, z),
 * point (i, j) (i along u, j along v) at points[3 * (i * nv + j)], as in the case files;
 * `weights` likewise nu * nv weights in [1e-6, 1e6], or NULL for all 1. Knot vectors hold
 * nu + degree_u + 1 (nv + degree_v + 1) finite, non-decreasing values with a non-empty domain,
 * no end knot repeated more than degree + 1 times and no interior knot more than degree times.
 * Unclamped knot vectors are accepted. */
int32_t smlib_surface_nurbs(int32_t degree_u, int32_t degree_v, int32_t nu, int32_t nv, const double *knots_u,
                            int64_t n_knots_u, const double *knots_v, int64_t n_knots_v, const double *points,
                            const double *weights, smlib_surface **out);

void smlib_surface_free(smlib_surface *s);

/* The surface's UV domain as [u0, u1, v0, v1]: what smlib_ssi intersects within. */
int32_t smlib_surface_domain(const smlib_surface *s, double out[4]);

/* The surface's point at (u, v) in its UV domain's parameterisation. */
int32_t smlib_surface_evaluate(const smlib_surface *s, double u, double v, double out[3]);

/* smlib_ssi methods. */
#define SMLIB_SSI_GLOBAL 0   /* SmSurface::GlobalSurfaceIntersect: analytic special cases, else 1 */
#define SMLIB_SSI_ADVANCED 1 /* SmAdvSurfaceIntersector (the general NURBS path) */
#define SMLIB_SSI_BASIC 2    /* SmSurfaceIntersector (the basic intersector) */

/* Curve types: SMLib's SmTsectCurveType. */
#define SMLIB_SSI_NONE (-1)          /* the method reported no type for this curve */
#define SMLIB_SSI_UNKNOWN 0
#define SMLIB_SSI_TOUCHING 1         /* a degenerate curve: one point where the surfaces touch */
#define SMLIB_SSI_CROSSING 2         /* a transversal curve */
#define SMLIB_SSI_TANGENT 3          /* a curve along which the surfaces are tangent */
#define SMLIB_SSI_COINCIDENT 4       /* where the surfaces butt up against each other */
#define SMLIB_SSI_NEAR_TANGENT 5     /* a curve with a small angle of intersection */
#define SMLIB_SSI_REGION_BOUNDARY 6  /* bounds a region where the surfaces coincide */

#define SMLIB_SSI_MAX_SAMPLES 100000

/* The intersection of `a` and `b`, each within its UV domain, owned by the caller
 * (smlib_ssi_result_free). */
typedef struct smlib_ssi_result smlib_ssi_result;

/* Intersects `a` and `b` (both surfaces' edges used to seed curves) with `method` at
 * approximation tolerance `approx_tol` (<= 0: SMLib's default, 1/1000 of the surfaces' size)
 * and samples each curve at `samples_per_curve` (2..SMLIB_SSI_MAX_SAMPLES) parameters, as the
 * corpus runner does (it uses 4001): the 3D curve at t_i = t0 + (t1 - t0) * i / (n - 1) over
 * its natural interval [t0, t1], each UV curve likewise over its own natural interval, and `a`
 * and `b` evaluated at the UV samples. A result over 10^7 samples in all returns SMLIB_ERR_ARG.
 * The status is SMLib's when the intersector fails. On failure the dump holds both surfaces as
 * an SSI corpus case file (`<n>-smlib_ssi-case.txt`) that the corpus runner replays. */
int32_t smlib_ssi(const smlib_surface *a, const smlib_surface *b, int32_t method, double approx_tol,
                  int32_t samples_per_curve, smlib_ssi_result **out);

void smlib_ssi_result_free(smlib_ssi_result *r);

int32_t smlib_ssi_curve_count(const smlib_ssi_result *r, int32_t *out);

/* smlib_ssi_curve_info flags: which UV curves the method returned. GlobalSurfaceIntersect's
 * analytic special cases (plane, cylinder, cone, sphere, torus pairs) return none. */
#define SMLIB_SSI_HAS_UV_A 1
#define SMLIB_SSI_HAS_UV_B 2

/* Curve `i`: its SMLIB_SSI_* type, SMLIB_SSI_HAS_UV_* flags, sample count (samples_per_curve),
 * SMLib's deviation estimate and the 3D curve's natural parameter interval. Any out-pointer may
 * be NULL. */
int32_t smlib_ssi_curve_info(const smlib_ssi_result *r, int32_t i, int32_t *type, int32_t *flags, int32_t *samples,
                             double *deviation, double interval[2]);

/* Copies curve `i`'s samples: `xyz` (3 per sample), `uv_a` and `uv_b` (2 per sample, the UV
 * curves on a and b), `a_at_uv` and `b_at_uv` (3 per sample, a and b evaluated at those UVs).
 * Each buffer may be NULL to skip it; `cap` is the buffers' capacity in samples, and less than
 * the sample count returns SMLIB_ERR_BUFFER. A value SMLib could not evaluate, and every UV
 * value of a UV curve the method did not return (with the surface at it), is NaN; the corpus
 * runner drops such samples instead (and writes empty arrays for a missing UV curve). */
int32_t smlib_ssi_curve_samples(const smlib_ssi_result *r, int32_t i, int64_t cap, double *xyz, double *uv_a, double *uv_b,
                                double *a_at_uv, double *b_at_uv);

/* Isolated points: the TOUCHING curves again, each as its first sample (UVs NaN where the
 * method returned no UV curve). */
int32_t smlib_ssi_point_count(const smlib_ssi_result *r, int32_t *out);
int32_t smlib_ssi_point(const smlib_ssi_result *r, int32_t i, double xyz[3], double uv_a[2], double uv_b[2]);

/* ---- Sweeps (src/smlib_sweeps.cpp) -------------------------------------------------------- */

/* Curve lists. Profiles, paths and sections cross as a list of 3D curves in two arrays:
 * `kinds` holds, for each curve in turn, its kind, and for a NURBS curve then its degree and its
 * number of control points n; `reals` holds each curve's numbers in turn:
 *   SMLIB_CURVE_LINE   start[3], end[3]                                                (6)
 *   SMLIB_CURVE_ARC    center[3], axis[3], ref[3], radius, start, end                   (12)
 *                      the arc counter-clockwise about the unit `axis` from angle `start`
 *                      to `end` (radians, 0 < end - start <= 2 pi), angle 0 along the unit
 *                      `ref` (perpendicular to `axis`)
 *   SMLIB_CURVE_NURBS  n control points (x, y, z, weight > 0), then n + degree + 1 knots,
 *                      non-decreasing, clamped (end multiplicity degree + 1, inner ones at
 *                      most degree)                                           (5n + degree + 1)
 * Both arrays must be used up exactly. At most SMLIB_CURVE_MAX_CURVES curves, coordinates within
 * 1e9 mm. A list that breaks any of these returns SMLIB_ERR_ARG. A profile's curves run end to
 * end in one or more closed loops (one closed planar loop and its holes for an extrusion). */
#define SMLIB_CURVE_LINE 1
#define SMLIB_CURVE_ARC 2
#define SMLIB_CURVE_NURBS 3
#define SMLIB_CURVE_MAX_CURVES 4096
#define SMLIB_CURVE_MAX_DEGREE 25
#define SMLIB_CURVE_MAX_POINTS 100000

/* Each sweep returns a new body. With `cap` 1 the profile's loops close it at both ends and the
 * result must be a closed manifold solid (SMLIB_ERR_NOT_SOLID otherwise); with `cap` 0 it is the
 * swept sheet alone. */

/* The planar profile extruded `distance` (> 0) along the unit `direction`, square to the
 * profile's plane. `draft` (radians, |draft| < pi/2) leans the walls: 0 straight, positive
 * inward (the far end smaller), negative outward; a drafted profile is one closed loop, its
 * curves in order, and its direction must be the profile's normal (either way round). */
int32_t smlib_extrude(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, const double direction[3],
                      double distance, double draft, int32_t cap, smlib_brep **out);

/* The planar profile turned `angle` radians (0 < angle <= 2 pi) right-handed about the axis
 * through `origin` along `direction`, which lies in the profile's plane; the profile may touch
 * the axis but not cross it. */
int32_t smlib_revolve(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, const double origin[3],
                      const double direction[3], double angle, int32_t cap, smlib_brep **out);

/* The profile turns with the path (keeping its attitude to the path's tangent). */
#define SMLIB_SWEEP_FOLLOW 0
/* The profile only moves along the path, keeping its orientation. */
#define SMLIB_SWEEP_TRANSLATE 1

/* The profile swept along a path of end-to-end curves (joined into one B-spline), carried with
 * the path's start: the profile keeps its place relative to it. The side walls are B-spline
 * approximations of the swept surface, one per stretch of profile between tangent breaks; along a
 * straight path (lines end to end in one direction) the sweep is an extrusion, and turning with a
 * path of arcs about one axis a revolution, with exact walls.
 * Round a closed path (ending where it starts) the walls close on each other: no caps. */
int32_t smlib_sweep(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, const int64_t *path_kinds,
                    int64_t n_path_kinds, const double *path_reals, int64_t n_path_reals, int32_t frame, int32_t cap,
                    smlib_brep **out);

/* smlib_sweep, also saying where it placed the profile: `placed` (36 numbers, or NULL for none)
 * gets three rigid motions, 12 numbers each as smlib_sweep_stations takes its maps (a 3 x 3
 * matrix by columns, then a translation), carrying the profile as given onto the sweep's first
 * section, a section about half way along (the middle one of those the sweep places) and its last
 * section. The curve sweep turns the profile on its own frame, computed from point to point along
 * the path, so only the sweep can say where it leaves the profile; a caller that carries on from
 * there (a path with corners swept piece by piece) takes the last motion. Along a straight path
 * the motions are the translations to the start, middle and end; turning with arcs about one
 * axis, the turns by none, half and all of the angle. Meaningful only on success. */
int32_t smlib_sweep_placed(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, const int64_t *path_kinds,
                           int64_t n_path_kinds, const double *path_reals, int64_t n_path_reals, int32_t frame, int32_t cap,
                           smlib_brep **out, double *placed);

/* A loft (skin) through `n_sections` (>= 2) sections, each a closed loop of curves: the curve
 * list holds the sections' curves one section after the other, `section_sizes[i]` curves for
 * section i. A section with holes holds its loops one after the other, its outline first; each loop
 * is lofted to the loop in the same place in the next section. `degree` (1 ruled, 2, 3) is the
 * walls' degree across the sections. With `cap` 1 the first and last sections close the ends. */
int32_t smlib_loft(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, const int64_t *section_sizes,
                   int64_t n_sections, int32_t degree, int32_t cap, smlib_brep **out);

/* A loft from a planar section to a point off its plane (a pyramid): the section, one closed loop
 * of lines in order, closes the base, and each line spans a triangle to `apex`. A section with a
 * curve, or an apex in its plane, is SMLIB_ERR_ARG. */
int32_t smlib_loft_to_point(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, const double apex[3], smlib_brep **out);

/* A sheet of planar faces bounded by the curve list's closed loops, all in one plane and in any
 * order (SmBrep::CreatePlanarFacesWith3DCurves): each outer loop one face, the loops directly
 * inside it its holes, a loop inside a hole the outer loop of another face. The loops must not
 * cross; loops meeting at an end point (SMLib's spine points) are SM_ERR. */
int32_t smlib_planar_faces(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, smlib_brep **out);

/* The profile carried by a screw motion about the axis through `origin` along `direction`:
 * `turns` (0.001 to SMLIB_SCREW_MAX_TURNS) turns right-handed about the axis, rising `pitch`
 * (signed, non-zero) per turn, capped at both ends: a thread's groove or a coil. The walls are
 * skinned (cubic) through the profile placed by the screw motion at stations close enough that
 * they stay within about `tolerance` (>= 1e-6 mm) of the true screw surfaces; at most
 * SMLIB_SCREW_MAX_STATIONS stations, else SMLIB_ERR_ARG. */
#define SMLIB_SCREW_MAX_TURNS 200.0
#define SMLIB_SCREW_MAX_STATIONS 12000.0
int32_t smlib_screw_sweep(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, const double origin[3],
                          const double direction[3], double pitch, double turns, double tolerance, smlib_brep **out);

/* The ends of smlib_sweep_stations: left open (the swept sheet), capped by the first and last
 * sections (a closed manifold solid), or joined, the last section lying on the first (a sweep
 * round a closed path: the walls' ends are glued to their starts into a closed manifold solid). */
#define SMLIB_ENDS_OPEN 0
#define SMLIB_ENDS_CAPPED 1
#define SMLIB_ENDS_JOINED 2
#define SMLIB_SWEEP_MAX_STATIONS 4096
#define SMLIB_SWEEP_MAX_STATION_CURVES 400000.0
/* The planar profile carried through `n_stations` (2 to SMLIB_SWEEP_MAX_STATIONS) stations, each
 * an affine map of 12 numbers in `maps`: a 3 x 3 matrix by columns, then a translation (a point x
 * goes to M x + t). Each section is the profile moved by its map, exactly (its curves' control
 * points mapped); the walls are skinned through the sections in order, cubic across them (a lower
 * degree with fewer than four), one wall per profile curve. A sweep whose profile twists, scales
 * or follows guide rails along its path is the profile placed at stations along it. Every map
 * must be finite, its matrix far from singular and its columns shorter than 1e6, its translation
 * within 1e9 mm; profile curves times stations at most SMLIB_SWEEP_MAX_STATION_CURVES; `ends` one
 * of SMLIB_ENDS_*. Otherwise SMLIB_ERR_ARG. */
int32_t smlib_sweep_stations(const int64_t *kinds, int64_t n_kinds, const double *reals, int64_t n_reals, const double *maps,
                             int64_t n_stations, int32_t ends, smlib_brep **out);

/* ---- Imprint (src/smlib_imprint.cpp) ------------------------------------------------------ */

/* The tool cut nothing: no face of the body meets it. */
#define SMLIB_ERR_NO_CUT (-1104)

/* A new body: `body` with each of its faces that `tool` (a sheet or a solid) meets cut along the
 * intersection, the pieces kept as faces of one body. Both inputs are left unchanged. A tool that
 * meets no face is SMLIB_ERR_NO_CUT. Trim and Split Face build on it: they choose which pieces to
 * keep. */
int32_t smlib_imprint(const smlib_brep *body, const smlib_brep *tool, smlib_brep **out);

/* ---- Split (src/smlib_split.cpp) --------------------------------------------------------- */

#define SMLIB_SPLIT_MAX_FACES 100000
#define SMLIB_SPLIT_MAX_TOOLS 4096

/* A new body: `body` with its faces `faces` (`n_faces` distinct indices in its face order) cut
 * where each tool surface (over its whole domain) meets them, the tools taken in turn, each one
 * cutting the pieces the earlier ones left. The curves are SMLib's intersections of each face's
 * surface with the tool (a plane's sections exactly), trimmed to the face; the edges they end on
 * are cut in every face that uses them. `made[i]` (n_faces of them) is how many faces face
 * `faces[i]` became besides itself: 0 where no tool split it. No face split at all is
 * SMLIB_ERR_NO_CUT. The input is left unchanged. */
int32_t smlib_split_faces(const smlib_brep *body, const int64_t *faces, int64_t n_faces, const smlib_surface *const *tools,
                          int64_t n_tools, smlib_brep **out, int64_t *made);

/* The manifold solid `body` cut in two by the tool surface (over its whole domain, which must
 * reach past the body): `back` the part on the side the surface's normal points away from,
 * `front` the part it points into, each closed by the tool's piece inside the body (SMLib's
 * slice). A part may hold several solids. A tool that misses the body leaves each part the whole
 * body, with no face of the tool (tag the body, smlib_history_tag, to tell); one that leaves a side
 * empty is SMLIB_ERR_NO_CUT; a body that is no manifold solid SMLIB_ERR_NOT_SOLID. */
int32_t smlib_split_body(const smlib_brep *body, const smlib_surface *tool, smlib_brep **back, smlib_brep **front);

/* ---- Extend (src/smlib_extend.cpp) -------------------------------------------------------- */

/* The input is a shape the call does not take (here: a sheet that is not one face bounded by
 * its surface's parameter lines, or an edge that is not a whole side of it). */
#define SMLIB_ERR_UNSUPPORTED (-1105)
/* The answer is not exact in the input's geometry (here: an analytic surface whose speed across
 * the side varies along it, or a box that would run past a pole or wrap round). */
#define SMLIB_ERR_NOT_EXACT (-1106)

/* How a surface continues past its side: its own geometry (a B-spline's end spans continued), a
 * G1 run-off, or a G2 run-off. */
#define SMLIB_EXTEND_NATURAL 0
#define SMLIB_EXTEND_TANGENT 1
#define SMLIB_EXTEND_CURVATURE 2

/* A new body: the one-face sheet `sheet`, whose face is the whole parameter box of its surface
 * (for a B-spline, all of the surface), continued `distance` (> 0) past the side its edge
 * nearest `point` (within `tol`) runs along, as `mode` says. An analytic surface keeps its
 * geometry and its box grows by the parameter step that crosses `distance` (the same all along
 * the side, else SMLIB_ERR_NOT_EXACT); only a plane takes the G1 and G2 modes (as itself). A
 * B-spline surface is SMLib's extended surface. The face keeps its orientation. */
int32_t smlib_extend_sheet(const smlib_brep *sheet, const double point[3], double tol, double distance, int32_t mode, smlib_brep **out);

/* ---- Untrim (src/smlib_untrim.cpp) -------------------------------------------------------- */

/* A new body: the face of `body` nearest `point` (within `tol`) restored towards its surface.
 * The holes taken away are the inner loops nearest each of the `n_holes` points `holes` (each
 * within `tol`, else SMLIB_ERR_PICK_NONE), or every inner loop when `n_holes` < 0. With `outer`,
 * the outer boundary goes back to the surface's parameter box (a free-form surface's own; an
 * analytic surface's the box of the face's outer loop in its parameters), bounded by the box's
 * iso-parameter edges (a box all the way round a periodic direction closes on the surface's
 * seam), the holes kept rebuilt in it from their own edges. The result is the face alone, or
 * with `in_place` the whole body with the face replaced (a hole another face uses is then
 * SMLIB_ERR_UNSUPPORTED). The face keeps its orientation and the attributes of the face and
 * edges it is made from (history tags). Nothing to do (no holes and not `outer`) is
 * SMLIB_ERR_ARG. */
int32_t smlib_untrim_face(const smlib_brep *body, const double point[3], double tol, const double *holes, int64_t n_holes, int32_t outer,
                          int32_t in_place, smlib_brep **out);

#ifdef __cplusplus
}
#endif

#endif /* SMLIB_C_H */
