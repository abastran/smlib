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
#define SMLIB_ABI_VERSION 2

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

/* An axis-aligned box from its lower corner `origin` with edge lengths `size`, in mm. */
int32_t smlib_make_box(const double origin[3], const double size[3], smlib_brep **out);

/* The enclosed volume at `relative_accuracy` (clamped by SMLib to [1e-8, 1e-1]). */
int32_t smlib_brep_volume(const smlib_brep *brep, double relative_accuracy, double *out);

/* Face, edge and vertex counts, in that order. */
int32_t smlib_brep_counts(const smlib_brep *brep, int64_t out[3]);

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
 * SMLib rebuilds them from the 3D edges on import; edge curves are always exported as NURBS
 * (lines, circles and ellipses as their exact NURBS form; other curve kinds approximated to the
 * B-rep's tolerance); surfaces that are neither analytic nor NURBS are approximated likewise.
 * Attributes, materials and colours are not carried.
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

/* FACE_SURFACE_TYPE values. */
#define SMLIB_FLAT_SURFACE_NURBS 0
#define SMLIB_FLAT_SURFACE_PLANE 1
#define SMLIB_FLAT_SURFACE_CYLINDER 2
#define SMLIB_FLAT_SURFACE_CONE 3
#define SMLIB_FLAT_SURFACE_SPHERE 4
#define SMLIB_FLAT_SURFACE_TORUS 5

/* EDGE_CURVE_TYPE and WIRE_EDGE_CURVE_TYPE values. Export writes only NURBS. */
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

#ifdef __cplusplus
}
#endif

#endif /* SMLIB_C_H */
