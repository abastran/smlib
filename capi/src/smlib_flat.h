// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// Flat, USD-free arrays for one SmBrep, and the conversions both ways.
//
// Adapted from source/BREP_SM_USD/src/SmuConvert.cpp (NVIDIA, Apache-2.0) without USD types;
// SolidCraft change. The arrays mirror one member of UsdBrepData::UsdBrepArrayData
// (source/BREP_USD_DATA/inc/UsdBrepArrayData.h); include/smlib.h documents each one.

#ifndef SMLIB_FLAT_H
#define SMLIB_FLAT_H

#include "smlib.h"

#include <SmBrep.h>
#include <SmTypes.h>

#include <cstdint>
#include <string>
#include <vector>

// Every array: X(id, member). I lists the int64 arrays, F the double ones.
#define SMLIB_FLAT_FIELDS(I, F)                                                                    \
    I(SMLIB_FLAT_REGION_SHELL_COUNT, region_shell_count)                                           \
    I(SMLIB_FLAT_REGION_TYPE, region_type)                                                         \
    I(SMLIB_FLAT_SHELL_FACEUSE_COUNT, shell_faceuse_count)                                         \
    I(SMLIB_FLAT_SHELL_WIRE_EDGE_COUNT, shell_wire_edge_count)                                     \
    I(SMLIB_FLAT_SHELL_POINT_TYPE, shell_point_type)                                               \
    I(SMLIB_FLAT_FACEUSE_FACE, faceuse_face)                                                       \
    I(SMLIB_FLAT_FACEUSE_ORIENTATION, faceuse_orientation)                                         \
    I(SMLIB_FLAT_FACE_LOOP_COUNT, face_loop_count)                                                 \
    I(SMLIB_FLAT_FACE_SURFACE_TYPE, face_surface_type)                                             \
    I(SMLIB_FLAT_FACE_TRIM_TYPE, face_trim_type)                                                   \
    I(SMLIB_FLAT_LOOP_EDGEUSE_COUNT, loop_edgeuse_count)                                           \
    I(SMLIB_FLAT_LOOP_VERTEX, loop_vertex)                                                         \
    I(SMLIB_FLAT_EDGEUSE_EDGE, edgeuse_edge)                                                       \
    I(SMLIB_FLAT_EDGEUSE_ORIENTATION, edgeuse_orientation)                                         \
    I(SMLIB_FLAT_EDGEUSE_NEXT_RADIAL, edgeuse_next_radial)                                         \
    I(SMLIB_FLAT_EDGEUSE_RADIAL_ENTRY, edgeuse_radial_entry)                                       \
    I(SMLIB_FLAT_EDGE_CURVE_TYPE, edge_curve_type)                                                 \
    I(SMLIB_FLAT_EDGE_VERTICES, edge_vertices)                                                     \
    I(SMLIB_FLAT_WIRE_EDGE_CURVE_TYPE, wire_edge_curve_type)                                       \
    I(SMLIB_FLAT_WIRE_EDGE_VERTICES, wire_edge_vertices)                                           \
    I(SMLIB_FLAT_VERTEX_POINT_TYPE, vertex_point_type)                                             \
    I(SMLIB_FLAT_EDGE_NURBS_VERTEX_COUNT, edge_nurbs_vertex_count)                                 \
    I(SMLIB_FLAT_EDGE_NURBS_ORDER, edge_nurbs_order)                                               \
    I(SMLIB_FLAT_WIRE_EDGE_NURBS_VERTEX_COUNT, wire_edge_nurbs_vertex_count)                       \
    I(SMLIB_FLAT_WIRE_EDGE_NURBS_ORDER, wire_edge_nurbs_order)                                     \
    I(SMLIB_FLAT_EDGEUSE_NURBS_VERTEX_COUNT, edgeuse_nurbs_vertex_count)                           \
    I(SMLIB_FLAT_EDGEUSE_NURBS_ORDER, edgeuse_nurbs_order)                                         \
    I(SMLIB_FLAT_FACE_NURBS_U_VERTEX_COUNT, face_nurbs_u_vertex_count)                             \
    I(SMLIB_FLAT_FACE_NURBS_V_VERTEX_COUNT, face_nurbs_v_vertex_count)                             \
    I(SMLIB_FLAT_FACE_NURBS_U_ORDER, face_nurbs_u_order)                                           \
    I(SMLIB_FLAT_FACE_NURBS_V_ORDER, face_nurbs_v_order)                                           \
    I(SMLIB_FLAT_SOURCE, source)                                                                   \
    F(SMLIB_FLAT_TOLERANCE, tolerance)                                                             \
    F(SMLIB_FLAT_EXTENT, extent)                                                                   \
    F(SMLIB_FLAT_FACE_RANGE, face_range)                                                           \
    F(SMLIB_FLAT_EDGE_RANGE, edge_range)                                                           \
    F(SMLIB_FLAT_WIRE_EDGE_RANGE, wire_edge_range)                                                 \
    F(SMLIB_FLAT_SHELL_POINT_POSITION, shell_point_position)                                       \
    F(SMLIB_FLAT_VERTEX_POSITION, vertex_position)                                                 \
    F(SMLIB_FLAT_EDGE_NURBS_CVS, edge_nurbs_cvs)                                                   \
    F(SMLIB_FLAT_EDGE_NURBS_KNOTS, edge_nurbs_knots)                                               \
    F(SMLIB_FLAT_EDGE_NURBS_WEIGHTS, edge_nurbs_weights)                                           \
    F(SMLIB_FLAT_EDGE_CIRCLE_CENTER, edge_circle_center)                                           \
    F(SMLIB_FLAT_EDGE_CIRCLE_AXIS, edge_circle_axis)                                               \
    F(SMLIB_FLAT_EDGE_CIRCLE_REF_DIRECTION, edge_circle_ref_direction)                             \
    F(SMLIB_FLAT_EDGE_CIRCLE_RADIUS, edge_circle_radius)                                           \
    F(SMLIB_FLAT_EDGE_LINE_ORIGIN, edge_line_origin)                                               \
    F(SMLIB_FLAT_EDGE_LINE_DIRECTION, edge_line_direction)                                         \
    F(SMLIB_FLAT_EDGE_ELLIPSE_CENTER, edge_ellipse_center)                                         \
    F(SMLIB_FLAT_EDGE_ELLIPSE_AXIS, edge_ellipse_axis)                                             \
    F(SMLIB_FLAT_EDGE_ELLIPSE_REF_DIRECTION, edge_ellipse_ref_direction)                           \
    F(SMLIB_FLAT_EDGE_ELLIPSE_X_RADIUS, edge_ellipse_x_radius)                                     \
    F(SMLIB_FLAT_EDGE_ELLIPSE_Y_RADIUS, edge_ellipse_y_radius)                                     \
    F(SMLIB_FLAT_WIRE_EDGE_NURBS_CVS, wire_edge_nurbs_cvs)                                         \
    F(SMLIB_FLAT_WIRE_EDGE_NURBS_KNOTS, wire_edge_nurbs_knots)                                     \
    F(SMLIB_FLAT_WIRE_EDGE_NURBS_WEIGHTS, wire_edge_nurbs_weights)                                 \
    F(SMLIB_FLAT_WIRE_EDGE_CIRCLE_CENTER, wire_edge_circle_center)                                 \
    F(SMLIB_FLAT_WIRE_EDGE_CIRCLE_AXIS, wire_edge_circle_axis)                                     \
    F(SMLIB_FLAT_WIRE_EDGE_CIRCLE_REF_DIRECTION, wire_edge_circle_ref_direction)                   \
    F(SMLIB_FLAT_WIRE_EDGE_CIRCLE_RADIUS, wire_edge_circle_radius)                                 \
    F(SMLIB_FLAT_WIRE_EDGE_LINE_ORIGIN, wire_edge_line_origin)                                     \
    F(SMLIB_FLAT_WIRE_EDGE_LINE_DIRECTION, wire_edge_line_direction)                               \
    F(SMLIB_FLAT_WIRE_EDGE_ELLIPSE_CENTER, wire_edge_ellipse_center)                               \
    F(SMLIB_FLAT_WIRE_EDGE_ELLIPSE_AXIS, wire_edge_ellipse_axis)                                   \
    F(SMLIB_FLAT_WIRE_EDGE_ELLIPSE_REF_DIRECTION, wire_edge_ellipse_ref_direction)                 \
    F(SMLIB_FLAT_WIRE_EDGE_ELLIPSE_X_RADIUS, wire_edge_ellipse_x_radius)                           \
    F(SMLIB_FLAT_WIRE_EDGE_ELLIPSE_Y_RADIUS, wire_edge_ellipse_y_radius)                           \
    F(SMLIB_FLAT_EDGEUSE_NURBS_CVS, edgeuse_nurbs_cvs)                                             \
    F(SMLIB_FLAT_EDGEUSE_NURBS_KNOTS, edgeuse_nurbs_knots)                                         \
    F(SMLIB_FLAT_EDGEUSE_NURBS_WEIGHTS, edgeuse_nurbs_weights)                                     \
    F(SMLIB_FLAT_FACE_NURBS_CVS, face_nurbs_cvs)                                                   \
    F(SMLIB_FLAT_FACE_NURBS_U_KNOTS, face_nurbs_u_knots)                                           \
    F(SMLIB_FLAT_FACE_NURBS_V_KNOTS, face_nurbs_v_knots)                                           \
    F(SMLIB_FLAT_FACE_NURBS_WEIGHTS, face_nurbs_weights)                                           \
    F(SMLIB_FLAT_FACE_SPHERE_CENTER, face_sphere_center)                                           \
    F(SMLIB_FLAT_FACE_SPHERE_AXIS, face_sphere_axis)                                               \
    F(SMLIB_FLAT_FACE_SPHERE_REF_DIRECTION, face_sphere_ref_direction)                             \
    F(SMLIB_FLAT_FACE_SPHERE_RADIUS, face_sphere_radius)                                           \
    F(SMLIB_FLAT_FACE_PLANE_ORIGIN, face_plane_origin)                                             \
    F(SMLIB_FLAT_FACE_PLANE_AXIS, face_plane_axis)                                                 \
    F(SMLIB_FLAT_FACE_PLANE_REF_DIRECTION, face_plane_ref_direction)                               \
    F(SMLIB_FLAT_FACE_CYLINDER_ORIGIN, face_cylinder_origin)                                       \
    F(SMLIB_FLAT_FACE_CYLINDER_AXIS, face_cylinder_axis)                                           \
    F(SMLIB_FLAT_FACE_CYLINDER_REF_DIRECTION, face_cylinder_ref_direction)                         \
    F(SMLIB_FLAT_FACE_CYLINDER_RADIUS, face_cylinder_radius)                                       \
    F(SMLIB_FLAT_FACE_CONE_ORIGIN, face_cone_origin)                                               \
    F(SMLIB_FLAT_FACE_CONE_AXIS, face_cone_axis)                                                   \
    F(SMLIB_FLAT_FACE_CONE_REF_DIRECTION, face_cone_ref_direction)                                 \
    F(SMLIB_FLAT_FACE_CONE_RADIUS, face_cone_radius)                                               \
    F(SMLIB_FLAT_FACE_CONE_SEMI_ANGLE, face_cone_semi_angle)                                       \
    F(SMLIB_FLAT_FACE_TORUS_ORIGIN, face_torus_origin)                                             \
    F(SMLIB_FLAT_FACE_TORUS_AXIS, face_torus_axis)                                                 \
    F(SMLIB_FLAT_FACE_TORUS_REF_DIRECTION, face_torus_ref_direction)                               \
    F(SMLIB_FLAT_FACE_TORUS_MAJOR_RADIUS, face_torus_major_radius)                                 \
    F(SMLIB_FLAT_FACE_TORUS_MINOR_RADIUS, face_torus_minor_radius)                                 \
    F(SMLIB_FLAT_EDGE_TOLERANCE, edge_tolerance)                                                   \
    F(SMLIB_FLAT_VERTEX_TOLERANCE, vertex_tolerance)

namespace smlib_flat_impl
{

// One B-rep's arrays, the per-brep slice of UsdBrepArrayData.
struct FlatBrep
{
#define SMLIB_FLAT_DECL_I(id, m) std::vector<int64_t> m;
#define SMLIB_FLAT_DECL_F(id, m) std::vector<double> m;
    SMLIB_FLAT_FIELDS(SMLIB_FLAT_DECL_I, SMLIB_FLAT_DECL_F)
#undef SMLIB_FLAT_DECL_I
#undef SMLIB_FLAT_DECL_F
};

// The array named by a public id: exactly one of `i` and `f` is set, or neither for an id
// that names no array.
struct Field
{
    int32_t id;
    const char *name; // the id without its SMLIB_FLAT_ prefix
    std::vector<int64_t> FlatBrep::*i;
    std::vector<double> FlatBrep::*f;
};

// Every array, in id order.
const std::vector<Field> &Fields();

// The field for `id`, or nullptr.
const Field *FindField(int32_t id);

// SmBrep -> arrays (port of SMU_BrepConvert::BrepAppend_SMLibToUsdBrep for one B-rep, with
// UV curves exported and face ranges left as SMLib has them, SmuConvert's defaults). The B-rep
// is only read.
SmStatus FlatFromBrep(const SmBrep &brep, FlatBrep &flat);

// Arrays -> new SmBrep in `context` (port of SMU_BrepConvert::BrepMove_OneUsdBrepToSMLibData
// and BrepMove_OneUsdBrepToSMLib). `flat` is validated first; malformed input returns
// SM_ERR_INVALID_INPUT. On failure `out` is NULL.
SmStatus BrepFromFlat(const FlatBrep &flat, const SmContext &context, SmBrep *&out);

// The text format of smlib_flat_write / smlib_flat_read and the import failure dumps.
bool WriteText(const FlatBrep &flat, const std::string &path);
bool ReadText(const std::string &path, FlatBrep &flat);

} // namespace smlib_flat_impl

#endif // SMLIB_FLAT_H
