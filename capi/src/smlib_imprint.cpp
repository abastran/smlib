// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_imprint: where a tool body meets a body, its faces cut along the intersection
// (SmApiBooleanWithOptions with SM_BO_IMPRINT, the tool kept). Trim and Split Face build on it.
// smlib_face_points: a point inside each face (SmFace::CalculateAnInternalPoint, which classifies
// parameter lines against the face's trims), for telling the pieces apart.

#include "smlib_shim.h"

#include <SmApiBrep.h>
#include <SmApiQueries.h>
#include <SmFace.h>

#include <memory>

using namespace smlib_shim;

namespace
{

using Owned = std::unique_ptr<SmBrep>;

Owned copy_of(const smlib_brep *b, int32_t &s)
{
    SmBrep *r = nullptr;
    s = status(SmApiBrepCopy(b->p, r));
    if (s != SMLIB_OK)
    {
        delete r;
        return nullptr;
    }
    if (r == nullptr)
    {
        s = SMLIB_ERR_EMPTY;
    }
    return Owned(r);
}

} // namespace

extern "C" int32_t smlib_imprint(const smlib_brep *body, const smlib_brep *tool, smlib_brep **out)
{
    if (!valid(body) || !valid(tool) || out == nullptr || body == tool)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_imprint", handle(body) + ", " + handle(tool));
    c.input("body", body);
    c.input("tool", tool);
    return c.run(
        [&]
        {
            int32_t s = SMLIB_OK;
            Owned a = copy_of(body, s);
            if (s != SMLIB_OK)
            {
                return s;
            }
            Owned t = copy_of(tool, s);
            if (s != SMLIB_OK)
            {
                return s;
            }
            ULONG before = a->GetNumFaces();
            SmBrep *r = nullptr;
            // Imprint the intersection on A and stop; keep the tool (freed here).
            s = status(SmApiBooleanWithOptions(a.get(), t.get(), SM_BO_IMPRINT, FALSE, TRUE, FALSE, FALSE, TRUE, r));
            if (s != SMLIB_OK)
            {
                return s;
            }
            if (r != nullptr && r != a.get())
            {
                // SMLib handed back another body and owns A now.
                a.release();
                a.reset(r);
            }
            if (a == nullptr || a->GetNumFaces() == 0)
            {
                return SMLIB_ERR_EMPTY;
            }
            if (a->GetNumFaces() == before)
            {
                return SMLIB_ERR_NO_CUT;
            }
            return take(a.release(), out);
        }
    );
}

extern "C" int32_t smlib_face_points(const smlib_brep *body, double *out, int64_t n_faces)
{
    if (!valid(body) || out == nullptr || n_faces < 0)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_face_points", handle(body) + ", n_faces=" + std::to_string(n_faces));
    c.input("body", body);
    return c.run(
        [&]
        {
            SmTArray<SmFace *> faces;
            body->p->GetFaces(faces);
            if (static_cast<int64_t>(faces.GetSize()) != n_faces)
            {
                return SMLIB_ERR_ARG;
            }
            for (ULONG i = 0; i < faces.GetSize(); ++i)
            {
                SmPoint3d p;
                int32_t s = status(faces[i]->CalculateAnInternalPoint(p));
                if (s != SMLIB_OK)
                {
                    return s;
                }
                out[3 * i] = p.x;
                out[3 * i + 1] = p.y;
                out[3 * i + 2] = p.z;
            }
            return SMLIB_OK;
        }
    );
}
