// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_imprint: where a tool body meets a body, its faces cut along the intersection
// (SmApiBooleanWithOptions with SM_BO_IMPRINT, the tool kept). Trim and Split Face build on it.

#include "smlib_shim.h"

#include <SmApiBrep.h>
#include <SmApiQueries.h>

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
