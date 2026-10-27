// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// libFuzzer target for smlib_brep_check: arbitrary levels, buffer sizes and NULLs, on bodies
// whose vertices the bytes may push anywhere (NaN and infinities included), so the check runs on
// broken geometry too. It may report anything; it must never crash, hang, leak or throw.
//
//   cmake -S capi -B build-fuzz -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DSMLIB_CAPI_FUZZ=ON
//   cmake --build build-fuzz --target smlib_fuzz_check && build-fuzz/smlib_fuzz_check -max_total_time=60

#include "../src/smlib_shim.h"

#include <SmApiPrimitives.h>
#include <SmVertex.h>

#include <cstddef>
#include <cstdint>
#include <cstring>

namespace
{

struct Input
{
    const uint8_t *data;
    size_t size;

    uint8_t byte()
    {
        if (size == 0)
        {
            return 0;
        }
        uint8_t b = *data;
        ++data;
        --size;
        return b;
    }

    double real()
    {
        double x = 0.0;
        size_t n = size < sizeof x ? size : sizeof x;
        std::memcpy(&x, data, n);
        data += n;
        size -= n;
        return x;
    }

    int64_t integer()
    {
        int64_t x = 0;
        size_t n = size < sizeof x ? size : sizeof x;
        std::memcpy(&x, data, n);
        data += n;
        size -= n;
        return x;
    }
};

SmBrep *primitive(uint8_t kind)
{
    const SmVector3d o(1.0, 2.0, 3.0);
    SmBrep *p = nullptr;
    switch (kind % 5)
    {
    case 0:
        SmApiCreateBox(o, 10.0, 20.0, 30.0, p);
        break;
    case 1:
        SmApiCreateCylinder(o, 15.0, 50.0, p);
        break;
    case 2:
        SmApiCreateSphere(o, 12.0, p);
        break;
    case 3:
        SmApiCreateCone(o, 10.0, 4.0, 20.0, p);
        break;
    default:
        SmApiCreateTorus(o, 22.0, 4.0, p);
        break;
    }
    return p;
}

} // namespace

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    Input in{data, size};
    if (smlib_init() != SMLIB_OK)
    {
        return 0;
    }
    SmBrep *p = primitive(in.byte());
    if (p == nullptr)
    {
        return 0;
    }
    smlib_brep *brep = new smlib_brep{p};
    // Move some vertices by arbitrary offsets.
    SmTArray<SmVertex *> vertices;
    p->GetVertices(vertices);
    for (uint8_t moves = in.byte() % 4; moves > 0 && vertices.GetSize() > 0; --moves)
    {
        SmVertex *v = vertices[in.byte() % vertices.GetSize()];
        SmPoint3d q = v->GetPoint();
        q.x += in.real();
        q.y += in.real();
        q.z += in.real();
        v->SetPoint(q);
    }
    char buf[256];
    while (in.size > 0)
    {
        int32_t level = static_cast<int32_t>(in.integer() % 5);
        int64_t cap = in.integer() % 300;
        uint8_t flags = in.byte();
        int64_t len = 0, failures = 0;
        smlib_brep_check(flags & 1 ? brep : nullptr, level, flags & 2 ? buf : nullptr, cap > 256 ? 256 : cap,
                         flags & 4 ? &len : nullptr, flags & 8 ? &failures : nullptr);
    }
    smlib_brep_free(brep);
    return 0;
}
