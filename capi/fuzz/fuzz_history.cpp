// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// libFuzzer target for the history shims (smlib_history_tag, smlib_history_read): two boxes
// from the input, tagged with arbitrary operand numbers (often out of range), a boolean of the
// tagged copies through SM_API (as the operation shims run it), and the result's history read
// back with arbitrary buffer sizes. Any status is fine; a crash, hang, leak or escaped exception
// is not.
//
//   cmake -S capi -B build-fuzz -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DSMLIB_CAPI_FUZZ=ON
//   cmake --build build-fuzz --target smlib_fuzz_history && build-fuzz/smlib_fuzz_history -max_total_time=60

#include "smlib_shim.h"

#include <SmApiBrep.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

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

    int32_t operand()
    {
        int32_t x = 0;
        for (int i = 0; i < 4; ++i)
        {
            x = static_cast<int32_t>(static_cast<uint32_t>(x) << 8 | byte());
        }
        // Mostly small and valid, sometimes anything.
        return byte() & 1 ? x : x & 3;
    }
};

// A box from the input; coordinates on a coarse grid half the time, so faces often coincide.
smlib_brep *box(Input &in)
{
    double o[3], s[3];
    bool grid = in.byte() & 1;
    for (int i = 0; i < 3; ++i)
    {
        o[i] = grid ? static_cast<double>(in.byte() % 8) * 2.5 : in.real();
        s[i] = grid ? static_cast<double>(1 + in.byte() % 8) * 2.5 : in.real();
    }
    smlib_brep *b = nullptr;
    smlib_make_box(o, s, &b);
    return b;
}

void read(Input &in, const smlib_brep *b)
{
    int64_t len = 0;
    std::vector<int64_t> buf(in.byte());
    smlib_history_read(in.byte() & 1 ? b : nullptr, buf.empty() ? nullptr : buf.data(), static_cast<int64_t>(buf.size()), &len);
    if (in.byte() & 1 && smlib_history_read(b, nullptr, 0, &len) == SMLIB_ERR_BUFFER && len > 0 && len < (1 << 24))
    {
        buf.resize(static_cast<size_t>(len));
        smlib_history_read(b, buf.data(), len, &len);
    }
}

} // namespace

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    Input in{data, size};
    smlib_init();
    smlib_brep *a = box(in), *b = box(in);
    smlib_brep *ta = nullptr, *tb = nullptr;
    smlib_history_tag(in.byte() & 1 ? a : nullptr, in.operand(), &ta);
    smlib_history_tag(b, in.operand(), &tb);
    // Tagging a tagged body again.
    if (ta != nullptr && in.byte() & 1)
    {
        smlib_brep *again = nullptr;
        if (smlib_history_tag(ta, in.operand(), &again) == SMLIB_OK)
        {
            smlib_brep_free(ta);
            ta = again;
        }
    }
    if (ta != nullptr && tb != nullptr)
    {
        static const SmBooleanOperationType ops[] = {SM_BO_UNION, SM_BO_DIFFERENCE, SM_BO_INTERSECTION};
        SmBrep *r = nullptr;
        SmApiStatus s = SM_ERR;
        try
        {
            s = SmApiBoolean(ta->p, tb->p, ops[in.byte() % 3], r);
        }
        catch (...)
        {
        }
        if (s == SM_SUCCESS)
        {
            // The kernel kept the first operand as the result and deleted the second.
            tb->p = nullptr;
            ta->p = r;
        }
    }
    while (in.size > 0)
    {
        read(in, in.byte() & 1 ? ta : a);
    }
    read(in, tb);
    smlib_brep_free(a);
    smlib_brep_free(b);
    smlib_brep_free(ta);
    smlib_brep_free(tb);
    return 0;
}
