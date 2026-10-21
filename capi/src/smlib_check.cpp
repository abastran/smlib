// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// smlib_brep_check: SMLib's own validity checks (SmBrep::ValidatePointers, then AssertValid over
// the whole topology graph), with every failed test as one line of text.

#include "smlib_shim.h"

#include <SmAssertArray.h>

using namespace smlib_shim;

// Defined in SmAssertArray.cpp, which declares it in no header.
TCHAR const *SM_GetObjTypeString(SM_TYPE lExeClassType);

namespace
{

static_assert(sizeof(TCHAR) == 1, "SMLib's report strings are taken as UTF-8 bytes");

const char *text(const TCHAR *s)
{
    return s != nullptr && s[0] != '\0' ? s : "?";
}

// `<reporting type>::Test #<n> <name>: <message, offset, tolerance> [owner <type>]`
std::string line(const SmAssertReport &r)
{
    char msg[SM_TBLOCK_SIZE];
    r.FormatLogMessage(msg, sizeof msg);
    msg[sizeof msg - 1] = '\0';
    char head[128];
    std::snprintf(head, sizeof head, "%s::Test #%lu ", text(SM_GetObjTypeString(r.m_lReportingType)),
                  static_cast<unsigned long>(r.m_lTestIndex));
    std::string s = std::string(head) + text(r.m_pName) + ": " + text(msg);
    if (r.m_sOwnerTypeString != nullptr && r.m_sOwnerTypeString[0] != '\0')
    {
        s += std::string(" [owner ") + r.m_sOwnerTypeString + "]";
    }
    return s;
}

// Appends the failed reports in `list` to `out`, counting them.
void collect(const SmAssertArray &list, std::string &out, int64_t &failures)
{
    for (ULONG i = 0; i < list.GetSize(); ++i)
    {
        const SmAssertReport *r = list[i];
        if (r == nullptr || r->m_bOK)
        {
            continue;
        }
        ++failures;
        out += line(*r);
        out += '\n';
    }
}

} // namespace

extern "C" int32_t smlib_brep_check(const smlib_brep *brep, int32_t level, char *buf, int64_t cap, int64_t *out_len,
                                    int64_t *out_failures)
{
    if (!valid(brep) || level < 0 || level > 2 || cap < 0 || (buf == nullptr && cap != 0) || out_len == nullptr ||
        out_failures == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_brep_check", handle(brep) + ", level=" + std::to_string(level) + ", cap=" + std::to_string(cap));
    c.input("brep", brep);
    return c.run(
        [&]
        {
            const SmAssertTestLevel lv = static_cast<SmAssertTestLevel>(level);
            std::string report;
            int64_t failures = 0;
            {
                SmAssertArray pointers;
                SmStatus s = brep->p->ValidatePointers(&pointers, lv);
                collect(pointers, report, failures);
                if (s != SM_SUCCESS && failures == 0)
                {
                    ++failures;
                    report += "SmBrep::ValidatePointers failed with status " + std::to_string(static_cast<int>(s)) + "\n";
                }
            }
            if (failures > 0)
            {
                // AssertValid walks the pointers ValidatePointers just found broken.
                report += "AssertValid skipped: ValidatePointers failed\n";
            }
            else
            {
                SmAssertArray asserts;
                SmBoolean ok = brep->p->AssertValid(&asserts, lv, SM_WALK);
                collect(asserts, report, failures);
                if (!ok && failures == 0)
                {
                    ++failures;
                    report += "SmBrep::AssertValid failed without a report\n";
                }
            }
            const int64_t n = static_cast<int64_t>(report.size());
            *out_len = n;
            if (n >= cap)
            {
                return SMLIB_ERR_BUFFER;
            }
            std::memcpy(buf, report.data(), report.size());
            buf[n] = '\0';
            *out_failures = failures;
            return SMLIB_OK;
        }
    );
}
