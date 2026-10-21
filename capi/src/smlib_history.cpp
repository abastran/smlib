// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// Face and edge history out of SMLib operations (include/smlib.h, "History").
//
// SMLib carries identity through its edits with attributes. Every edit that changes topology
// calls SmObject::Notify on what it changes, and SmAObject::Notify hands each attribute of the
// object to a callback: Split when a face or edge is cut in two (both children get it), Merge
// when two become one (the survivor gets it), Copy when an object is copied (a boolean copies
// the second operand's faces and edges into the first), Coincident when a boolean finds a face,
// edge or vertex of one operand lying on one of the other's. New topology (a boolean's
// intersection edges) is made without any.
//
// So a Tag attribute put on each operand face and edge, holding the (operand, index) pairs it
// stands for, ends up on every result face and edge that came from it. Its callbacks join the
// pairs where SMLib's defaults would keep only one: on a merge or a coincidence the survivor
// stands for both. Reading the tags back gives each result face's origins; whether a face is
// Kept, Modified, Split or Merged follows from how many origins it has, how many result faces
// share them and whether its boundary is still exactly its origin's edges.

#include "smlib_shim.h"

#include <SmApiQueries.h>
#include <SmAttribute.h>
#include <SmEdge.h>
#include <SmFace.h>

#include <map>
#include <set>
#include <unordered_map>

using namespace smlib_shim;

namespace
{

// Above SMLib's own ids (all below SM_AI_USER_1, SmAttribute.h).
constexpr ULONG TAG_ID = SM_AI_USER_1 + 7;

// Its own type id, so a tag is told apart without RTTI (SmTypes.h: applications start here).
constexpr SM_TYPE TAG_TYPE = APPLICATION_BASE_TYPE + 4711;

// Operand numbers and indices a tag accepts (keeps every value inside a 32-bit long).
constexpr long TAG_MAX = 1L << 30;

using Origin = std::pair<long, long>;

// What a face or edge stands for: the operand entities it came from, and for a face tagged
// directly, its own boundary edges then (to tell a face that is still whole).
class Tag : public SmAttribute
{
public:
    std::set<Origin> origins;
    std::set<long> edges;

    Tag(std::set<Origin> o, std::set<long> e) : SmAttribute(TAG_ID, SM_AB_COPY), origins(std::move(o)), edges(std::move(e))
    {
    }

    SmAttribute *MakeCopy(const SmContext &c) const override
    {
        return new (c) Tag(origins, edges);
    }

    SM_TYPE GetType() const override
    {
        return TAG_TYPE;
    }
    SmBoolean IsKindOf(SM_TYPE t) const override
    {
        return t == TAG_TYPE ? TRUE : SmAttribute::IsKindOf(t);
    }

    // A split: both children stand for what the split object stood for.
    void Split(SmAObject *split, SmAObject *a, SmAObject *b) override
    {
        if (a != split)
        {
            join(a, origins, edges);
        }
        if (b != split)
        {
            join(b, origins, edges);
        }
    }

    // A merge: the result stands for both.
    void Merge(SmAObject *a, SmAObject *b, SmAObject *result) override
    {
        std::set<Origin> o;
        std::set<long> e;
        for (SmAObject *x : {a, b})
        {
            if (Tag *t = on(x); t != nullptr && x != result)
            {
                o.insert(t->origins.begin(), t->origins.end());
                e.insert(t->edges.begin(), t->edges.end());
            }
        }
        join(result, o, e);
    }

    // One operand's entity found lying on the other's, which survives: it stands for both.
    void Coincident(SmAObject *, SmAObject *to, SmAObject *, SmAObject *) override
    {
        join(to, origins, edges);
    }

    // The tag on `x`, if it has one this library made.
    static Tag *on(const SmAObject *x)
    {
        SmAttribute *a = x == nullptr ? nullptr : x->FindAttribute(TAG_ID);
        return a != nullptr && a->GetType() == TAG_TYPE ? static_cast<Tag *>(a) : nullptr;
    }

    // Adds `o` and `e` to `x`'s tag, giving it one if it has none.
    static void join(SmAObject *x, const std::set<Origin> &o, const std::set<long> &e)
    {
        if (x == nullptr || o.empty())
        {
            return;
        }
        if (Tag *t = on(x))
        {
            t->origins.insert(o.begin(), o.end());
            t->edges.insert(e.begin(), e.end());
            return;
        }
        // An attribute with this id that is no Tag (someone else's) is replaced.
        if (SmAttribute *old = x->FindAttribute(TAG_ID))
        {
            x->RemoveAttribute(old);
        }
        x->AddAttribute(new (*x->GetContext()) Tag(o, e));
    }

    // Replaces whatever tag `x` has with exactly `o` and `e`.
    static void set(SmAObject *x, std::set<Origin> o, std::set<long> e)
    {
        if (SmAttribute *old = x->FindAttribute(TAG_ID))
        {
            x->RemoveAttribute(old);
        }
        x->AddAttribute(new (*x->GetContext()) Tag(std::move(o), std::move(e)));
    }
};

// What `x`'s tag says; false when it has none.
bool read(const SmAObject *x, std::set<Origin> &o, std::set<long> &e)
{
    Tag *t = Tag::on(x);
    o = t != nullptr ? t->origins : std::set<Origin>();
    e = t != nullptr ? t->edges : std::set<long>();
    return t != nullptr;
}

// Each face's and edge's position in `b`'s face and edge order.
template <typename T>
std::unordered_map<const T *, long> positions(const SmTArray<T *> &list)
{
    std::unordered_map<const T *, long> at;
    for (ULONG i = 0; i < list.GetSize(); ++i)
    {
        at.emplace(list[i], static_cast<long>(i));
    }
    return at;
}

} // namespace

extern "C" int32_t smlib_history_tag(const smlib_brep *brep, int32_t operand, smlib_brep **out)
{
    if (!valid(brep) || out == nullptr || operand < 0 || operand >= TAG_MAX)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_history_tag", handle(brep) + ", operand=" + std::to_string(operand));
    c.input("brep", brep);
    return c.run(
        [&]
        {
            SmBrep *r = nullptr;
            int32_t s = status(SmApiBrepCopy(brep->p, r));
            if (s != SMLIB_OK)
            {
                delete r;
                return s;
            }
            SmTArray<SmFace *> faces;
            SmTArray<SmEdge *> edges;
            r->GetFaces(faces);
            r->GetEdges(edges);
            if (faces.GetSize() >= static_cast<ULONG>(TAG_MAX) || edges.GetSize() >= static_cast<ULONG>(TAG_MAX))
            {
                delete r;
                return SMLIB_ERR_ARG;
            }
            auto edge_at = positions(edges);
            for (ULONG i = 0; i < faces.GetSize(); ++i)
            {
                SmTArray<SmEdge *> bound;
                faces[i]->GetEdges(bound);
                std::set<long> own;
                for (ULONG j = 0; j < bound.GetSize(); ++j)
                {
                    auto it = edge_at.find(bound[j]);
                    if (it != edge_at.end())
                    {
                        own.insert(it->second);
                    }
                }
                Tag::set(faces[i], {{operand, static_cast<long>(i)}}, std::move(own));
            }
            for (ULONG i = 0; i < edges.GetSize(); ++i)
            {
                Tag::set(edges[i], {{operand, static_cast<long>(i)}}, {});
            }
            return take(r, out);
        }
    );
}

extern "C" int32_t smlib_history_read(const smlib_brep *result, int64_t *buf, int64_t cap, int64_t *out_len)
{
    if (!valid(result) || out_len == nullptr || cap < 0 || (buf == nullptr && cap > 0))
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_history_read", handle(result) + ", cap=" + std::to_string(cap));
    c.input("result", result);
    return c.run(
        [&]
        {
            SmTArray<SmFace *> faces;
            SmTArray<SmEdge *> edges;
            result->p->GetFaces(faces);
            result->p->GetEdges(edges);
            auto face_at = positions(faces);
            auto edge_at = positions(edges);
            ULONG nf = faces.GetSize(), ne = edges.GetSize();

            std::vector<std::set<Origin>> fo(nf), eo(ne);
            std::vector<std::set<long>> fe(nf);
            std::set<long> unused;
            std::map<Origin, long> face_uses, edge_uses;
            for (ULONG i = 0; i < nf; ++i)
            {
                read(faces[i], fo[i], fe[i]);
                for (const Origin &o : fo[i])
                {
                    ++face_uses[o];
                }
            }
            for (ULONG i = 0; i < ne; ++i)
            {
                read(edges[i], eo[i], unused);
                for (const Origin &o : eo[i])
                {
                    ++edge_uses[o];
                }
            }

            std::vector<int64_t> h{SMLIB_HISTORY_VERSION, static_cast<int64_t>(nf), static_cast<int64_t>(ne)};
            for (ULONG i = 0; i < nf; ++i)
            {
                int64_t change = SMLIB_CHANGE_GENERATED;
                if (fo[i].size() > 1)
                {
                    change = SMLIB_CHANGE_MERGED;
                }
                else if (fo[i].size() == 1)
                {
                    const Origin o = *fo[i].begin();
                    if (face_uses[o] > 1)
                    {
                        change = SMLIB_CHANGE_SPLIT;
                    }
                    else
                    {
                        // Kept when its boundary is still exactly the edges it had, each whole
                        // (one origin, under no other result edge) and its own operand's.
                        SmTArray<SmEdge *> bound;
                        faces[i]->GetEdges(bound);
                        std::set<long> now;
                        bool whole = true;
                        for (ULONG j = 0; j < bound.GetSize() && whole; ++j)
                        {
                            auto it = edge_at.find(bound[j]);
                            const std::set<Origin> *e = it == edge_at.end() ? nullptr : &eo[static_cast<size_t>(it->second)];
                            whole = e != nullptr && e->size() == 1 && e->begin()->first == o.first && edge_uses[*e->begin()] == 1;
                            if (whole)
                            {
                                now.insert(e->begin()->second);
                            }
                        }
                        change = whole && now == fe[i] ? SMLIB_CHANGE_KEPT : SMLIB_CHANGE_MODIFIED;
                    }
                }
                h.push_back(change);
                h.push_back(static_cast<int64_t>(fo[i].size()));
                for (const Origin &o : fo[i])
                {
                    h.push_back(o.first);
                    h.push_back(o.second);
                }
            }
            for (ULONG i = 0; i < ne; ++i)
            {
                SmTArray<SmFace *> on;
                edges[i]->GetFaces(on);
                std::set<long> bounds;
                for (ULONG j = 0; j < on.GetSize(); ++j)
                {
                    auto it = face_at.find(on[j]);
                    if (it != face_at.end())
                    {
                        bounds.insert(it->second);
                    }
                }
                h.push_back(static_cast<int64_t>(bounds.size()));
                h.insert(h.end(), bounds.begin(), bounds.end());
                h.push_back(static_cast<int64_t>(eo[i].size()));
                for (const Origin &o : eo[i])
                {
                    h.push_back(o.first);
                    h.push_back(o.second);
                }
            }

            *out_len = static_cast<int64_t>(h.size());
            if (static_cast<int64_t>(h.size()) > cap)
            {
                return SMLIB_ERR_BUFFER;
            }
            std::copy(h.begin(), h.end(), buf);
            return SMLIB_OK;
        }
    );
}
