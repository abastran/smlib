// SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors
// SPDX-License-Identifier: Apache-2.0

// The flat B-rep arrays' entry points (include/smlib.h, "Flat B-rep arrays"): the opaque
// smlib_flat, array access by id, export, import and the text form used by failure dumps. The
// conversion itself is src/smlib_flat.cpp.

#include "smlib_flat.h"
#include "smlib_shim.h"

#include <SmApiQueries.h>

using namespace smlib_shim;

struct smlib_flat
{
    smlib_flat_impl::FlatBrep f;
};

namespace
{

std::string handle(const smlib_flat *f)
{
    char t[32];
    std::snprintf(t, sizeof t, "%p", static_cast<const void *>(f));
    return f == nullptr ? std::string("null") : std::string("flat@") + t;
}

// Dumps `f` in smlib_flat_write's text form when the call fails.
void input(Call &c, const char *name, const smlib_flat *f)
{
    if (f != nullptr)
    {
        const smlib_flat_impl::FlatBrep *p = &f->f;
        c.input_text(name, [p](const std::string &path) { smlib_flat_impl::WriteText(*p, path); });
    }
}

} // namespace

extern "C" int32_t smlib_brep_is_manifold_solid(const smlib_brep *brep, int32_t *out)
{
    if (!valid(brep) || out == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_brep_is_manifold_solid", handle(brep));
    c.input("brep", brep);
    return c.run(
        [&]
        {
            SmBoolean m = FALSE;
            int32_t s = status(SmApiBrepIsManifoldSolid(brep->p, m));
            if (s == SMLIB_OK)
            {
                *out = m == TRUE ? 1 : 0;
            }
            return s;
        }
    );
}

/* ---- Flat B-rep arrays ------------------------------------------------------------------------ */

namespace
{

// The member array for `id` when it holds `T`, else nullptr.
template <typename T>
std::vector<T> *array_of(smlib_flat_impl::FlatBrep &f, int32_t id);

template <>
std::vector<int64_t> *array_of<int64_t>(smlib_flat_impl::FlatBrep &f, int32_t id)
{
    const smlib_flat_impl::Field *fd = smlib_flat_impl::FindField(id);
    return fd != nullptr && fd->i != nullptr ? &(f.*(fd->i)) : nullptr;
}

template <>
std::vector<double> *array_of<double>(smlib_flat_impl::FlatBrep &f, int32_t id)
{
    const smlib_flat_impl::Field *fd = smlib_flat_impl::FindField(id);
    return fd != nullptr && fd->f != nullptr ? &(f.*(fd->f)) : nullptr;
}

template <typename T>
int32_t flat_get(const char *fn, const smlib_flat *flat, int32_t id, const T **data, int64_t *len)
{
    if (flat == nullptr || data == nullptr || len == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    return Call(fn, handle(flat) + ", id=" + std::to_string(id)).run(
        [&]
        {
            std::vector<T> *v = array_of<T>(const_cast<smlib_flat *>(flat)->f, id);
            if (v == nullptr)
            {
                return SMLIB_ERR_ARG;
            }
            *data = v->empty() ? nullptr : v->data();
            *len = static_cast<int64_t>(v->size());
            return SMLIB_OK;
        }
    );
}

template <typename T>
int32_t flat_set(const char *fn, smlib_flat *flat, int32_t id, const T *data, int64_t len)
{
    if (flat == nullptr || len < 0 || len > SMLIB_FLAT_MAX_LEN || (data == nullptr && len > 0))
    {
        return SMLIB_ERR_ARG;
    }
    return Call(fn, handle(flat) + ", id=" + std::to_string(id) + ", len=" + std::to_string(len)).run(
        [&]
        {
            std::vector<T> *v = array_of<T>(flat->f, id);
            if (v == nullptr)
            {
                return SMLIB_ERR_ARG;
            }
            v->assign(data, data + len);
            return SMLIB_OK;
        }
    );
}

} // namespace

extern "C" smlib_flat *smlib_flat_new(void)
{
    try
    {
        return new smlib_flat();
    }
    catch (...)
    {
        return nullptr;
    }
}

extern "C" void smlib_flat_free(smlib_flat *flat)
{
    delete flat;
}

extern "C" int32_t smlib_flat_get_i64(const smlib_flat *flat, int32_t id, const int64_t **data, int64_t *len)
{
    return flat_get<int64_t>("smlib_flat_get_i64", flat, id, data, len);
}

extern "C" int32_t smlib_flat_get_f64(const smlib_flat *flat, int32_t id, const double **data, int64_t *len)
{
    return flat_get<double>("smlib_flat_get_f64", flat, id, data, len);
}

extern "C" int32_t smlib_flat_set_i64(smlib_flat *flat, int32_t id, const int64_t *data, int64_t len)
{
    return flat_set<int64_t>("smlib_flat_set_i64", flat, id, data, len);
}

extern "C" int32_t smlib_flat_set_f64(smlib_flat *flat, int32_t id, const double *data, int64_t len)
{
    return flat_set<double>("smlib_flat_set_f64", flat, id, data, len);
}

extern "C" int32_t smlib_brep_export(const smlib_brep *brep, smlib_flat **out)
{
    if (!valid(brep) || out == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_brep_export", handle(brep));
    c.input("brep", brep);
    return c.run(
        [&]
        {
            smlib_flat *r = new smlib_flat();
            int32_t s = status(smlib_flat_impl::FlatFromBrep(*brep->p, r->f));
            if (s != SMLIB_OK)
            {
                delete r;
                return s;
            }
            *out = r;
            return SMLIB_OK;
        }
    );
}

extern "C" int32_t smlib_brep_import(const smlib_flat *flat, smlib_brep **out)
{
    if (flat == nullptr || out == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_brep_import", handle(flat));
    input(c, "flat", flat);
    return c.run(
        [&]
        {
            SmBrep *r = nullptr;
            int32_t s = status(smlib_flat_impl::BrepFromFlat(flat->f, *SmApiGetOrCreateContext(), r));
            if (s != SMLIB_OK)
            {
                delete r;
                return s;
            }
            return take(r, out);
        }
    );
}

extern "C" int32_t smlib_flat_write(const smlib_flat *flat, const char *path)
{
    if (flat == nullptr || path == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    Call c("smlib_flat_write", handle(flat) + ", path=" + path);
    return c.run([&] { return smlib_flat_impl::WriteText(flat->f, path) ? SMLIB_OK : static_cast<int32_t>(SM_ERR); });
}

extern "C" int32_t smlib_flat_read(const char *path, smlib_flat **out)
{
    if (path == nullptr || out == nullptr)
    {
        return SMLIB_ERR_ARG;
    }
    return Call("smlib_flat_read", std::string("path=") + path).run(
        [&]
        {
            smlib_flat *r = new smlib_flat();
            if (!smlib_flat_impl::ReadText(path, r->f))
            {
                delete r;
                return static_cast<int32_t>(SM_ERR_INVALID_INPUT);
            }
            *out = r;
            return SMLIB_OK;
        }
    );
}
