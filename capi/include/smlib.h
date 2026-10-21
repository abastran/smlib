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
#define SMLIB_ABI_VERSION 1

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
 * `<n>-<function>.txt` with the call's line, and a native SMLib `<n>-<function>-<input>.smb` for
 * each B-rep it was given. NULL or "" turns this off (the default). */
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

#ifdef __cplusplus
}
#endif

#endif /* SMLIB_C_H */
