/* SPDX-FileCopyrightText: Copyright (c) 2026 SolidCraft contributors */
/* SPDX-License-Identifier: Apache-2.0 */

/* Links the static library from plain C and checks a box end to end, plus the error paths. */

#include "smlib.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;

static int log_lines = 0;
static int log_failures = 0;
static char in_flight[256];

static void on_log(void *user, int32_t level, const char *line)
{
    (void)user;
    log_lines++;
    log_failures += level != 0;
    /* While the callback runs, the call is still the one in flight. */
    smlib_current_call(in_flight, sizeof in_flight);
    (void)line;
}

#define CHECK(cond)                                                                                \
    do                                                                                             \
    {                                                                                              \
        if (!(cond))                                                                               \
        {                                                                                          \
            fprintf(stderr, "%s:%d: CHECK(%s) failed\n", __FILE__, __LINE__, #cond);               \
            failures++;                                                                            \
        }                                                                                          \
    } while (0)

int main(void)
{
    CHECK(smlib_abi_version() == SMLIB_ABI_VERSION);
    CHECK(smlib_init() == SMLIB_OK);

    const double origin[3] = {1.0, 2.0, 3.0};
    const double size[3] = {10.0, 20.0, 30.0};
    smlib_brep *box = NULL;
    CHECK(smlib_make_box(origin, size, &box) == SMLIB_OK);
    CHECK(box != NULL);

    double volume = 0.0;
    CHECK(smlib_brep_volume(box, 1e-7, &volume) == SMLIB_OK);
    CHECK(fabs(volume - 6000.0) < 1e-6);

    int64_t counts[3] = {0, 0, 0};
    CHECK(smlib_brep_counts(box, counts) == SMLIB_OK);
    CHECK(counts[0] == 6 && counts[1] == 12 && counts[2] == 8);

    smlib_brep *copy = NULL;
    CHECK(smlib_brep_copy(box, &copy) == SMLIB_OK);
    smlib_brep_free(box);
    volume = 0.0;
    CHECK(smlib_brep_volume(copy, 1e-7, &volume) == SMLIB_OK);
    CHECK(fabs(volume - 6000.0) < 1e-6);
    smlib_brep_free(copy);

    /* Rejected arguments never reach SMLib and leave the output alone. */
    const double bad[3] = {10.0, -1.0, 30.0};
    const double nan_size[3] = {NAN, 1.0, 1.0};
    smlib_brep *none = NULL;
    CHECK(smlib_make_box(origin, bad, &none) == SMLIB_ERR_ARG && none == NULL);
    CHECK(smlib_make_box(origin, nan_size, &none) == SMLIB_ERR_ARG && none == NULL);
    CHECK(smlib_make_box(NULL, size, &none) == SMLIB_ERR_ARG);
    CHECK(smlib_brep_volume(NULL, 1e-7, &volume) == SMLIB_ERR_ARG);
    CHECK(smlib_brep_counts(NULL, counts) == SMLIB_ERR_ARG);
    smlib_brep_free(NULL);

    /* Boxes SMLib cannot make, so far from the origin for their size that their faces come out
       degenerate, fail with no body and leave nothing behind (tangency-labs/solidcore#377; LSan
       finds a leak in an ASan build, `leaks --atExit -- smlib_smoke` on macOS). */
    const double far_origins[3][3] = {{0.0, 0.0, 1e200}, {1e95, 0.0, 0.0}, {1e294, 1e294, 1e294}};
    const double far_sizes[3][3] = {{1.0, 1.0, 1.0}, {0.5, 0.5, 0.5}, {1e-3, 1e-3, 1e-3}};
    for (int i = 0; i < 3; ++i)
    {
        CHECK(smlib_make_box(far_origins[i], far_sizes[i], &none) != SMLIB_OK && none == NULL);
    }

    /* A box far smaller than the tolerance (smlib_fuzz's leak-000bffce): its top face closes up
       on the bottom one and SMLib's test of the closed shell fails (SM_ERR, 1001) after the
       face is made. The face keeps its plane, so freeing the failed body frees the plane too. */
    const double tiny_origin[3] = {0x1.e3e3e2a292902p-28, 0x1.e3e3e3e3e3e3ep-28, 0x1.e3e3e3e3e3e3ep-28};
    const double tiny_size[3] = {0x1.e3e3e3e3e3e3cp-28, 0x1.fe0f5053e3e3ep-16, 0x1.e3f3e3e3e3e3ep-28};
    CHECK(smlib_make_box(tiny_origin, tiny_size, &none) == 1001 && none == NULL);

    /* One log line per call; the call in flight is visible to a crash handler. */
    smlib_set_log(on_log, NULL);
    smlib_brep *logged = NULL;
    CHECK(smlib_make_box(origin, size, &logged) == SMLIB_OK);
    CHECK(log_lines == 1 && log_failures == 0);
    CHECK(strncmp(in_flight, "smlib_make_box(origin=[1,2,3], size=[10,20,30])", 48) == 0);
    char none_now[8];
    CHECK(smlib_current_call(none_now, sizeof none_now) == 0);

    /* A failing call dumps its line and its B-rep inputs; the dumped B-rep reads back. */
    const char *dir = getenv("SMLIB_SMOKE_DUMP_DIR");
    if (dir != NULL)
    {
        CHECK(smlib_set_dump_dir(dir) == SMLIB_OK);
        CHECK(smlib_brep_write(logged, "/nonexistent-dir/x.smb", 0) != SMLIB_OK);
        CHECK(log_failures == 1);
        char path[1024];
        snprintf(path, sizeof path, "%s/1-smlib_brep_write-brep.smb", dir);
        smlib_brep *replayed = NULL;
        CHECK(smlib_brep_read(path, 0, &replayed) == SMLIB_OK);
        volume = 0.0;
        CHECK(smlib_brep_volume(replayed, 1e-7, &volume) == SMLIB_OK);
        CHECK(fabs(volume - 6000.0) < 1e-6);
        smlib_brep_free(replayed);
        CHECK(smlib_set_dump_dir(NULL) == SMLIB_OK);
    }
    smlib_set_log(NULL, NULL);
    smlib_brep_free(logged);

    /* ---- Primitives ----------------------------------------------------------------------- */
    {
        const double pi = 3.14159265358979323846;
        const double at[3] = {0.0, 0.0, 0.0};
        smlib_brep *b = NULL;
        CHECK(smlib_make_cylinder(at, 5.0, 10.0, &b) == SMLIB_OK);
        CHECK(smlib_brep_volume(b, 1e-7, &volume) == SMLIB_OK && fabs(volume - pi * 250.0) < 1e-4);
        smlib_brep_free(b);
        b = NULL;
        CHECK(smlib_make_sphere(at, 5.0, &b) == SMLIB_OK);
        CHECK(smlib_brep_volume(b, 1e-7, &volume) == SMLIB_OK && fabs(volume - 4.0 / 3.0 * pi * 125.0) < 1e-3);
        smlib_brep_free(b);
        b = NULL;
        CHECK(smlib_make_cone(at, 5.0, 2.0, 10.0, &b) == SMLIB_OK);
        CHECK(smlib_brep_volume(b, 1e-7, &volume) == SMLIB_OK && fabs(volume - pi * 10.0 / 3.0 * 39.0) < 1e-3);
        smlib_brep_free(b);
        b = NULL;
        CHECK(smlib_make_torus(at, 10.0, 3.0, &b) == SMLIB_OK);
        CHECK(smlib_brep_volume(b, 1e-7, &volume) == SMLIB_OK && fabs(volume - 2.0 * pi * pi * 10.0 * 9.0) < 1e-2);
        smlib_brep_free(b);

        smlib_brep *none2 = NULL;
        CHECK(smlib_make_cylinder(at, -1.0, 10.0, &none2) == SMLIB_ERR_ARG && none2 == NULL);
        CHECK(smlib_make_sphere(at, NAN, &none2) == SMLIB_ERR_ARG);
        CHECK(smlib_make_cone(at, 5.0, 5.0, 10.0, &none2) == SMLIB_ERR_ARG);
        CHECK(smlib_make_cone(at, 5.0, -1.0, 10.0, &none2) == SMLIB_ERR_ARG);
        CHECK(smlib_make_torus(at, 3.0, 3.0, &none2) == SMLIB_ERR_ARG);
        CHECK(smlib_make_torus(NULL, 10.0, 3.0, &none2) == SMLIB_ERR_ARG && none2 == NULL);
    }

    if (failures == 0)
    {
        printf("smlib C API smoke test: ok\n");
    }
    return failures == 0 ? 0 : 1;
}
