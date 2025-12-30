// routines for wigner seitz cells
//
//
#include "wigner_seitz.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

#include "common/ELPH_timers.h"
#include "common/constants.h"
#include "common/error.h"
#include "common/kdtree/kdtree.h"
#include "common/numerical_func.h"
#include "elphC.h"

#define WS_SUPERCELL_SEARCH_SIZE 6
// we do searching from [-WS_SUPERCELL_SEARCH_SIZE, WS_SUPERCELL_SEARCH_SIZE]
//

static struct kdtree* setup_ws_tree(const ND_int* grid,
                                    const ELPH_float* lat_vecs,
                                    const ND_int ws_ssize);

static ND_int get_ws_nearest_superlat(struct kdtree* tree, double* query_pnt,
                                      const double eps, ELPH_float* pts_buf);

ND_int build_wigner_seitz_vectors(const ND_int* grid,
                                  const ELPH_float* lat_vecs, double eps,
                                  const ELPH_float* rvec_m, ND_int nrvec_m,
                                  const ELPH_float* rvec_n, ND_int nrvec_n,
                                  ND_int** ws_vecs, ND_int** ws_degen)
{
    // Note. you need to free buffer allocated outside of this functions
    // find all vectors T such that |r_m-(r_n+R+T)|(or |(r_n+R+T)-r_m|) is
    // is mim for a given grid, R vecs are miller indices given by
    // "get_miller_idx" function "common/numerical_func.c"
    // R+T are the wigner seitz vectors for particular atomic position.
    //
    // Eq. 46 of G Pizzi et al 2020 J. Phys.: Condens. Matter 32 165902
    //
    // NM : we assume that order of R vecs is (X,Y,Z) in row major i.e
    // Rz variying fast followed by Ry and then Rx
    //
    // rvec_m and rvec_n are in cart units
    // must free ws_vecs and ws_degen outside of the function
    // ws_degen has a size = grid[0]*grid[1]*grid[2]*nrvec_n*nrvec_m with
    // ws_degen[i] (Rx,Ry,Rz,nrvec_m, nrvec_n) give the degeneracy of R point.
    // ws_vecs has as size = 3*sum(ws_degen)
    // return sum(ws_degen)
    //
    ELPH_start_clock("Wigner seitz");

    struct kdtree* tree =
        setup_ws_tree(grid, lat_vecs, WS_SUPERCELL_SEARCH_SIZE);
    if (!tree)
    {
        error_msg("Buildling KD Tree for wigner seitz cell failed.");
    }

    ELPH_float blat[9];
    reciprocal_vecs(lat_vecs, blat);
    // remove 2*pi
    for (ND_int xi = 0; xi < 9; ++xi)
    {
        blat[xi] /= (2 * ELPH_PI);
    }

    const ELPH_float zerovec[3] = {0.0, 0.0, 0.0};
    // create a tmp buffer to read points.
    size_t tree_count = 2 * WS_SUPERCELL_SEARCH_SIZE + 1;
    tree_count = tree_count * tree_count * tree_count;
    ELPH_float* pts_buf = malloc(3 * tree_count * sizeof(*pts_buf));
    CHECK_ALLOC(pts_buf);
    //
    //
    if (!rvec_m)
    {
        nrvec_m = 1;
        rvec_m = zerovec;
    }
    if (!rvec_n)
    {
        nrvec_n = 1;
        rvec_n = zerovec;
    }
    //
    const ND_int nRpts = grid[0] * grid[1] * grid[2];
    const ND_int nRmnpts = nRpts * nrvec_m * nrvec_n;
    //
    // intially allocated 2*nRmnpts ws vectors. if more
    // needed, we can realloc.
    ND_int ws_vec_size = nRmnpts * 2;

    ND_int* ws_vec_buf = malloc(3 * ws_vec_size * sizeof(*ws_vec_buf));
    CHECK_ALLOC(ws_vec_buf);

    ND_int* ws_degen_buf = calloc(nRmnpts, sizeof(*ws_degen_buf));
    CHECK_ALLOC(ws_degen_buf);
    *ws_degen = ws_degen_buf;
    //
    ND_int nws_vec_found = 0;
    //
    ND_int Gridyz = grid[1] * grid[2];
    //
    for (ND_int i = 0; i < nRpts; ++i)
    {
        ND_int Rx = i / Gridyz;
        ND_int Ry = (i % Gridyz) / grid[2];
        ND_int Rz = (i % Gridyz) % grid[2];
        //
        Rx = get_miller_idx(Rx, grid[0]);
        Ry = get_miller_idx(Ry, grid[1]);
        Rz = get_miller_idx(Rz, grid[2]);
        //
        ELPH_float Rvec[3];
        Rvec[0] = lat_vecs[0] * Rx + lat_vecs[1] * Ry + lat_vecs[2] * Rz;
        Rvec[1] = lat_vecs[3] * Rx + lat_vecs[4] * Ry + lat_vecs[5] * Rz;
        Rvec[2] = lat_vecs[6] * Rx + lat_vecs[7] * Ry + lat_vecs[8] * Rz;
        //
        for (ND_int im = 0; im < nrvec_m; ++im)
        {
            for (ND_int in = 0; in < nrvec_n; ++in)
            {
                ND_int ishift = in + im * nrvec_n + i * nrvec_m * nrvec_n;
                double query_pnt[3];
                query_pnt[0] =
                    rvec_m[3 * im + 0] - rvec_n[3 * in + 0] - Rvec[0];
                query_pnt[1] =
                    rvec_m[3 * im + 1] - rvec_n[3 * in + 1] - Rvec[1];
                query_pnt[2] =
                    rvec_m[3 * im + 2] - rvec_n[3 * in + 2] - Rvec[2];
                // find the norm of query point
                double query_norm2 = sqrt(dot3_macro(query_pnt, query_pnt));
                //
                ND_int i_ws_found = get_ws_nearest_superlat(
                    tree, query_pnt, eps * query_norm2 + 1e-8, pts_buf);
                if (0 == i_ws_found)
                {
                    error_msg("No Wigner Seitz vector found.");
                }
                //
                ws_degen_buf[ishift] = i_ws_found;
                if (nws_vec_found + i_ws_found > ws_vec_size)
                {
                    // realloc.
                    ws_vec_size = ws_vec_size + (MAX(nRmnpts, i_ws_found));
                    ND_int* realloc_ptr = realloc(
                        ws_vec_buf, 3 * ws_vec_size * sizeof(*ws_vec_buf));
                    if (!realloc_ptr)
                    {
                        free(ws_vec_buf);
                    }
                    CHECK_ALLOC(realloc_ptr);
                    ws_vec_buf = realloc_ptr;
                }
                // Copy the wigner seitz vectors.
                for (ND_int ii = 0; ii < i_ws_found; ++ii)
                {
                    //
                    ELPH_float tmp_out[3];
                    tmp_out[0] = blat[0] * pts_buf[3 * ii] +
                                 blat[3] * pts_buf[3 * ii + 1] +
                                 blat[6] * pts_buf[3 * ii + 2];
                    tmp_out[1] = blat[1] * pts_buf[3 * ii] +
                                 blat[4] * pts_buf[3 * ii + 1] +
                                 blat[7] * pts_buf[3 * ii + 2];
                    tmp_out[2] = blat[2] * pts_buf[3 * ii] +
                                 blat[5] * pts_buf[3 * ii + 1] +
                                 blat[8] * pts_buf[3 * ii + 2];
                    //
                    ws_vec_buf[3 * nws_vec_found] = rint(tmp_out[0]);
                    ws_vec_buf[3 * nws_vec_found + 1] = rint(tmp_out[1]);
                    ws_vec_buf[3 * nws_vec_found + 2] = rint(tmp_out[2]);
                    // Do a sanity check
                    tmp_out[0] -= ws_vec_buf[3 * nws_vec_found];
                    tmp_out[1] -= ws_vec_buf[3 * nws_vec_found + 1];
                    tmp_out[2] -= ws_vec_buf[3 * nws_vec_found + 2];
                    if (fabs(tmp_out[0]) > 1e-3 || fabs(tmp_out[1]) > 1e-3 ||
                        fabs(tmp_out[2]) > 1e-3)
                    {
                        error_msg(
                            "Something wrong with kdtree. not a lattice "
                            "vector.");
                    }
                    //
                    ++nws_vec_found;
                }
            }
        }
    }

    // shrink the buffer
    if (nws_vec_found)
    {
        ND_int* realloc_ptr =
            realloc(ws_vec_buf, 3 * nws_vec_found * sizeof(*ws_vec_buf));
        if (realloc_ptr)
        {
            ws_vec_buf = realloc_ptr;
        }
    }
    // set the ws_vec pointer
    *ws_vecs = ws_vec_buf;
    //

    free(pts_buf);
    kd_free(tree);
    //
    ELPH_stop_clock("Wigner seitz");
    return nws_vec_found;
}

/* Static functions */

static struct kdtree* setup_ws_tree(const ND_int* grid,
                                    const ELPH_float* lat_vecs,
                                    const ND_int ws_ssize)
{
    // build a kdree for superlattice search
    struct kdtree* tree = kd_create(3);
    if (!tree)
    {
        return tree;
    }
    //
    ELPH_float superlat_vecs[9];
    // a[:,i] are latvecs
    for (ND_int i = 0; i < 3; ++i)
    {
        for (ND_int j = 0; j < 3; ++j)
        {
            superlat_vecs[j * 3 + i] = lat_vecs[j * 3 + i] * grid[i];
        }
    }
    //
    for (ND_int i = -ws_ssize; i < ws_ssize + 1; ++i)
    {
        for (ND_int j = -ws_ssize; j < ws_ssize + 1; ++j)
        {
            for (ND_int k = -ws_ssize; k < ws_ssize + 1; ++k)
            {
                double treepoint[3];
                treepoint[0] = superlat_vecs[0] * i + superlat_vecs[1] * j +
                               superlat_vecs[2] * k;
                treepoint[1] = superlat_vecs[3] * i + superlat_vecs[4] * j +
                               superlat_vecs[5] * k;
                treepoint[2] = superlat_vecs[6] * i + superlat_vecs[7] * j +
                               superlat_vecs[8] * k;
                //
                kd_insert(tree, treepoint, NULL);
            }
        }
    }
    //
    return tree;
}
static ND_int get_ws_nearest_superlat(struct kdtree* tree, double* query_pnt,
                                      const double eps, ELPH_float* pts_buf)
{
    // here we assume that the user gives
    // enough buffer to the function.
    // the size is atmost 3*get_ws_tree_size()
    // return number of point found
    // epsilon is the thresould
    // The return points will be in cart units
    if (!tree)
    {
        return 0;
    }
    // first get the nearest neighbour.
    struct kdres* nearest = kd_nearest(tree, query_pnt);
    if (!nearest)
    {
        error_msg("Nearest neighbor search failed");
    }
    //
    // Get the distance and position of the nearest neighbor
    double nearest_pos[3];
    kd_res_item(nearest, nearest_pos);
    kd_res_free(nearest);
    //
    nearest_pos[0] -= query_pnt[0];
    nearest_pos[1] -= query_pnt[1];
    nearest_pos[2] -= query_pnt[2];
    //
    double nearest_dist = sqrt(dot3_macro(nearest_pos, nearest_pos)) + 1e-10;
    // NM : we add a small absoulte tolerance so to protect eps being 0

    // now find all within radius
    struct kdres* range =
        kd_nearest_range(tree, query_pnt, nearest_dist + fabs(eps));
    if (!range)
    {
        error_msg("Cannot find wigner seitz vectors.");
    }
    //
    kd_res_rewind(range);
    ND_int Nfound = 0;
    while (!kd_res_end(range))
    {
        double pos[3];
        kd_res_item(range, pos);
        pts_buf[3 * Nfound] = pos[0];
        pts_buf[3 * Nfound + 1] = pos[1];
        pts_buf[3 * Nfound + 2] = pos[2];
        kd_res_next(range);
        ++Nfound;
    }
    kd_res_free(range);

    return Nfound;
}
