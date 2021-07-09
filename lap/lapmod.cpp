#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lapjv.h"

// #define ORIGINAL

/** Column-reduction and reduction transfer for a sparse cost matrix.
 */
int_t _ccrrt_sparse(const int_t n, cost_t *cc, int_t *ii, int_t *kk,
                      int_t *free_rows, int_t *x, int_t *y, cost_t *v,
                      cost_t large)
{
    int_t n_free_rows;
    boolean *unique;

    for (int_t i = 0; i < n; i++) {
        x[i] = -1;
        v[i] = large;
        y[i] = 0;
    }
    for (int_t i = 0; i < n; i++) {
        for (int_t k = ii[i]; k < ii[i+1]; k++) {
            const int_t j = kk[k];
            const cost_t c = cc[k];
            if (c < v[j]) {
                v[j] = c;
                y[j] = i;
            }
            PRINTF("i=%d, k=%d, j=%d, c[i,j]=%f, v[j]=%f y[j]=%d\n", i, k, j, c, v[j], y[j]);
        }
    }
    PRINT_COST_ARRAY(v, n);
    PRINT_INDEX_ARRAY(y, n);
    NEW(unique, boolean, n);
    memset(unique, TRUE, n);
    {
        int_t j = n;
        do {
            j--;
            const int_t i = y[j];
            if (x[i] < 0) {
                x[i] = j;
            } else {
                unique[i] = FALSE;
                y[j] = -1;
            }
        } while (j > 0);
    }
    n_free_rows = 0;
    for (int_t i = 0; i < n; i++) {
        if (x[i] < 0) {
            free_rows[n_free_rows++] = i;
        } else if (unique[i]) {
            const int_t j = x[i];
            cost_t min = LARGE;
            
            // <<
            // for (uint_t k = ii[i]; k < ii[i+1]; k++) {
            //     const int_t j2 = kk[k];
            //     if (j2 == j) continue;
                
            //     const cost_t c = cc[k] - v[j2];
            //     if(c < min)
            //         min = c;
            // }
            // --
            int_t k    = ii[i];
            cost_t cj2 = 0;
            
            for (int_t j2 = 0; j2 < n; j2++) {
                
                if(k < ii[i + 1] && kk[k] == j2) {
                    cj2 = cc[k];
                    k++;
                } else {
                    cj2 = large;
                }
                
                if (j2 == (int_t)j) continue;
                
                const cost_t c = cj2 - v[j2];
                if (c < min)
                    min = c;
            }
            // >>
            
            v[j] -= min;
        }
    }
    FREE(unique);
    return n_free_rows;
}


/** Augmenting row reduction for a sparse cost matrix.
 */
int_t _carr_sparse(
    const int_t n, cost_t *cc, int_t *ii, int_t *kk,
    const int_t n_free_rows,
    int_t *free_rows, int_t *x, int_t *y, cost_t *v,
    cost_t large)
{
    int_t current = 0;
    int_t new_free_rows = 0;
    PRINT_INDEX_ARRAY(x, n);
    PRINT_INDEX_ARRAY(y, n);
    PRINT_COST_ARRAY(v, n);
    PRINT_INDEX_ARRAY(free_rows, n_free_rows);
    while (current < n_free_rows) {
        int_t i0;
        int_t j1, j2;
        cost_t v1, v2, v1_new;
        boolean v1_lowers;

        PRINTF("current = %d\n", current);
        const int_t free_i = free_rows[current++];
        PRINTF("free_i = %d\n", free_i);
        // if (ii[free_i+1] - ii[free_i] > 0) {
        //     const int_t k = ii[free_i];
        //     j1 = kk[k];
        //     v1 = cc[k] - v[j1];
        // } else {
        //     // This means the next row is empty
        //     j1 = 0;
        //     v1 = large;
        // }
        j1 = 0;
        if (ii[free_i+1] - ii[free_i] > 0 && kk[ii[free_i]] == 0) {
            const int_t k = ii[free_i];
            v1 = cc[k] - v[j1];
        } else {
            v1 = large - v[0];
        }
        j2 = -1;
        v2 = LARGE;
        // meaning this loop doesn't run
        // for (int_t k = ii[free_i]+1; k < ii[free_i+1]; k++) {
        //     PRINTF("%d = %f %d = %f\n", j1, v1, j2, v2);
        //     const int_t j = kk[k];
        //     const cost_t c = cc[k] - v[j];

        int_t k = ii[free_i];
        if(kk[k] == 0) k++;
        
        // << HOTSPOT

#ifdef ORIGINAL
        cost_t cj = 0;
        for (int_t j = 1; j < n; j++) {
            
            if(k < ii[free_i + 1] && kk[k] == j) {
                cj = cc[k];
                k++;
            } else {
                cj = large;
            }
            
            const cost_t c = cj - v[j];
            
            if (c < v2) {
                if (c >= v1) {
                    v2 = c;
                    j2 = j;
                } else {
                    v2 = v1;
                    v1 = c;
                    j2 = j1;
                    j1 = j;
                }
            }
        }
#else
        for (int_t k = ii[free_i]; k < ii[free_i + 1]; k++) {
            int_t j        = kk[k];            
            const cost_t c = cc[k] - v[j]; // find smallest values of cc[k] - v[j]
            
            if (c < v2) {
                if (c >= v1) {
                    v2 = c;
                    j2 = j;
                } else {
                    v2 = v1;
                    v1 = c;
                    j2 = j1;
                    j1 = j;
                }
            }
        }
                
        for (int_t j = 1; j < n; j++) {
            const cost_t c = large - v[j]; // find largest values of v[j]
            if (c < v2) {
                if (c >= v1) {
                    v2 = c;
                    j2 = j;
                } else {
                    v2 = v1;
                    v1 = c;
                    j2 = j1;
                    j1 = j;
                }
            }
        }
#endif
        
        i0 = y[j1];
        v1_new = v[j1] - (v2 - v1);
        v1_lowers = v1_new < v[j1];

        if (v1_lowers) {
            v[j1] = v1_new;
        } else if (i0 >= 0 && j2 >= 0) {
            j1 = j2;
            i0 = y[j2];
        }
        x[free_i] = j1;
        y[j1] = free_i;
        if (i0 >= 0) {
            if (v1_lowers) {
                free_rows[--current] = i0;
            } else {
                free_rows[new_free_rows++] = i0;
            }
        }
    }
    return new_free_rows;
}


/** Find columns with minimum d[j] and put them on the SCAN list.
 */
int_t _find_sparse_1(const int_t n, int_t lo, cost_t *d, int_t *cols, int_t *y)
{
    
    int_t hi    = lo + 1;
    cost_t mind = d[cols[lo]];
    
    for (int_t k = hi; k < n; k++) {
        int_t j = cols[k];
        if (d[j] <= mind) {
            if (d[j] < mind) {
                hi   = lo;
                mind = d[j];
            }
            cols[k]    = cols[hi];
            cols[hi++] = j;
        }
    }
    return hi;
}

/** Scan all columns in TODO starting from arbitrary column in SCAN and try to
 * decrease d of the TODO columns using the SCAN column.
 */
int_t _scan_sparse_1(
    const int_t n, cost_t *cc, int_t *ii, int_t *kk,
    int_t *plo, int_t *phi,
    cost_t *d, int_t *cols, int_t *pred,
    int_t *y, cost_t *v, cost_t large, int_t* rev_kk)
{
    
    int_t lo = *plo;
    int_t hi = *phi;
    cost_t h, cred_ij;

    while (lo != hi) {
        int_t kj;
        int_t j           = cols[lo++];
        const int_t i     = y[j];
        const cost_t mind = d[j];
        
        for (int_t k = 0; k < n; k++) rev_kk[k] = -1;
        
        for (int_t k = ii[i]; k < ii[i+1]; k++) {
            const int_t j = kk[k];
            rev_kk[j] = k;
        }

        kj = rev_kk[j];
        if (kj == -1) {
            h = large - v[j] - mind;
        } else {
            h = cc[kj] - v[j] - mind;
        }
        
        // For all columns in TODO
        for (int_t k = hi; k < n; k++) {
            j = cols[k];

            if ((kj = rev_kk[j]) == -1) {
                cred_ij = large - v[j] - h;
            } else {
                cred_ij = cc[kj] - v[j] - h;
            }
            
            if (cred_ij < d[j]) {
                d[j]    = cred_ij;
                pred[j] = i;
                
                if (cred_ij == mind) {
                    if (y[j] < 0) {
                        return j;
                    }
                    
                    cols[k] = cols[hi];
                    cols[hi++] = j;
                }
            }
        }
    }
    *plo = lo;
    *phi = hi;
    
    return -1;
}

/** Single iteration of modified Dijkstra shortest path algorithm as explained in the JV paper.
 *
 * This version loops over all column indices (some of which might be inf).
 *
 * \return The closest free column index.
 */
__inline__ int_t find_path_sparse_1(
    const int_t n,
    cost_t *cc,
    int_t *ii,
    int_t *kk,
    const int_t start_i,
    int_t *y,
    cost_t *v,
    int_t *pred,
    cost_t large,
    int_t* cols,
    cost_t* d
) {
    int_t lo = 0, hi = 0;
    int_t final_j = -1;
    int_t n_ready = 0;

#ifdef ORIGINAL
    for (int_t i = 0; i < n; i++) {
        cols[i] = i;       // range(n)
        d[i]    = large;   // [large] * n
        pred[i] = start_i; // [start_i] * n
    }

    // <<
    // for (int_t i = ii[start_i]; i < ii[start_i + 1]; i++) {
    //     const int_t j = kk[i];
    //     d[j] = cc[i] - v[j];
    // }
    // --
    cost_t cci = 0;
    int_t k    = ii[start_i];
    for(int_t j = 0; j < n; j++) {
        
        if(k < ii[start_i + 1] && j == kk[k]) {
            cci = cc[k];
            k++;
        } else {
            cci = large;
        }
        
        d[j] = cci - v[j];
    }
    // >>
#else
    for (int_t i = 0; i < n; i++) {
        cols[i] = i;
        d[i]    = large - v[i];
        pred[i] = start_i;
    }
    for (int_t k = ii[start_i]; k < ii[start_i + 1]; k++) {
        int_t j = kk[k];
        d[j]    = cc[k] - v[j];
    }
#endif

    int_t *rev_kk;
    NEW(rev_kk, int_t, n);
    
    while (final_j == -1) {
        if (lo == hi) {

            n_ready = lo;
            hi      = _find_sparse_1(n, lo, d, cols, y);
            
            for (int_t k = lo; k < hi; k++) {
                const int_t j = cols[k];
                if (y[j] < 0) {
                    final_j = j;
                }
            }
        }
        
        if (final_j == -1) {
            final_j = _scan_sparse_1(n, cc, ii, kk, &lo, &hi, d, cols, pred, y, v, large, rev_kk);
        }
    }
    
    const cost_t mind = d[cols[lo]];
    for (int_t k = 0; k < n_ready; k++) {
        const int_t j = cols[k];
        v[j] += d[j] - mind;
    }

    FREE(rev_kk);
    return final_j;
}

typedef int_t (*fp_function_t)(
        const int_t, cost_t *, int_t *, int_t *, const int_t, int_t *, cost_t *, int_t *, cost_t large);


/** Augment for a sparse cost matrix. */
int_t _ca_sparse(
    const int_t n, cost_t *cc, int_t *ii, int_t *kk,
    const int_t n_free_rows,
    int_t *free_rows, int_t *x, int_t *y, cost_t *v,
    int fp_version, cost_t large)
{
    int_t *pred;
    int_t *cols;
    cost_t *d;
    NEW(pred, int_t, n);
    NEW(cols, int_t, n);
    NEW(d,    cost_t, n);
    
    for (int_t *pfree_i = free_rows; pfree_i < free_rows + n_free_rows; pfree_i++) {
        int_t i = -1;
        int_t k = 0;
        
        int_t j = find_path_sparse_1(n, cc, ii, kk, *pfree_i, y, v, pred, large, cols, d);
        
        while (i != *pfree_i) {
            i    = pred[j];
            y[j] = i;
            SWAP_INDICES(j, x[i]);
            k++;
        }
        
    }
    FREE(cols);
    FREE(d);
    FREE(pred);
    return 0;
}


/** Solve square sparse LAP.
 */
int_t lapmod_internal(
    const int_t n,
    cost_t *cc,
    int_t *ii,
    int_t *kk,
    int_t *x,
    int_t *y,
    fp_t fp_version,
    cost_t large
) {

#ifdef ORIGINAL
    printf("original\n");
#else
    printf("optimized\n");
#endif

    int_t ret;
    int_t *free_rows;
    cost_t *v;

    NEW(free_rows, int_t, n);
    NEW(v, cost_t, n);
    
    printf("_ccrrt_sparse\n");
    ret = _ccrrt_sparse(n, cc, ii, kk, free_rows, x, y, v, large);
    int i = 0;
    
    while (ret > 0 && i < 2) {
        printf("_carr_sparse\n");
        ret = _carr_sparse(n, cc, ii, kk, ret, free_rows, x, y, v, large);
        i++;
    }
    
    if (ret > 0) {
        printf("_ca_sparse\n");
        ret = _ca_sparse(n, cc, ii, kk, ret, free_rows, x, y, v, fp_version, large);
    }
    
    FREE(v);
    FREE(free_rows);
    return ret;
}
