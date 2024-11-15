#include "header.h"

SEXP compute_CUSUM_and_Threshold(SEXP cumsumsI, SEXP cumsumI, SEXP tI, SEXP gridI,
                                 SEXP lenI, SEXP pI, SEXP asI, SEXP nu_asI, SEXP len_asI, SEXP estimate_meanI)
{
    PROTECT(cumsumsI);       // matrix of cumulative sums
    PROTECT(cumsumI);        // vector of last cumulative sum Y_1 + ... + Y_t
    PROTECT(tI);             // t
    PROTECT(gridI);          // the grid G^{(t)}
    PROTECT(lenI);           // length of G^{(t)}
    PROTECT(pI);             // dim p
    PROTECT(asI);            // treshold values a
    PROTECT(nu_asI);         // nu_a
    PROTECT(len_asI);        // length of the a vector
    PROTECT(estimate_meanI); // length of the a vector
    // cumsumsI: p \times len

    double *cumsums = REAL(cumsumsI);
    double *cumsum = REAL(cumsumI);
    int t = *(INTEGER(tI));
    int *grid = (INTEGER(gridI));
    int len = *(INTEGER(lenI));
    int p = *(INTEGER(pI));
    double *as = REAL(asI);
    double *nu_as = REAL(nu_asI);
    int len_as = *(INTEGER(len_asI));
    int estimate_mean = *(INTEGER(estimate_meanI));

    // array of maximum of thresholded cusum (over sparsity levels)
    // for each specific candidate changepoint location
    SEXP teststatSEXP = PROTECT(allocVector(REALSXP, len * len_as));
    double *teststat = REAL(teststatSEXP);

    for (int i = 0; i < len * len_as; ++i)
    {
        teststat[i] = 0;
    }

    // CUSUM matrix:
    SEXP CUSUMSEXP = PROTECT(allocVector(REALSXP, len * p));
    double *CUSUM = REAL(CUSUMSEXP);

    // compute CUSUM:
    int g = 0;
    for (int i = 0; i < len; ++i)
    {
        g = grid[i];

        for (int j = 0; j < p; ++j)
        {

            if (estimate_mean)
            {
                CUSUM[cord_spec(j, i, p)] = sqrt(((double)g) / (t * (t - g))) * cumsums[cord_spec(j, i, p)] -
                                            sqrt(((double)(t - g)) / (t * g)) * (cumsum[j] - cumsums[cord_spec(j, i, p)]);
            }
            else
            {
                CUSUM[cord_spec(j, i, p)] = (cumsum[j] - cumsums[cord_spec(j, i, p)]) / sqrt(((double)(g)));
            }
        }
    }

    double a = 0;
    double nu_a = 0;
    double r = 0;
    double prev_nu_a = 1;

    for (int h = 0; h < len_as; ++h)
    {
        a = as[h];
        nu_a = nu_as[h];
        internal_threshold_matrix(CUSUM, p, len, a, nu_a, (int)(h > 0), prev_nu_a);
        prev_nu_a = nu_a;

        internal_colSum(CUSUM, p, len, &(teststat[len * h]));
    }

    UNPROTECT(12);
    return teststatSEXP;
}

void internal_threshold_matrix(double *matrix, int r1, int c1, double a, double nu_a, int previously_tresholded,
                               double prev_nu_a)
{

    // should be possible to optimize this code a bit
    double a_sq = a * a;
    double true_val = 0;
    if (previously_tresholded)
    {
        for (int j = 0; j < c1; ++j)
        {
            for (int i = 0; i < r1; ++i)
            {
                if (fabs(matrix[cord_spec(i, j, r1)]) > 1e-10)
                {
                    true_val = matrix[cord_spec(i, j, r1)] + prev_nu_a;

                    if (true_val > a_sq)
                    {
                        matrix[cord_spec(i, j, r1)] = true_val - nu_a;
                    }
                    else
                    {
                        matrix[cord_spec(i, j, r1)] = 0.0;
                    }
                }
            }
        }
    }
    else
    {
        for (int j = 0; j < c1; ++j)
        {
            for (int i = 0; i < r1; ++i)
            {
                if (fabs(matrix[cord_spec(i, j, r1)]) > a)
                {
                    matrix[cord_spec(i, j, r1)] = matrix[cord_spec(i, j, r1)] * matrix[cord_spec(i, j, r1)] - nu_a;
                }
                else
                {
                    matrix[cord_spec(i, j, r1)] = 0.0;
                }
            }
        }
    }
}

void internal_colSum(double *matrix, int r1, int c1, double *vector)
{
    // memset(vector, 0, c1);
    for (int j = 0; j < c1; ++j)
    {
        vector[j] = 0.0;
        for (int i = 0; i < r1; ++i)
        {
            vector[j] += matrix[cord_spec(i, j, r1)];
        }
    }
}
