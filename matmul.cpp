#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "perf_utilities.h"

#define TEST_FUNC(name) void name( void )
typedef void (* testFunc) ( void );


// Note: All matrices are assumed column-major.


#define M 100
#define N 50
#define P 80

// double x[N];
double c[M*N];
double cc[M*N] = { 0 };
double c0[M*N];
double c1[M*N];
double c2[M*N];
double c3[M*N];

double a[M*P];
double b[P*N];


// Matrix utilities
void mat_init( pcg32_random_t *pcg, double *x, int m, int n )
{
    int dim = m*n;
    for ( int i = 0; i < dim; i++ ) {
        x[i] = (double) pcg32_random_r_rangef( pcg, -100.0f, 100.0f );
    }
}


void mat_print( double *x, int m, int n, const char *name )
{
    printf("Mat (%s) = {\n", name);
    printf("\trows = %u\n", m);
    printf("\tcols = %u\n", n);
    printf("\tdata = {\n\n");
    unsigned int dim;
    unsigned int prec = 1 + (unsigned int) log10( (double) (m * n) );
    for (unsigned int i=0; i<m; ++i) {
        for (unsigned int j=0; j<n; ++j) {
            dim = i + j*m;
            printf("[%.*d] ", prec, dim);
            if ( x[dim] == 0.0  ) { printf("     .0   "); } else { printf("%10.4f",  x[dim]); }
            printf("  ");
        }
        printf("\n");
    }
    printf("\n\t}\n");
    printf("}\n");
}


int mat_equal( double *a, double *b, int m, int n )
{
    for ( int i = 0; i < m*n; i++ ) {
        if ( a[i] != b[i] ) {
            return 0;
        }
    }
    return 1;
}



// Matrix multiplication: c += a*b.
// Dimensions: c is (mxn), a is (mxp), b is (pxn).

void matmul_0( double *c, double *a, double *b, int m, int n, int p )
{
    for ( int i = 0; i < m; i++ ) {
        for ( int j = 0; j < n; j++ ) {
            for ( int k = 0; k < p; k++ ) {
                c[i + j*m] += a[i + k*m] * b[k + j*p];
            }
        }
    }
}

void matmul_1( double *c, double *a, double *b, int m, int n, int p )
{
    for ( int j = 0; j < n; j++ ) {
        for ( int i = 0; i < m; i++ ) {
            for ( int k = 0; k < p; k++ ) {
                c[i + j*m] += a[i + k*m] * b[k + j*p];
            }
        }
    }
}


void matmul_2( double *c, double *a, double *b, int m, int n, int p )
{
    for ( int j = 0; j < n; j++ ) {
        for ( int k = 0; k < p; k++ ) {
            for ( int i = 0; i < m; i++ ) {
                c[i + j*m] += a[i + k*m] * b[k + j*p];
            }
        }
    }
}



void matmul_3( double *c, double *a, double *b, int m, int n, int p )
{

}


// reference test function for matrix initialization and memset
TEST_FUNC(test_reference)
{
    pcg32_random_t pcg;
    pcg32_srandom_r( &pcg, 1234, 5678 );

    memset(c, 0, sizeof(double)*M*N);
    mat_init( &pcg, a, M, P );
    mat_init( &pcg, b, P, N );
}


TEST_FUNC(test_matmul_0)
{
    pcg32_random_t pcg;
    pcg32_srandom_r( &pcg, 1234, 5678 );

    mat_init( &pcg, a, M, P );
    mat_init( &pcg, b, P, N );

    memset(c0, 0, sizeof(double)*M*N);
    matmul_0(c0, a, b, M, N, P);
}


TEST_FUNC(test_matmul_1)
{
    pcg32_random_t pcg;
    pcg32_srandom_r( &pcg, 1234, 5678 );

    mat_init( &pcg, a, M, P );
    mat_init( &pcg, b, P, N );

    memset(c1, 0, sizeof(double)*M*N);
    matmul_1(c1, a, b, M, N, P);
}


TEST_FUNC(test_matmul_2)
{
    pcg32_random_t pcg;
    pcg32_srandom_r( &pcg, 1234, 5678 );

    mat_init( &pcg, a, M, P );
    mat_init( &pcg, b, P, N );

    memset(c2, 0, sizeof(double)*M*N);
    matmul_2(c2, a, b, M, N, P);
}


int main( int argn, const char *argv[] )
{
    Test tests[] = {
        (Test){ test_reference, "Reference Test" },
        (Test){ test_matmul_0,  "matmul_0" },
        (Test){ test_matmul_1,  "matmul_1" },
        (Test){ test_matmul_2,  "matmul_2" },
    };

    RunTests( tests );


    if ( !mat_equal(c, c, M, N) ) {
        printf("ERROR: matrices c and cc not equal\n");
    }

    if ( !mat_equal(c0, c1, M, N) ) {
        printf("ERROR: matrices c0 and c1 not equal\n");
    }

    if ( !mat_equal(c1, c2, M, N) ) {
        printf("ERROR: matrices c1 and c2 not equal\n");
    }

    // mat_print(c0, M, N, "c0");
    // mat_print(c1, M, N, "c1");
    // mat_print(c2, M, N, "c2");

    return 0;
}
