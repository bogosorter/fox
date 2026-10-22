#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define MIN(A, B) (A < B) ? A : B
#define INF 1e9


void fox(int n, int *matrix);
void fox_step(int n, int *a, int *b, int *c);

void read_matrix(int n, int *matrix);
void print_matrix(int n, int *matrix);


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int n;
    scanf("%d", &n);

    int *full_matrix = malloc(sizeof(int) * n * n);
    read_matrix(n, full_matrix);

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    fox(n, full_matrix);

    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();
    
    printf("%f\n", end_time - start_time);
    print_matrix(n, full_matrix);
    
    MPI_Finalize();
    return 0;
}


void fox(int n, int *matrix) {
    // See report for a detailed description of Fox's algorithm
    size_t size = n * n * sizeof(int);

    int *result = malloc(sizeof(int) * n * n);    
    for (int i = 0; i < n * n; i++) result[i] = INF;
    for (int matrix_power = 1; matrix_power < n; matrix_power *= 2) {
        fox_step(n, matrix, matrix, result);
        memcpy(matrix, result, size);
    }
}

void fox_step(int n, int *a, int *b, int *c) {
    // The normal order of the loops (i, j, k) is altered to (i, k, j) so as to
    // reduce the number of row-major changes, since j is only used to address
    // the column
    int i, j, k;
    for (i = 0; i < n; i++) {
        for (k = 0; k < n; k++) {
            for (j = 0; j < n; j++) {
                c[i * n + j] = MIN(c[i * n + j], a[i * n + k] + b[k * n + j]);
            }
        }
    }
}

void read_matrix(int n, int *matrix) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // The input format uses 0 when refering to the distance to a
            // node to itself (which is actually 0) and when referring to
            // a non-existing path (where we consider the distance to be
            // infinite)
            scanf("%d", &matrix[i * n + j]);
            if (matrix[i * n + j] == 0 && i != j) matrix[i * n + j] = INF;
        }
    }
}

void print_matrix(int n, int *matrix) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j) printf(" ");

            // The input format uses 0 when refering to the distance to a
            // node to itself (which is actually 0) and when referring to
            // a non-existing path (where we consider the distance to be
            // infinite)
            int v = matrix[i * n + j];
            if (v == INF) v = 0;
            printf ("%d", v);
        }
        printf ("\n");
    }
}