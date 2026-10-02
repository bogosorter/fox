#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define INFINITY 1e9

void swap(int **a, int **b);
int min(int a, int b);
void square_matrix(int n, int *matrix, int *result);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int processes, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    // Read matrix size and check that the number of processes is as expected

    int n;
    if (rank == 0) scanf("%d", &n);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (processes != n * n) {
        MPI_Finalize();
        if (rank == 0) printf("for now, there must be as many processes as entries in the matrix\n");
        exit(0);
    }

    // Create the square network topology

    MPI_Comm grid;
    int dimension_count = 2;
    int dimensions[2] = {n, n};
    int wrap_around[2] = {0, 1};
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, dimension_count, dimensions, wrap_around, reorder, &grid);

    MPI_Comm_rank(grid, &rank);
    int row = rank / n;
    int column = rank % n;


    // Process zero reads the input and sends to the corresponding processes

    int process_value;

    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int v;
                scanf("%d", &v);

                // When no edge exists, we assume an infinite weight, except if we
                // are considering the path from a vertex to itself, whose length
                // is always zero
                if (i != j && v == 0) v = INFINITY;

                // Send to the corresponding process, or keep the value if it
                // belongs to process 1.
                if (i == 0 && j == 0) process_value = v;
                else MPI_Send(&v, 1, MPI_INT, i * n + j, 0, grid);
            }
        }
    } else {
        MPI_Recv(&process_value, 1, MPI_INT, 0, 0, grid, NULL);
    }


    // Perform a single vertical shift of the numbers

    int destination_row = (row + n - 1) % n; // we add n to avoid negative values on the remainder
    int source_row = (row + 1) % n;
    MPI_Send(&process_value, 1, MPI_INT, destination_row * n + column, 0, grid);
    MPI_Recv(&process_value, 1, MPI_INT, source_row * n + column, 0, grid, NULL);


    // Gather the result from all processes and print it to stdout

    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j) printf(" ");

                int v;
                if (i == 0 && j == 0) v = process_value;
                else MPI_Recv(&v, 1, MPI_INT, i * n + j, 0, grid, NULL);

                // When no path exists, we output 0
                if (v == INFINITY) printf("%d", 0);
                else printf("%d", v);
            }
            printf("\n");
        }
    } else {
        MPI_Send(&process_value, 1, MPI_INT, 0, 0, grid);
    }

    MPI_Finalize();
}

void square_matrix(int n, int *matrix, int *result) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i * n + j] = matrix[i * n + j];
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                result[i * n + j] = min(result[i * n + j], matrix[i * n + k] + matrix[k * n + j]);
            }
        }
    }
}

void swap(int **a, int **b) {
    int *temporary = *b;
    *b = *a;
    *a = temporary;
}

int min(int a, int b) {
    if (a < b) return a;
    return b;
}