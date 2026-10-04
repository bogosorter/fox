#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define INF 1e9

typedef struct {
    int q;
    int row;
    int column;
    MPI_Comm grid;
    MPI_Comm row_communicator;
    MPI_Comm column_communicator;
    MPI_Datatype submatrix_type;
} TopologyInfo;

void swap(int **a, int **b);
int min(int a, int b);
void square_matrix(int n, int *matrix, int *result);
int distributed_matrix_squaring(int process_value, TopologyInfo *topology);

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int processes, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    // Read matrix size and validate the number of processes

    int n;
    if (rank == 0) scanf("%d", &n);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int q = sqrt(processes);
    if (q * q != processes) {
        MPI_Finalize();
        if (rank == 0) printf("the number of processes must be a perfect square\n");
        exit(0);
    }

    if (n % q != 0) {
        MPI_Finalize();
        if (rank == 0) printf("sqrt(processes) must divide the size of the matrix\n");
        exit(0);
    }
    int nq = n / q;


    // Create the grid network topology

    TopologyInfo topology;
    int dimension_count = 2;
    int dimensions[2] = {q, q};
    int wrap_around[2] = {0, 1};
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, dimension_count, dimensions, wrap_around, reorder, &topology.grid);

    // Create row and column sub-communicators
    int free_axis[2] = {0, 1};
    MPI_Cart_sub(topology.grid, free_axis, &topology.row_communicator);
    free_axis[0] = 1;
    free_axis[1] = 0;
    MPI_Cart_sub(topology.grid, free_axis, &topology.column_communicator);

    MPI_Comm_rank(topology.grid, &rank);
    topology.q = q;
    topology.row = rank / q;
    topology.column = rank % q;

    // Create a datatype for identifying a submatrix nq * nq within a larger
    // n * n matrix
    MPI_Type_vector(nq, nq, n, MPI_INT, &topology.submatrix_type);
    MPI_Type_commit(&topology.submatrix_type);


    // Process zero reads the input and sends to the corresponding processes

    int *process_matrix = malloc(sizeof(int) * nq * nq);

    if (rank == 0) {
        int *full_matrix = malloc(sizeof(int) * n * n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                scanf("%d", &full_matrix[i * n + j]);

                // When no edge exists, we assume an infinite weight, except if we
                // are considering the path from a vertex to itself, whose length
                // is always zero
                if (i != j && full_matrix[i * n + j] == 0) full_matrix[i * n + j] = INF;
            }
        }

        // Save process zero's own submatrix
        for (int i = 0; i < nq; i++) {
            for (int j = 0; j < nq; j++) {
                process_matrix[i * nq + j] = full_matrix[i * n + j];
            }
        }

        for (int i = 0; i < q; i++) {
            for (int j = 0; j < q; j++) {
                // Send the submatrix to the corresponding process
                if (i == 0 && j == 0) continue;
                else MPI_Send(&full_matrix[i * nq * n + j * nq], 1, topology.submatrix_type, i * q + j, 0, topology.grid);
            }
        }

        free(full_matrix);
    } else {
        MPI_Recv(process_matrix, nq * nq, MPI_INT, 0, 0, topology.grid, NULL);
    }

    
    // Use the repeated squaring method to calculate the distance matrix

    // for (int step = 1; step < topology.n; step *= 2) {
    //     process_value = distributed_matrix_squaring(process_value, &topology);
    // }


    // Gather the result from all processes and print it to stdout
    
    if (rank == 0) {
        int *full_matrix = malloc(sizeof(int) * n * n);

        // Copy process zero's own submatrix to the full matrix
        for (int i = 0; i < nq; i++) {
            for (int j = 0; j < nq; j++) {
                full_matrix[i * n + j] = process_matrix[i * nq + j];
            }
        }

        for (int i = 0; i < q; i++) {
            for (int j = 0; j < q; j++) {
                // Send the submatrix to the corresponding process
                if (i == 0 && j == 0) continue;
                else MPI_Recv(&full_matrix[i * nq * n + j * nq], 1, topology.submatrix_type, i * q + j, 0, topology.grid, NULL);
            }
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j) printf(" ");

                // When no path exists, we output 0
                if (full_matrix[i * n + j] == INF) printf("%d", 0);
                else printf("%d", full_matrix[i * n + j]);
            }
            printf("\n");
        }
    } else {
        MPI_Send(process_matrix, nq * nq, MPI_INT, 0, 0, topology.grid);
    }
    
    free(process_matrix);
    MPI_Finalize();
}

int distributed_matrix_squaring(int process_value, TopologyInfo *topology) {
    int result = INF;
    int a = process_value;

    for (int i = 0; i < topology->q; i++) {
        // At the i-th step, the value b on the (i + row)-th column is broadcast accross the row
        int b = process_value;
        int source = (i + topology->row) % topology->q;
        MPI_Bcast(&b, 1, MPI_INT, source, topology->row_communicator);

        result = min(result, a + b);

        // After having started with its own value, each process broadcasts its value accross the column
        if (i != topology->q - 1) {
            int destination = (topology->row + topology->q - 1) % topology->q; // we add n to avoid negative values on the remainder
            int source = (topology->row + 1) % topology->q;
            MPI_Send(&a, 1, MPI_INT, destination, 0, topology->column_communicator);
            MPI_Recv(&a, 1, MPI_INT, source, 0, topology->column_communicator, NULL);
        }
    }

    return result;
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