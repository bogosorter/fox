#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

// We define infinity to be a distance large enough to not occur in the data,
// but still allowing infinity + infinity without overflow
#define INF 1e9


// Information about the problem's topology and distribution accross multiple
// processors
typedef struct {
    int processes;
    int rank;
    int n;
    int q;
    int nq;
    int row;
    int column;
    MPI_Comm grid;
    MPI_Comm row_communicator;
    MPI_Comm column_communicator;
    MPI_Datatype submatrix_type;
} TopologyInfo;

// Helper functions
void initialize_mpi(TopologyInfo *topology, int *argc, char ***argv);
void read_size(TopologyInfo *topology);
void create_grid(TopologyInfo *topology);
int *read_input(TopologyInfo *topology);
void write_output(TopologyInfo *topology, int *process_submatrix);
void distributed_matrix_squaring(int *process_submatrix, TopologyInfo *topology);
void matrix_multiplication_step(int n, int *a, int *b, int *result);
int min(int a, int b);


int main(int argc, char **argv) {
    TopologyInfo topology;

    initialize_mpi(&topology, &argc, &argv);
    read_size(&topology);
    create_grid(&topology);
    int *process_matrix = read_input(&topology);

    
    // Use the repeated squaring method to calculate the distance matrix

    for (int step = 1; step < topology.n; step *= 2) {
        distributed_matrix_squaring(process_matrix, &topology);
    }


    write_output(&topology, process_matrix);
    
    free(process_matrix);
    MPI_Finalize();
}

void initialize_mpi(TopologyInfo *topology, int *argc, char ***argv) {
    MPI_Init(argc, argv);
    MPI_Comm_size(MPI_COMM_WORLD, &topology->processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &topology->rank);
}

void read_size(TopologyInfo *topology) {
    if (topology->rank == 0) scanf("%d", & topology->n);
    MPI_Bcast(& topology->n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    topology->q = sqrt(topology->processes);
    if (topology->q * topology->q != topology->processes) {
        MPI_Finalize();
        if (topology->rank == 0) printf("the number of processes must be a perfect square\n");
        exit(0);
    }

    if (topology->n % topology->q != 0) {
        MPI_Finalize();
        if (topology->rank == 0) printf("sqrt(processes) must divide the size of the matrix\n");
        exit(0);
    }
    topology->nq = topology->n / topology->q;
}

void create_grid(TopologyInfo *topology) {
    int dimension_count = 2;
    int dimensions[2] = {topology->q, topology->q};
    int wrap_around[2] = {0, 1};
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, dimension_count, dimensions, wrap_around, reorder, &topology->grid);

    // Create row and column sub-communicators
    
    int free_axis[2] = {0, 1};
    MPI_Cart_sub(topology->grid, free_axis, &topology->row_communicator);
    free_axis[0] = 1;
    free_axis[1] = 0;
    MPI_Cart_sub(topology->grid, free_axis, &topology->column_communicator);

    MPI_Comm_rank(topology->grid, &topology->rank);
    topology->row = topology->rank / topology->q;
    topology->column = topology->rank % topology->q;

    // Create a datatype for identifying a submatrix nq * nq within a larger
    // topology->n * topology->n matrix
    MPI_Type_vector(topology->nq, topology->nq, topology->n, MPI_INT, &topology->submatrix_type);
    MPI_Type_commit(&topology->submatrix_type);
}

int *read_input(TopologyInfo *topology) {
    int *process_matrix = malloc(sizeof(int) * topology->nq * topology->nq);

    if (topology->rank == 0) {
        int *full_matrix = malloc(sizeof(int) * topology->n * topology->n);

        for (int i = 0; i < topology->n; i++) {
            for (int j = 0; j < topology->n; j++) {
                scanf("%d", &full_matrix[i * topology->n + j]);

                // When no edge exists, we assume an infinite weight, except if we
                // are considering the path from a vertex to itself, whose length
                // is always zero
                if (i != j && full_matrix[i * topology->n + j] == 0) full_matrix[i * topology->n + j] = INF;
            }
        }

        // Save process zero's own submatrix
        for (int i = 0; i < topology->nq; i++) {
            for (int j = 0; j < topology->nq; j++) {
                process_matrix[i * topology->nq + j] = full_matrix[i * topology->n + j];
            }
        }

        for (int i = 0; i < topology->q; i++) {
            for (int j = 0; j < topology->q; j++) {
                // Send the submatrix to the corresponding process
                if (i == 0 && j == 0) continue;
                else MPI_Send(&full_matrix[i * topology->nq * topology->n + j * topology->nq], 1, topology->submatrix_type, i * topology->q + j, 0, topology->grid);
            }
        }

        free(full_matrix);
    } else {
        MPI_Recv(process_matrix, topology->nq * topology->nq, MPI_INT, 0, 0, topology->grid, NULL);
    }

    return process_matrix;
}

void write_output(TopologyInfo *topology, int *process_submatrix) {
    if (topology->rank == 0) {
        int *full_matrix = malloc(sizeof(int) * topology->n * topology->n);

        // Copy process zero's own submatrix to the full matrix
        for (int i = 0; i < topology->nq; i++) {
            for (int j = 0; j < topology->nq; j++) {
                full_matrix[i * topology->n + j] = process_submatrix[i * topology->nq + j];
            }
        }

        for (int i = 0; i < topology->q; i++) {
            for (int j = 0; j < topology->q; j++) {
                // Send the submatrix to the corresponding process
                if (i == 0 && j == 0) continue;
                else MPI_Recv(&full_matrix[i * topology->nq * topology->n + j * topology->nq], 1, topology->submatrix_type, i * topology->q + j, 0, topology->grid, NULL);
            }
        }

        for (int i = 0; i < topology->n; i++) {
            for (int j = 0; j < topology->n; j++) {
                if (j) printf(" ");

                // When no path exists, we output 0
                if (full_matrix[i * topology->n + j] == INF) printf("%d", 0);
                else printf("%d", full_matrix[i * topology->n + j]);
            }
            printf("\n");
        }
    } else {
        MPI_Send(process_submatrix, topology->nq * topology->nq, MPI_INT, 0, 0, topology->grid);
    }
}

void distributed_matrix_squaring(int *process_submatrix, TopologyInfo *topology) {
    int submatrix_size = topology->nq * topology->nq;
    
    int *result = malloc(sizeof(int) * submatrix_size);
    for (int i = 0; i < topology->nq; i++) {
        for (int j = 0; j < topology->nq; j++) {
            result[i * topology->nq + j] = INF;
        }
    }

    int *a = malloc(sizeof(int) * topology->nq * topology->nq);
    memcpy(a, process_submatrix, sizeof(int) * submatrix_size);
    
    int *b = malloc(sizeof(int) * submatrix_size);

    for (int step = 0; step < topology->q; step++) {
        // At the step-th iteration, the value b on the (step + row)-th column is broadcast accross the row
        int source = (step + topology->row) % topology->q;
        if (source == topology->column) memcpy(b, process_submatrix, sizeof(int) * submatrix_size);
        MPI_Bcast(b, submatrix_size, MPI_INT, source, topology->row_communicator);

        matrix_multiplication_step(topology->nq, b, a, result);

        // After having started with its own value, each process broadcasts the value of b accross the column
        if (step != topology->q - 1) {
            int destination = (topology->row + topology->q - 1) % topology->q; // we add topology->q to avoid negative values on the remainder
            int source = (topology->row + 1) % topology->q;
            MPI_Sendrecv_replace(a, submatrix_size, MPI_INT, destination, 0, source, 0, topology->column_communicator, NULL);
        }
    }

    memcpy(process_submatrix, result, sizeof(int) * submatrix_size);
}

void matrix_multiplication_step(int n, int *a, int *b, int *result) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                result[i * n + j] = min(result[i * n + j], a[i * n + k] + b[k * n + j]);
            }
        }
    }
}

int min(int a, int b) {
    if (a < b) return a;
    return b;
}