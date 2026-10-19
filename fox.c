#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define MIN(A, B) (A < B) ? A : B
#define ROOT 0
#define MPI_TAG 1
#define INF 1e9


// Contains information about the grid topology
typedef struct {
    int n; // size of the original matrix
    int q; // size of the matrix of processes
    int nq; // size of each submatrix
    int rank;
    int row;
    int col;
    int root; // grid rank of the world root process 
    MPI_Comm communicator;
    MPI_Comm row_communicator;
    MPI_Comm column_communicator;
} GridInfo;

// Contains information about a custom MPI datatype used to transmit
// submatrices
typedef struct {
    MPI_Datatype blocktype; // The actual derived type
    int *counts; // The number of blocks to send to each process (always 1)
    int *displacements; // The offsets of those blocks relative to the full matrix
} BlockLayout;


// Ensures that the number of processes is valid and respects the size of the
// matrix
int validate_grid(int n, int processes, int q, int nq);

// Creates the grid of processes
GridInfo create_grid(int n, int rank, int processes);

// Creates a derived MPI type for a submatrix
BlockLayout create_block_layout(int n, int q, int rank);

// Performs Fox's algorithm
void fox(int *prossess_matrix, GridInfo *grid);

// Performs a single sub-matrix multiplication step
void fox_step(int n, int *a, int *b, int *c);

void read_matrix(int n, int *matrix);
void print_matrix(int n, int *matrix);


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);

    int n;
    if (rank == ROOT) scanf("%d", &n);
    MPI_Bcast(&n, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    int q = sqrt(processes);
    int nq = n / q;
    if (validate_grid(n, processes, q, nq) == 1) {
        if (rank == ROOT) printf("ERROR: Invalid configuration!\n");
        MPI_Finalize();
        exit(1);
    }

    int *full_matrix = malloc(sizeof(int) * n * n);
    if (rank == ROOT) {
        read_matrix(n, full_matrix);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    BlockLayout layout = create_block_layout(n, q, rank);
    GridInfo grid = create_grid(n, rank, processes);
    
    // The root process of scatterv must not be ROOT, as the process with rank
    // ROOT on the grid communicator may no longer be the process that read the
    // input, since reorder is set to 1
    int *sub_matrix = malloc(sizeof(int) * nq * nq);
    MPI_Scatterv(full_matrix, layout.counts, layout.displacements, layout.blocktype, sub_matrix, nq * nq, MPI_INT, grid.root, grid.communicator);

    fox(sub_matrix, &grid);

    // See above comment for a justification on why ROOT is not used
    MPI_Gatherv(sub_matrix, nq * nq, MPI_INT, full_matrix, layout.counts, layout.displacements, layout.blocktype, grid.root, grid.communicator);

    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();
    
    if (rank == ROOT) {
        printf("%f\n", end_time - start_time);
        print_matrix(n, full_matrix);
    }

    free(layout.counts);
    free(layout.displacements);
    MPI_Comm_free(&grid.communicator);
    MPI_Comm_free(&grid.row_communicator);
    MPI_Comm_free(&grid.column_communicator);
    MPI_Type_free(&layout.blocktype);
    MPI_Finalize();

    return 0;
}

GridInfo create_grid(int n, int rank, int processes) {
    GridInfo grid;
    grid.n = n;
    grid.q = sqrt(processes);
    grid.nq = n / grid.q;

    int dimensions[2] = {grid.q, grid.q};
    int periods[2] = {1, 1};
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, periods, reorder, &(grid.communicator));

    int coordinates[2];
    MPI_Comm_rank(grid.communicator, &(grid.rank));
    MPI_Cart_coords(grid.communicator, grid.rank, 2, coordinates);
    grid.row = coordinates[0];
    grid.col = coordinates[1];

    int remain_row[2] = {0, 1};
    int remain_col[2] = {1, 0};
    MPI_Cart_sub(grid.communicator, remain_row, &(grid.row_communicator));
    MPI_Cart_sub(grid.communicator, remain_col, &(grid.column_communicator));

    if (rank == ROOT) grid.root = grid.rank;
    MPI_Bcast(&grid.root, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    return grid;
}

BlockLayout create_block_layout(int n, int q, int rank) {
    BlockLayout layout = {};

    int nq = n / q;

    // Due to the memory layout of the full matrix, the submatrices are not
    // contiguous. As such, sending the second matrix requires starting at a
    // position that is in the middle of the first matrix. MPI only accepts
    // displacements that are multiples of the datasize, though, which requires
    // us to shrink the virtual size of our created type.
    MPI_Datatype block;
    MPI_Type_vector(nq, nq, n, MPI_INT, &block);
    MPI_Type_create_resized(block, 0, sizeof(int), &layout.blocktype);
    MPI_Type_commit(&layout.blocktype);
    MPI_Type_free(&block);

    if (rank == 0) {
        layout.counts = malloc(q * q * sizeof(int));
        layout.displacements = malloc(q * q * sizeof(int));
        for (int i = 0; i < q; i++) {
            for (int j = 0; j < q; j++) {
                layout.counts[i * q + j] = 1;
                layout.displacements[i * q + j] = i * nq * n + j * nq;
            }
        }
    }
    return layout;
}

void fox(int *submatrix, GridInfo *grid) {
    // See report for a detailed description of Fox's algorithm

    int nq = grid->nq;
    size_t size = nq * nq * sizeof(int);

    int *result = malloc(sizeof(int) * nq * nq);
    int *row_matrix = malloc(sizeof(int) * nq * nq);
    int *column_matrix = malloc(sizeof(int) * nq * nq);
    memcpy(result, submatrix, size);
    
    int source = (grid->row + 1) % grid->q;
    int dest = (grid->row + grid->q - 1) % grid->q; // we add grid->q to avoid negative remainder values
    
    for (int matrix_power = 1; matrix_power < grid->n; matrix_power *= 2) {
        memcpy(column_matrix, submatrix, size);
        
        for (int step = 0; step < grid->q; step++) {
            int row_root = (grid->row + step) % grid->q;
            if (row_root == grid->col) memcpy(row_matrix, submatrix, size);
            MPI_Bcast(row_matrix, nq * nq, MPI_INT, row_root, grid->row_communicator);
            
            fox_step(nq, row_matrix, column_matrix, result);
            
            if (step < grid->q - 1) {
                MPI_Sendrecv_replace(column_matrix, nq * nq, MPI_INT, dest, MPI_TAG, source, MPI_TAG, grid->column_communicator, MPI_STATUS_IGNORE);
            }
        }

        memcpy(submatrix, result, size);
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

int validate_grid(int n, int processes, int q, int nq) {
    if (q * q != processes) return 1;
    if (n % q != 0) return 1;
    return 0;
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