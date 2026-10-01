#include <stdio.h>
#include <stdlib.h>

#define INFINITY 1e9

void swap(int **a, int **b);
int min(int a, int b);
void square_matrix(int n, int *matrix, int *result);

int main() {
    int n;
    scanf("%d", &n);

    // Read input
    int *matrix = malloc(sizeof(int) * n * n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            scanf("%d", &matrix[i * n + j]);

            // When no edge exists, we assume an infinite weight, except if we
            // are considering the path from a vertex to itself, whose length
            // is always zero
            if (i != j && matrix[i * n + j] == 0) matrix[i * n + j] = INFINITY;
        }
    }

    // Apply the repeated squaring method
    int *temporary = malloc(sizeof(int) * n * n);
    for (int i = 1; i < n; i *= 2) {
        square_matrix(n, matrix, temporary);
        swap(&matrix, &temporary);
    }

    // Write resulting matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j) printf(" ");

            // When no path exists, we output 0
            if (matrix[i * n + j] == INFINITY) printf("%d", 0);
            else printf("%d", matrix[i * n + j]);
        }
        printf("\n");
    }

    free(temporary);
    free(matrix);
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