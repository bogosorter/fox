#include <stdio.h>

void square_matrix(int n, int matrix[n][n], int result[n][n]);

int main() {
    int n;
    scanf("%d", &n);

    int matrix[n][n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            scanf("%d", &matrix[i][j]);
        }
    }

    int result[n][n];
    square_matrix(n, matrix, result);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j) printf(" ");
            printf("%d", result[i][j]);
        }
        printf("\n");
    }
}

void square_matrix(int n, int matrix[n][n], int result[n][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = 0;
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                result[i][j] += matrix[i][k] * matrix[k][j];
            }
        }
    }
}