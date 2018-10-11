#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "mkl_lapacke.h"

//Function to get a column of a matrix
//Takes a 2d array (matrix), the dimension (size), and the number of the column to get (col)
//Return an array which contain matrix[:,col]
double *get_column_matrix(double **matrix, int size, int col){
    double *vector = malloc( size * sizeof(double));
    for (int i = 0; i < size; ++i)
    {
        vector[i] = matrix[i][col];
    }
    return vector;
}

//Funchion which computes the scalar product of two vectors
//Takes two arrays (v1 and v2) and the size of the arrays (size)
//Returns a double which corresponds to the scalar product of v1 and v2
double dot_product_vectors(double *v1, double *v2, int size){
    double res = 0;
    for (int i = 0; i < size; ++i)
    {
        res += v1[i]*v2[i];
    }
    return res;
}

double *multiply_vector_scalar(double *v, double scalar, int size){
    for (int i = 0; i < size; ++i)
    {
        v[i] *= scalar;
    }
    return v;
}

//Function which divide a vector by its norm
//Takes an array v, and the size of the vector (size)
//Returns an anrray which containes each element of v divided by norm(v)
double *norm_vector(double *v, int size){
    double norm = sqrt(dot_product_vectors(v,v,size));
    for (int i = 0; i < size; ++i)
    {
        v[i]/=norm;
    }
    return v;
}


double **create_matrix(double *matrix, int size){
    double **new_matrix = NULL;
    new_matrix = malloc( size * sizeof(double *));
    for(int i = 0; i < size; i++)
    {
        new_matrix[i] = malloc( size * sizeof(double));
    }
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            new_matrix[i][j] = matrix[i*size+j];
        }
    }
    return new_matrix;
}

double **multiply_2matrix(double **A, double **B, int size){
    double **res = NULL;
    res = malloc( size * sizeof(double *));
    for(int i = 0; i < size; i++)
    {
        res[i] = malloc( size * sizeof(double));
    }
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){ 
            res[i][j]=0; 
            for(int k = 0; k < size; ++k){
                res[i][j] = res[i][j] + A[i][k] * B[k][j]; 
            } 
        }
    } 
    return res;
}

double **transpose_matrix(double **mat, int size){
    double **res = NULL;
    res = malloc( size * sizeof(double *));
    for(int i = 0; i < size; i++)
    {
        res[i] = malloc( size * sizeof(double));
    }
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            res[i][j] = mat[j][i];
        }
    }
    return res;
}


void print_2Dmatrix(const char *name, double **matrix, int size)
{
    int i, j;
    printf("matrix: %c \n", matrix);

    for (i = 0; i < size; i++)
    {
            for (j = 0; j < size; j++)
            {
                printf("%f ", matrix[i][j]);
            }
            printf("\n");
    }
}


double ** QR_GramSchmidt(double **A, double **B, int size){
    //Initialisation of Q, R, X and QTB = Q^T*B matrix
    double **Q = NULL;
    Q = malloc( size * sizeof(double *));
    double **R = NULL;
    R = malloc( size * sizeof(double *));
    double **QTB = NULL;
    QTB = malloc( size * sizeof(double *));
    double **X = NULL;
    X = malloc( size * sizeof(double *));
    for(int i = 0; i < size; i++)
    {
        Q[i] = malloc( size * sizeof(double));
        R[i] = malloc( size * sizeof(double));
        QTB[i] = malloc( size * sizeof(double));
        X[i] = malloc( size * sizeof(double));
    }
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            Q[i][j] = 0;
            R[i][j] = 0;
            X[i][j] = 0;
        }
    }

    //Calculation of Q and R matrix with Gram Schmidt method
    double *v = NULL;
    v = malloc(sizeof(double)*size);
    v = get_column_matrix(A,size,0);
    R[0][0] = sqrt(dot_product_vectors(v,v,size));
    double *temp_q = (double *)malloc(sizeof(double)*size);
    temp_q = norm_vector(v,size);
    for (int i = 0; i < size; ++i)
    {
        Q[i][0] = temp_q[i];
    }

    double norm;
    double *q = NULL;
    q = (double *)malloc(sizeof(double)*size); 
    for (int i = 0; i < size-1; ++i)
    {
        v = get_column_matrix(A,size,i+1);
        for (int k = 0; k < i+1; ++k)
        {
            q = get_column_matrix(Q, size, k);
            R[k][i+1] = dot_product_vectors(v, q, size);
            double *v_temp = multiply_vector_scalar(q, R[k][i+1], size);
            for (int j = 0; j < size; ++j)
            {
                v[j] -= v_temp[j];
            }
        }
        norm = sqrt(dot_product_vectors(v,v,size));
        R[i+1][i+1] = norm;
        temp_q = norm_vector(v,size);
        for (int j = 0; j < size; ++j)
        {
            Q[j][i+1] = temp_q[j];
        }
        
    }
    //print_2Dmatrix("Q", Q, size);
    //print_2Dmatrix('R', R, size);

    //Solve of AX = B system :
    Q = transpose_matrix(Q, size);
    QTB = multiply_2matrix(Q, B, size);
    for (int j = size - 1; j >= 0; --j){
        for(int i = size - 1; i >= 0; --i){ 
            double sum = 0;
            for(int k = size - 1; k >= 0; --k){ 
                sum += R[i][k] * X[k][j]; 
            } 
            if (R[i][i] != 0){
                X[i][j] = (QTB[i][j] - sum)/R[i][i];
            } 
        }
     }
     //print_2Dmatrix('X', X, size);
    return X;
}

int check_result_matrix(double **bref, double **b, int size) {
    int i;
    for(i=0;i<size;i++) {
        for (int j = 0; j < size; ++j)
        {
            if (bref[i][j]!=b[i][j]) return 0;
        }
    }
    return 1;
}


double *generate_matrix(int size)
{
    int i;
    double *matrix = (double *)malloc(sizeof(double) * size * size);
    srand(1);

    for (i = 0; i < size * size; i++)
    {
        matrix[i] = rand() % 100;
    }

    return matrix;
}

void print_matrix(const char *name, double *matrix, int size)
{
    int i, j;
    printf("matrix: %s \n", matrix);

    for (i = 0; i < size; i++)
    {
            for (j = 0; j < size; j++)
            {
                printf("%f ", matrix[i * size + j]);
            }
            printf("\n");
    }
}

int check_result(double *bref, double *b, int size) {
    int i;
    for(i=0;i<size*size;i++) {
        if (bref[i]!=b[i]) return 0;
    }
    return 1;
}

int my_dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb) {

    //Replace with your implementation
    LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipi
v, b, ldb);
    
}


    void main(int argc, char *argv[])
    {

        int size = atoi(argv[1]);

        double *a, *aref;
        double *b, *bref;

        a = generate_matrix(size);
        aref = generate_matrix(size);        
        b = generate_matrix(size);
        bref = generate_matrix(size);
        double **a_new = create_matrix(a, size);
        double **b_new = create_matrix(b, size);
        clock_t tStart = clock();
        double **X = QR_GramSchmidt(a_new, b_new, size);

        //print_matrix("A", a, size);
        //print_matrix("B", b, size);

        // Using MKL to solve the system
        //MKL_INT n = size, nrhs = size, lda = size, ldb = size, info;
        //MKL_INT *ipiv = (MKL_INT *)malloc(sizeof(MKL_INT)*size);

   
        //info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
        //printf("Time taken by MKL: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

            
        //MKL_INT *ipiv2 = (MKL_INT *)malloc(sizeof(MKL_INT)*size);        
        //my_dgesv(n, nrhs, a, lda, ipiv2, b, ldb);
        printf("Time taken by my implementation: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
        //bref = create_matrix(bref, size);
        /*if (check_result_matrix(bref,X,size)==1)
            printf("Result is ok!\n");
        else    
            printf("Result is wrong!\n");
        */
        //print_matrix("X", b, size);
        //print_matrix("Xref", bref, size);
    }
