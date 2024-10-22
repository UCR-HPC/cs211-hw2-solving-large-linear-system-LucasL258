#ifndef __MY_C__
#define __MY_C__

#include "include.h"

int mydgetrf(double *A,int *ipiv,int n)
{
    //TODO
    //The return value (an integer) can be 0 or 1
    //If 0, the matrix is irreducible and the result will be ignored
    //If 1, the result is valid
    int i, j, k;

    for (i = 0; i < n-1; i++) {
        // pivoting
        int maxind = i;
        double temp_max = fabs(A[i * n + i]);
        for (int t = i + 1; t < n; t++) {
            if (fabs(A[t * n + i]) > temp_max) {
                maxind = t;
                temp_max = fabs(A[t * n + i]);
            }
        } 

        if (temp_max == 0) {
            return 0;
        }
        else {
            if (maxind != i) {
                // save pivoting information
                int temps = ipiv[i];
                ipiv[i] = ipiv[maxind];
                ipiv[maxind] = temps;
                // swap rows
                double tempv[n];
                for (int x = 0; x < n; x++) {
                    tempv[x] = A[i * n + x];
                    A[i * n + x] = A[maxind * n + x];
                    A[maxind * n + x] = tempv[x];
                }

            }
        }

        // factorization
        for (j = i + 1; j < n; j++) {
            A[j * n + i] = A[j * n + i] / A[i * n + i];
        }

        for (j = i + 1; j < n; j++) {
            for (k = i + 1; k < n; k++) {
                A[j * n + k] = A[j * n + k] - A[j * n + i] * A[i * n + k];;
            }
        }
    }
    return 1;
}

void mydtrsv(char UPLO, double *A,double *B,int n,int *ipiv)
{
    double y[n], x[n]; 

    if (UPLO == 'L') {
        y[0] = B[ipiv[0]];
        for (int i = 1; i < n; i++) {
            y[i] = B[ipiv[i]];
            for (int j = 0; j < i; j++) {
                y[i] -= A[i * n + j] * y[j];
            }
        }

        for (int i = 0; i < n; i++) {
            B[i] = y[i];
        }
    } 

    else if (UPLO == 'U') {
        for (int i = n - 1; i >= 0; i--) {
            for (int j = i + 1; j < n; j++) {
                B[i] -= A[i * n + j] * B[j];
            }
            
            B[i] /= A[i * n + i];
        }
    } 
}

void my_f(double *A,double *B,int n)
{
    int *ipiv=(int*)malloc(n*sizeof(int));
    for (int i=0;i<n;i++)
        ipiv[i]=i;
    if (mydgetrf(A,ipiv,n)==0) 
    {
        printf("LU factoration failed: coefficient matrix is singular.\n");
        return;
    }
    mydtrsv('L',A,B,n,ipiv);
    mydtrsv('U',A,B,n,ipiv);
}

#endif