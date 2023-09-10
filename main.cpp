#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>

// Funciones para la construccion del Sistema de Ecuaciones 
double **getMesh(double, double, double, double, unsigned int , unsigned int );
double *getBvector(double(*)(double,double), double(*)(double, double, unsigned int, unsigned int, 
                   int, int), double **, unsigned int, unsigned int, double, double, double);
double **getAmatrix(unsigned int, unsigned int, double, double, double);
double **getTridiag(unsigned int, double, double, double);
double **getDiag(unsigned int, double, double);
double f(double, double);
double g(double, double, unsigned int, unsigned int, int, int);

// Lectura y escritura de matrices CSR
void readMatrix_CSR( const char*, unsigned int* , unsigned int* ,
                     unsigned int** , double*** , unsigned int*** );
void writeMatrix_CSR(const char *, double **, unsigned int, unsigned int);

// Impresion de matriz y liberacion de memoria
void printMatrix(double **, unsigned int, unsigned int);
void printVector(double *, unsigned int);
void deleteMatrix(double **, unsigned int);
void deleteIntMatrix(unsigned int **, unsigned int);

// Funciones necesarias para los Métodos directos
void gaussElim(double **, double *, unsigned int);
void factLU(double **, unsigned int);
void cholesky(double **, unsigned int);
double **qr_gramSchmidt(double ***, unsigned int, unsigned int);
double **transposedMatrix(double **, unsigned int, unsigned int);
double scalarProd(double *, double *, unsigned int);
double *vectMatrixProd(double **, double *, unsigned int);
double *solve_ATS(double **, double *, unsigned int );
double *solve_ATI(double **, double *, unsigned int );
double *solve_ATI_LU(double **, double *, unsigned int );

// Funciones necesarias para los Métodos iterativos
double *jacobi(double **, unsigned int **, unsigned int *, double *, double *, double, unsigned int);
double *gauss_seidel(double **, unsigned int **, unsigned int *, double *, double *, double, unsigned int);
double *grad_conj(double **, unsigned int **, unsigned int *, double *, double *, double, unsigned int);
int getColumnJA(unsigned int **, unsigned int *, unsigned int, unsigned int);
double get_aij(double **, unsigned int **, unsigned int *, unsigned int, unsigned int);
double vectNorm(double *, unsigned int);
double *substractVect(double *, double *, unsigned int);
double *startVector(unsigned int);
double *vectMatrixProd_CSR(double **, unsigned int **, unsigned int *, double *, unsigned int);

// Guardar la solucion, necesario para graficar
void writeSolve(const char *, double **, double *, unsigned int, unsigned int);

using namespace std;
using namespace std::chrono; 

int main(int argc, char const *argv[]){
    double **Mesh, **Amatrix, a, b, c, d, h, k, D, *b_vector, *u_Vector;
    double *yVector, **Rmatrix, *qt_matrix, *vector_QTb;  /* necesarios para los metodos directos. */
    unsigned int N, M, n;
    double **AAmatrix, *x0Vector, tol;
    unsigned int **JAmatrix, *n_Zmatrix, n_rows, n_cols; /* necesario para matriz CSR */

    D = 45;
    a = 0;
    b = 1;
    c = 0;
    d = 1;
    h = 0.25;
    k = 0.25;
    N = (b-a)/h + 2;
    M = (d-c)/k + 1;
    n = (N - 3) * (M - 2);  /* tamaño de la matriz : nxn, tamaño del vector: n */

    Mesh = getMesh(a, c, h, k, N, M);
    b_vector = getBvector(f, g, Mesh, N, M, h, k, D);
    Amatrix = getAmatrix(N, M, D, h, k);

    writeMatrix_CSR("A_matrix_CSR.txt", Amatrix, n, n);

    cout << "Discretización del dominio: " << endl;
    printMatrix(Mesh, M, N);

    cout << "\nVector de términos independientes:" << endl;
    printVector(b_vector,n);

    cout << "\nMatriz de coeficientes: " << endl;
    printMatrix(Amatrix, n, n);

    /* ---------------------- MÉTODOS DIRECTOS ----------------------------------*/

    // cout << "\n Solución por eliminación Gaussiana: " << endl;
    // gaussElim(Amatrix, b_vector, n);
    // u_Vector = solve_ATS(Amatrix,b_vector, n);
    // printVector(u_Vector,n);

    // cout << "\nSolución Factorización LU: " << endl;
    // factLU(Amatrix, n);
    // yVector = solve_ATI(Amatrix, b_vector, n);
    // u_Vector = solve_ATS(Amatrix,yVector, n);
    // printVector(u_Vector,n);

    cout << "\n Solución Factorización de Cholesky: " << endl;
    cholesky(Amatrix, n);
    yVector = solve_ATI(Amatrix, b_vector, n);
    u_Vector = solve_ATS(Amatrix,yVector, n);
    printVector(u_Vector,n);

    
    // cout << "\n Solución Factorización QR: " << endl;
    // Rmatrix = qr_gramSchmidt(&Amatrix, n, n);
    // vector_QTb = vectMatrixProd( transposedMatrix(Amatrix, n, n), b_vector, n );
    // u_Vector = solve_ATS(Rmatrix, vector_QTb, n);
    // printVector(u_Vector,n);

    
    /* ------------------------------ MÉTODOS ITERATIVOS ----------------------------------*/

    tol = pow(10,-6);
    readMatrix_CSR("A_matrix_CSR.txt", &n_rows, &n_cols, &n_Zmatrix, &AAmatrix, &JAmatrix);
    x0Vector = startVector(n_rows);

    // cout << "\nSolución Método de Jacobi: " << endl;
    // u_Vector = jacobi(AAmatrix, JAmatrix, n_Zmatrix, b_vector, x0Vector, tol, n_rows);
    // printVector(u_Vector, n_rows);

    // cout << "\nSolución Método de Gauss-Seidel: " << endl;
    // u_Vector = gauss_seidel(AAmatrix, JAmatrix, n_Zmatrix, b_vector, x0Vector, tol, n_rows);
    // printVector(u_Vector, n_rows);

    // cout << "\nSolución Método del Gradiente Conjugado: " << endl;
    // u_Vector = grad_conj(AAmatrix, JAmatrix, n_Zmatrix, b_vector, x0Vector, tol, n_rows);
    // printVector(u_Vector, n_rows);

    
    // auto start = high_resolution_clock::now();           //Para medir el tiempo de ejecucion de un 
    // auto stop = high_resolution_clock::now();            //determinado bloque de codigo 
    // auto duration = duration_cast<nanoseconds>(stop - start);
    // cout << "\nTiempo de ejecución: " << duration.count() * 1.0000e-9 << " segundos" << endl;

    
    writeSolve("Solution.csv", Mesh, u_Vector, N, M);


    deleteMatrix(Mesh, M);
    deleteMatrix(Amatrix, n);
    // deleteMatrix(Rmatrix, n);
    deleteMatrix(AAmatrix, n_rows);
    deleteIntMatrix(JAmatrix, n_rows);
    delete [] n_Zmatrix;
    delete [] x0Vector;
    // delete [] vector_QTb;
    delete [] b_vector;
    delete [] u_Vector;
    delete [] yVector;
    exit(EXIT_SUCCESS);
}


/* ---------------------------------Generando la malla--------------------------------------- */

double **getMesh(double a, double c, double h, double k, unsigned int N, unsigned int M){
    /* x en [a,b], y en [c,d], h, k son los tamaños de paso para generar los xi, yj, 
    N y M son el numero de nodos */
    double **xyMesh = new double *[M];    /* Creación de la malla (arreglo bidimensional) */
    for (int i = 0; i < M; i++){
        xyMesh[i] = new double [N];
        xyMesh[i][0] = c + i*k;           /* Almacenando los yj en la 1er posicion de cada fila */
        for (int j = 1; j < N; j++){
            xyMesh[i][j] = a + (j-1)*h;
        } 
    }
    return xyMesh;
}

/* ---------------------------- Matriz A de coeficientes ---------------------- */

/*Bloque tridiagonal*/
double **getTridiag(unsigned int N, double D, double h, double k){
    double **tridiag_matrix = new double *[N-3];
    for (int i = 0; i < N-3; i++){
        tridiag_matrix[i] = new double [N-3];
        for (int j = 0; j < N-3; j++){
            if (i == j){
                tridiag_matrix[i][j] = 2*D*(1/(h*h)+1/(k*k));
            } else if ( j == i + 1 || j == i - 1 ){
                tridiag_matrix[i][j] = -D/(h*h);
            } else {
                tridiag_matrix[i][j] = 0.;
            }
        }
    }
    return tridiag_matrix;
}

/*Bloque diagonal*/
double **getDiag(unsigned int N, double D, double k){
    double **diag_matrix = new double *[N-3];
    for (int i = 0; i < N-3; i++){
        diag_matrix[i] = new double [N-3];
        for (int j = 0; j < N-3; j++){
            if (j == i){
                diag_matrix[i][j] = -D/(k*k);
            } else {
                diag_matrix[i][j] = 0.;
            }
        }
    }
    return diag_matrix;
}

/* Función para combinar bloques tridiagonales y diagonales creando así la matriz 
   de coeficientes deseada.*/

double **getAmatrix(unsigned int N, unsigned int M, double D, double h, double k) {
    int n = (N - 3) * (M - 2);
    double **Amatrix = new double *[n];
    double **tridiag = getTridiag(N, D, h, k);
    double **diag = getDiag(N, D, k);
    
    for (int i = 0; i < n; i++){
        Amatrix[i] = new double [n];
        for (int j = 0; j < n; j++){
            if ( i / (N-3) == j / (N-3) ){                   /* condicion para los bloques tridiag */
                Amatrix[i][j] = tridiag[i % (N-3)][j % (N-3)];
            } else if (i / (N-3) == j / (N-3) + 1 ){         /* rellena los bloques diagonal inferior */
                Amatrix[i][j] = diag[i % (N-3)][j % (N-3)];
            } else if (i / (N-3) == j / (N-3) - 1 ){         /* rellena los bloques diagonal superior */
                Amatrix[i][j] = diag[i % (N-3)][j % (N-3)];
            } else {
                Amatrix[i][j] = 0.;
            }
        }
    }
    deleteMatrix(tridiag, N-3);
    deleteMatrix(diag, N-3);

    return Amatrix;
}

/* ---------------------------- Vector b de terminos independientes ---------------------- */

double *getBvector(double(*f_function)(double,double), 
                   double(*g_function)(double, double, unsigned int, unsigned int, int, int), 
                   double ** mesh, unsigned int N, unsigned int M, double h, double k, double D) {
    int n = (N-3)*(M-2);            /* Tamaño del vector */
    double *bVector = new double[n];

    /* Valores para los extremos del bVector */
    bVector[0] = f_function(mesh[1][2], mesh[1][0]) 
                + D*(g_function(mesh[1][1],mesh[1][0], N, M, 0, 1)/(h*h) 
                + g_function(mesh[0][2],mesh[0][0], N, M, 1, 0)/(k*k));
    bVector[n-1] = f_function(mesh[M-2][N-2],mesh[M-2][0]) 
                + D*(g_function(mesh[M-2][N-1],mesh[M-2][0], N, M, N-2, M-2)/(h*h) 
                + g_function(mesh[M-1][N-2],mesh[M-1][0], N, M, N-3, M-1)/(k*k));

    /* Calcular valores para los elementos intermedios del bVector */
    int comp = 0;  /* Comp: componente ... es un contador para moverse en las posiciones del arreglo */
    for (int i = 1; i < N-2; i++) {
        for (int j = 1; j < M-1; j++) {
            if (i == 1 && j == M-2){
                bVector[comp] = f_function(mesh[M-2][2], mesh[M-2][0]) 
                                + D*(g_function(mesh[M-2][1],mesh[M-2][0], N, M, 0, M-2)/(h*h) 
                                + g_function(mesh[M-1][2],mesh[M-1][0], N, M, 1, M-1)/(k*k));
            }else if (i == 1 && j > 1 && j < M-2){
                bVector[comp] = f_function(mesh[j][2], mesh[j][0]) 
                                + D*(g_function(mesh[j][1],mesh[j][0], N, M, 0, j)/(h*h));
            }else if (i == N-3 && j == 1 ){
                bVector[comp] = f_function(mesh[1][N-2], mesh[1][0]) 
                                + D*(g_function(mesh[1][N-1],mesh[1][0], N, M, N-2, 1)/(h*h) 
                                + g_function(mesh[0][N-2],mesh[0][0], N, M, N-3, 0)/(k*k));
            }else if (i == N-3 && j > 1 && j < M-2) {
                bVector[comp] = f_function(mesh[j][N-2], mesh[j][0]) 
                                + D*(g_function(mesh[j][N-1],mesh[j][0], N, M, N-2, j)/(h*h));
            }else if (i > 1 && i < N-3 && j == 1){
                bVector[comp] = f_function(mesh[1][i+1], mesh[1][0]) 
                                + D*(g_function(mesh[0][i+1],mesh[0][0], N, M, i, 0)/(k*k));
            }else if (i > 1 && i < N-3 && j == M-2){
                bVector[comp] = f_function(mesh[M-2][i+1], mesh[M-2][0]) 
                                + D*(g_function(mesh[M-1][i+1],mesh[M-1][0], N, M, i, M-1)/(k*k));
            }else if (i > 1 && i < N-3 && j > 1 && j < M-2){
                bVector[comp] = f_function(mesh[j][i+1], mesh[j][0]);
            }  
            comp++;
        }
    }
    return bVector;
}

/* --------------------------------- Función f(x,y) --------------------------------------- */
double f(double x, double y){
    // return x*x + y*y;
    return exp(-x-y);
}

/* -------------------------- Función g(x,y) para las fronteras --------------------------- */
double g(double x, double y, unsigned int N, 
        unsigned int M, int i, int j){
    if (i == N-2){
        return 200*y;
    } else if (j == M-1){
        return 200*x;
    } else {
        return 0.;
    }   
}

/* ---------------------------------Imprimir matriz--------------------------------------- */
void printMatrix(double **ptrMatrix, unsigned int rows, unsigned int cols){
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            cout << setw(8) << ptrMatrix[i][j];
        }
        cout << endl;
    }
}

/* ---------------------------------Imprimir vector--------------------------------------- */
void printVector(double *vector, unsigned int dim){
    for( int i=0 ; i<dim ; i++){
        cout << vector[i] << endl;
    }
}

/* -------------------------------Eliminar memoria de matriz----------------------------- */
void deleteMatrix(double **ptrMatrix, unsigned int dim){
    for (int i = 0; i < dim; i++){
        delete [] ptrMatrix[i];
    }
    delete []ptrMatrix;
}

/* ----------------------------- Eliminar memoria matriz de tipo entero ------------------------ */
void deleteIntMatrix(unsigned int ** ptrMatrix, unsigned int n_row){
    for (int i = 0; i < n_row; i++){
        delete [] ptrMatrix[i];
    }
    delete [] ptrMatrix;
        
}

/* -------------------------------Guardar matriz en formato CSR----------------------------- */
void writeMatrix_CSR(const char *fileName, double ** ptrMatrix, unsigned int n_rows, unsigned int n_cols){
    int cont;     /* contador */
    ofstream fcout(fileName);
    fcout << n_rows << " ";
    fcout << n_cols << endl;

    for (int i = 0; i < n_rows; i++){
        cont = 0;
        for (int j = 0; j < n_cols; j++){
            if (ptrMatrix[i][j] != 0.){
                cont++;
            }
        }
        fcout << cont << endl;
        for (int j = 0; j < n_cols; j++){
            if (ptrMatrix[i][j] != 0){
                fcout << j << " ";
                fcout << ptrMatrix[i][j] << " ";
            }
        }
        fcout << endl;
    }
    fcout.close();
}

/* --------------------------------- Lectura de Matriz CSR ------------------------------------ */

void readMatrix_CSR( const char *fileName, unsigned int *n_rows , unsigned int *n_cols ,
                     unsigned int **n_Z , double ***AA , unsigned int ***JA ){
    unsigned int *temp_n_Z , **temp_JA;
    double **temp_AA;
    ifstream fcin( fileName );
    fcin>> *n_rows;
    fcin>> *n_cols;
    temp_n_Z = new unsigned int[ *n_rows ];
    temp_AA = new double *[ *n_rows ];
    temp_JA = new unsigned int *[ *n_rows ];
    for( int i=0; i < *n_rows ; i++ ){
        fcin>> temp_n_Z[ i ];
        temp_AA[i] =  new double[ temp_n_Z[ i ] ];
        temp_JA[i] =  new unsigned int[ temp_n_Z[ i ] ];
        for( int j=0 ; j < temp_n_Z[ i ] ; j++ ){
            fcin>> temp_JA[i][j];
            fcin>> temp_AA[i][j];
        }
    }
    fcin.close();
    *n_Z = temp_n_Z;
    *JA = temp_JA;
    *AA = temp_AA;
}

/*---------------------------- Guardar datos para graficarla solución-------------------------*/
void writeSolve(const char *fileName, double **mesh, double *solveVector, unsigned int N, unsigned int M){
    int count = 0;
    ofstream fcout(fileName);
    for (int i = 0; i < M; i++){
        for (int j = 0; j < N-1; j++){
            if (i == 0 || j == 0){
                fcout << mesh[j][i+1] << " " << mesh[j][0] << " " << g(mesh[j][i+1],mesh[j][0],N, M, i, j) << endl;
                // count++;
            }else if (i == M-1 || j == N-2){
                fcout << mesh[j][i+1] << " " << mesh[j][0] << " " << g(mesh[j][i+1],mesh[j][0],N, M, i, j) << endl;
                // count++;
            } else{
                fcout << mesh[j][i+1] << " " << mesh[j][0] << " " << solveVector[count] << endl;
                count++;
            }
        }
    }
    fcout.close();
}

/* ------------------------------------------------------------------------------------------ */
                                /* METODOS DE SOLUCION DIRECTOR */
/* ------------------------------------------------------------------------------------------ */

/*----------------------------------- ELIMINACIÓN GAUSSIANA -----------------------------------*/
void gaussElim(double **Amatrix, double *bvector, unsigned int dim){
    for (int k = 0; k < dim - 1; k++){
        for (int i = k+1; i < dim; i++){
            double temp = Amatrix[i][k]/Amatrix[k][k];
            for (int j = k+1; j < dim; j++){
                Amatrix[i][j] -= temp*Amatrix[k][j];
            }
            bvector[i] -= temp*bvector[k];
        }
    }
}
/*------------------------------------ Sustitución hacia atrás -------------------------------*/
double *solve_ATS(double **Amatrix, double *bvector, unsigned int n ){
    double *ptrXvector = new double [n];
    ptrXvector[n-1] = bvector[n-1]/Amatrix[n-1][n-1];
    for (int i = n-2; i >= 0; i--){
        ptrXvector[i] = 0;
        for (int k = i+1; k < n; k++){
            ptrXvector[i] += Amatrix[i][k]*ptrXvector[k];
        }
        ptrXvector[i] = (bvector[i] - ptrXvector[i])/Amatrix[i][i];
    }

    return ptrXvector;
}

/*-------------------------------- FACTORIZACIÓN LU ---------------------------------*/
void factLU(double **Amatrix, unsigned int dim){
    for (int k = 0; k < dim - 1; k++){
        for (int i = k+1; i < dim; i++){
            double temp = Amatrix[i][k]/Amatrix[k][k];
            for (int j = k+1; j < dim; j++){
                Amatrix[i][j] -= temp*Amatrix[k][j];
            }
            Amatrix[i][k] = temp;
        }
    }
}

/*--------------------- Sustitución hacia adelante con 1s en la diagonal --------------------*/
double *solve_ATI_LU(double **Amatrix, double *bvector, unsigned int n ){
    double *ptrXvector = new double [n];
    ptrXvector[0] = bvector[0];
    for (int i = 1; i < n; i++){
        ptrXvector[i] = 0;
        for (int k = 0; k < i; k++){
            ptrXvector[i] += Amatrix[i][k]*ptrXvector[k];
        }
        ptrXvector[i] = bvector[i] - ptrXvector[i];
    }
    return ptrXvector;
}

/*-------------------------------- FACTORIZACIÓN DE CHOLESKY ---------------------------------*/
void cholesky(double **Amatrix, unsigned int dim){
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < i; j++){
            for (int k = 0; k < j; k++){
                Amatrix[i][j] -= Amatrix[i][k]*Amatrix[j][k];
            }
            Amatrix[i][j] /= Amatrix[j][j];
            Amatrix[j][i] = Amatrix[i][j]; // para guardar L transpuesta en la parte superior de A
        }
        for (int k = 0; k < i; k++){
            Amatrix[i][i] -= Amatrix[i][k]*Amatrix[i][k];
        }
        Amatrix[i][i] = sqrt(Amatrix[i][i]);
    }
}

/*-------------------------------- Sustitución hacia adelante ---------------------------------*/
double *solve_ATI(double **Amatrix, double *bvector, unsigned int n ){
    double *ptrXvector = new double [n];
    ptrXvector[0] = bvector[0]/Amatrix[0][0];
    for (int i = 1; i < n; i++){
        ptrXvector[i] = 0;
        for (int k = 0; k < i; k++){
            ptrXvector[i] += Amatrix[i][k]*ptrXvector[k];
        }
        ptrXvector[i] = (bvector[i] - ptrXvector[i])/Amatrix[i][i];
    }
    return ptrXvector;
}

/*------------------------ Factorización QR mediante Gram-Schmidt------------------------*/
double **qr_gramSchmidt(double ***Amatrix, unsigned int m_rows, unsigned int n_cols){
    double **Rmatrix = new double *[n_cols];
    double **matrixTemp;
    int k = 0;    /* contador */
    
    matrixTemp = transposedMatrix(*Amatrix, m_rows, n_cols); /* transpuesta de A  para acceder*/
                                                             /* facilmente a las columnas de A */
    for (int j = 0; j < n_cols; j++){
        Rmatrix[j] = new double [n_cols];
        for (int i = 0; i < n_cols; i++){
            Rmatrix[j][i] = 0.;
        }
        
        for (k = 0; k < j; k++){
            Rmatrix[j][k] = scalarProd(matrixTemp[k], matrixTemp[j], m_rows);
            for (int n = 0; n < m_rows; n++){
                matrixTemp[j][n] -= Rmatrix[j][k]*matrixTemp[k][n];
            } 
        }

        Rmatrix[k][k] = sqrt(scalarProd(matrixTemp[j], matrixTemp[j], m_rows));
        for (int n = 0; n < m_rows; n++){
            matrixTemp[j][n] = matrixTemp[j][n]/Rmatrix[k][k];
        }  
    }
    Rmatrix = transposedMatrix(Rmatrix, n_cols, n_cols);
    *Amatrix = transposedMatrix(matrixTemp, n_cols, m_rows);
    deleteMatrix(matrixTemp, m_rows);
    return Rmatrix;
}

/*----------------------------- transpuesta de una matriz -----------------------------------*/
double **transposedMatrix(double ** ptrMatrix, unsigned int m_rows, unsigned int n_cols){
    double **transMatrix = new double *[n_cols];
    for (int i = 0; i < n_cols; i++){
        transMatrix[i] = new double [m_rows];
        for (int j = 0; j < m_rows; j++){
            transMatrix[i][j] = ptrMatrix[j][i];
        }
    }
    return transMatrix;
}

/*----------------------------- Producto escalar ---------------------------------------*/
double scalarProd(double *xVector, double *yVector, unsigned int dim){
    double prod = 0.;
    for (int i = 0; i < dim; i++){
        prod += xVector[i]*yVector[i];
    }
    return prod;
}

/*--------------------------- Producto Matriz-Vector -----------------------*/
double *vectMatrixProd(double **Amatrix, double * xVector, unsigned int dim){
    double * vectTemp = new double [dim];
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            vectTemp[i] += Amatrix[i][j]*xVector[j];
        }
    }
    return vectTemp;
}

/* ------------------------------------------------------------------------------------------ */
                               /* METODOS ITERATIVOS CLASICOS */
/* ------------------------------------------------------------------------------------------ */

/* ------------------------------- Método de Jacobi----------------------------------------- */
double *jacobi(double **AAmatrix, unsigned int **JAmatrix, unsigned int *n_Zmatrix, double *bVector, 
                double *x0Vector, double tol, unsigned int dim){
    double *x_kVector = new double [dim];
    double *ptrVectorTemp = new double [dim];
    // int contador = 1;

    do{
        // contador++;
        /* Para hacer una copia del vector inicial. Esto me evita fuga de memoria */
        for (int i = 0; i < dim; i++){
            ptrVectorTemp[i] = x0Vector[i];
        }
        
        for (int i = 0; i < dim; i++){
            x_kVector[i] = 0.;
            for (int j = 0; j < dim; j++){
                if (j != i){
                    x_kVector[i] -= get_aij(AAmatrix, JAmatrix, n_Zmatrix, i, j) * ptrVectorTemp[j];
                }
            }
            x_kVector[i] = (x_kVector[i] + bVector[i]) / get_aij(AAmatrix, JAmatrix,n_Zmatrix, i,i);
        }
        
        /* Copia del vector x_k en x_0. Evito fuga de memoria y ahora x_0 es mi iteración k
            para luego almacenar la iteración k+1 en x_k*/
        for (int i = 0; i < dim; i++){
            x0Vector[i] = x_kVector[i];
        }

    } while (vectNorm( substractVect(x_kVector, ptrVectorTemp, dim), dim ) >= tol);
    
    // cout << "Iteraciones: " << contador << endl;
    delete [] ptrVectorTemp;
    return x_kVector;
}

/* ------------------------------- Método de Gauss-Seidel ----------------------------------------- */
double *gauss_seidel(double **AAmatrix, unsigned int **JAmatrix, unsigned int *n_Zmatrix, double *bVector,
                double *x0Vector, double tol, unsigned int dim){
    double *x_kVector = new double [dim];
    double *ptrVectorTemp = new double [dim];
    double temp;
    // int contador = 1;

    do{
        // contador++;
        /* Para hacer una copia del vector inicial. Esto me evita fuga de memoria */
        for (int i = 0; i < dim; i++){
            ptrVectorTemp[i] = x0Vector[i];
        }

        for (int i = 0; i < dim; i++){
            x_kVector[i] = 0.;
            for (int j = i+1; j < dim; j++){
                x_kVector[i] -= get_aij(AAmatrix, JAmatrix, n_Zmatrix, i, j) * ptrVectorTemp[j];
            }
            for (int j = 0; j < i; j++){
                temp = x_kVector[j];
                x_kVector[i] -= get_aij(AAmatrix, JAmatrix, n_Zmatrix, i, j) * temp;
            }
            x_kVector[i] = (x_kVector[i] + bVector[i]) / get_aij(AAmatrix, JAmatrix,n_Zmatrix, i,i);
        }

        /* Copia del vector x_k en x_0. Evito fuga de memoria y ahora x_0 es mi iteración k
            para luego almacenar la iteración k+1 en x_k*/
        for (int i = 0; i < dim; i++){
            x0Vector[i] = x_kVector[i];
        }

    } while (vectNorm( substractVect(x_kVector, ptrVectorTemp, dim), dim ) >= tol);

    // cout << "Iteraciones: " << contador << endl;
    delete [] ptrVectorTemp;
    return x_kVector;
}

/* ------------------------------- Método del Gradiente conjugado ---------------------------------- */
double *grad_conj(double **AAmatrix, unsigned int **JAmatrix, unsigned int *n_Zmatrix, double *bVector, 
                double *x0Vector, double tol, unsigned int dim){
    double *r_vector = new double [dim];   /* sirve para almacenar los p_i y los r_i */
    double *p_vector = new double [dim];  
    double *w_vector = new double [dim];
    double *x_vector = new double [dim];
    double alpha = 0., beta = 0.;
    // int contador = 0;

    /* Asignando r0 vector a r_vector */
    r_vector = substractVect(bVector, vectMatrixProd_CSR(AAmatrix, JAmatrix, n_Zmatrix, x0Vector, dim), dim);
    while (scalarProd(r_vector, r_vector, dim) >= tol){
        for (int i = 0; i < dim; i++){  /* copia de r0_vector a p_vector */
            p_vector[i] = r_vector[i];
        }
        w_vector = vectMatrixProd_CSR(AAmatrix,JAmatrix, n_Zmatrix, p_vector, dim);
        alpha = scalarProd(r_vector, r_vector, dim) / scalarProd(p_vector, w_vector, dim);
        
        for (int i = 0; i < dim; i++){
            x_vector[i] = x0Vector[i] + alpha*p_vector[i];  /* calculando el x_i+1  */
            r_vector[i] = p_vector[i] - alpha*w_vector[i];  /* calculando el r_i+1 */
        }
        
        beta = scalarProd(r_vector, w_vector, dim) / scalarProd(p_vector, w_vector, dim);
        for (int i = 0; i < dim; i++){
            p_vector[i] = r_vector[i] - beta*p_vector[i];  /* calculando el p_i+1  */
            x0Vector[i] = x_vector[i];                     /* copiando el x_i en x0 */
        }
        // contador++;
    }
    // cout << "Iteraciones: " << contador << endl;
    delete [] r_vector;
    delete [] w_vector;
    delete [] p_vector;
    
    return x_vector;
}

/* ------------------------------- Vector inicial x0----------------------------------------- */
double *startVector(unsigned int n_rows){
    double *x0Vector = new double [n_rows];
    for (int i = 0; i < n_rows; i++){
        x0Vector[i] = 0.;
    }
    return x0Vector;
}

/* --------------------------- Obtener el elemento Aij por medio de AA---------------------- */
double get_aij(double **AAmatrix, unsigned int **JAmatrix, unsigned int *n_Zmatrix, 
                    unsigned int i_row, unsigned int j_col){
    double a_ij;  /* elemento a_ij de la matriz */
    int col_k;    /* columna k que se obtendrá de JA */ 
    if (n_Zmatrix[i_row] != 0){
        col_k = getColumnJA(JAmatrix, n_Zmatrix, i_row, j_col);
        if (col_k != -1){
            a_ij = AAmatrix[i_row][col_k];
        } else {
            return 0.;
        }
    } else {
        return 0.;
    }
    return a_ij;
}

/* --------------------------------- Obtener la columna de AA ------------------------------------ */
int getColumnJA(unsigned int **JAmatrix, unsigned int *n_Zmatrix , unsigned int i_row,  unsigned int j){
    /*  i_row: sirve para posicionarse en la fila donde se desea buscar 
        j: elemento a comparar en JA  */

    unsigned int colAA;  /* colAA: es la columna donde se buscará el elemento en AA */
    for (int i = 0; i < n_Zmatrix[i_row]; i++){
        if (JAmatrix[i_row][i] == j){
            colAA = i;
            return colAA;
        } else {
            colAA = -1;  /* Esto nos indica que el valor de la columna buscada no existe */
        }
    }
    return colAA;
}

/* ------------------------------- Norma euclidea ----------------------------------------- */
double vectNorm(double *xVector, unsigned int dim){
    double temp = 0;
    for (int i = 0; i < dim; i++){
        temp += xVector[i]*xVector[i];
    }
    return sqrt(temp);
}

/*------------------------------- Resta de dos vectores -----------------------------------*/
double *substractVect(double *xVector, double *yVector, unsigned int dim){
    double *zVector = new double [dim];
    for (int i = 0; i < dim; i++){
        zVector[i] = xVector[i] - yVector[i];
    }
    return zVector;
}

/*--------------------------- Producto Matriz-Vector en CSR -----------------------*/
double *vectMatrixProd_CSR(double **AAmatrix, unsigned int **JAmatrix, unsigned int *n_Zmatrix, 
                       double * xVector, unsigned int dim){
    double *vectTemp = new double [dim];
    for (int i = 0; i <dim; i++){
        vectTemp[i] = 0.;
        for (int j = 0; j < dim; j++){
            vectTemp[i] += get_aij(AAmatrix, JAmatrix, n_Zmatrix, i, j)* xVector[j];
        }
    }
    return vectTemp;
}