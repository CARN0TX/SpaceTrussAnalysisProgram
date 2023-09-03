/* Matrix Toolbox

  NOTES
   (1) Illustrates the development of a matrix toolbox
       based on the CVector and CMatrix classes
       developed earlier.
   (2) Exceptions should be thrown for invalid operations

*/
#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <exception>
#include "ArrayContainersEXH.h"

template <class T>
class CMatToolBox
{
    const int ELEMENTSPERLINE = 4;  // # of vector/matrix elements per line
    const int FW = 16;              // field width

    public:
        enum class Error
        { VECERR_ADD, VECERR_SUBTRACT, VECERR_DOTPRODUCT, VECERR_CROSSPRODUCT,
          VECERR_NORMALIZE,
          MATERR_ADD, MATERR_SUBTRACT, MATERR_MULTIPLY, MATERR_DETERMINANT,
          MATERR_TRANSPOSE, MATERR_MATMULTVEC, MATERR_LUFACTORIZATION,
          MATERR_LUSOLVE, MATERR_LDLTFACTORIZATION, MATERR_LDLTSOLVE,
          MATERR_GAUSSELIMINATION, MATERR_SINGULARMATRIX, MATERR_NOTPOSDEFMATRIX,
          MATERR_RESIDUALVECTOR, UNSUPPORTEDOPERATION, MATERR_NTSQRMAT
        };
        CMatToolBox ();
        ~CMatToolBox ();

        // vector-related functions
        void Display (const std::string& strMessage,
                      const CVector<T>& A) const;
        void Add (const CVector<T>& A, const CVector<T>& B, 
                  CVector<T>& C);
        void Subtract (const CVector<T>& A, const CVector<T>& B,
                       CVector<T>& C);
        T DotProduct (const CVector<T>& A,
                      const CVector<T>& B);
        void Normalize (CVector<T>& A);
        void Scale (CVector<T>& A, const T factor);
        T    MaxValue (const CVector<T>& A) const;
        T    MinValue (const CVector<T>& A) const;
        T    TwoNorm (const CVector<T>& A);
        T    MaxNorm (const CVector<T>& A) const;
        void CrossProduct (const CVector<T>& A,
                           const CVector<T>& B, CVector<T>& C);

        // matrix-related functions
        void Display (const std::string& strMessage,
                      const CMatrix<T>& A) const;
        void Multiply (const CMatrix<T>& A,
                       const CMatrix<T>& B, CMatrix<T>& C);
        T    MaxNorm (const CMatrix<T>& A) const;
        void Transpose (const CMatrix<T>& A,
                        CMatrix<T>& B);
        void MatMultVec (const CMatrix<T>& A,
                         const CVector<T>& x,
                         CVector<T>& b);
        void AxEqb (CMatrix<T>& A,
                    CVector<T>& x,
                    CVector<T>& b,
                    const T TOL);
        void LDLTFactorization (CMatrix<T>& A, const T TOL);
        void LDLTSolve (const CMatrix<T>& A, CVector<T>& x,
                        const CVector<T>& b);

        // helper function
        void ResidualVector (const CMatrix<T>& A, const CVector<T>& x,
                             const CVector<T>& b, CVector<T>& R,
                             T& AbsError, T& RelError);
        bool IsEqual (const T d1, const T d2, T TOL = TOLDEF) const;
        bool IsEqual (const CMatrix<T>& dMA,
                      const CMatrix<T>& dMB, T TOL = TOLDEF) const;
        bool IsEqual (const CVector<T>& dVA,
                      const CVector<T>& dVB, T TOL = TOLDEF) const;
        void GetFLOPStats (double& dAS, double& dM, double& dD) const;

    private:
        // these are declared double to avoid integer overflow
        double m_dASOP; // # of floating point additions and subtractions
        double m_dMOP;  // # of floating point multiplications
        double m_dDOP;  // # of floating point divisions
        void ErrorHandler (CMatToolBox<T>::Error err) const;

    protected:
};

// ctor
template <class T>
CMatToolBox<T>::CMatToolBox ()
// ==================================================================
// Function: default constructor
//    Input: none
//   Output: none
// ==================================================================
{
    m_dASOP = m_dMOP = m_dDOP = 0.0;
}

// dtor
template <class T>
CMatToolBox<T>::~CMatToolBox ()
// ==================================================================
// Function: destructor
//    Input: none
//   Output: none
// ==================================================================
{
}

// ---------------------------------------------------------------
// ------------------------ vector functions ---------------------
// ---------------------------------------------------------------
template <class T>
void CMatToolBox<T>::Display (const std::string& strMessage,
                              const CVector<T>& A) const
// ==================================================================
// Function: displays a message and the elements of a vector
//    Input: message string and the vector 
//   Output: None
// ==================================================================
{
    std::cout << '\n' << strMessage << '\n';
    std::cout.setf(std::ios::left);
    for (int i=1; i <= A.GetSize(); i++)
    {
        std::cout << "(" << i << ") "
                  << std::setw(FW) << A(i) << " ";
        if ((i % ELEMENTSPERLINE) == 0)
            std::cout << '\n';
    }
}

template <class T>
void CMatToolBox<T>::Add (const CVector<T>& A, const CVector<T>& B,
                          CVector<T>& C)
// ==================================================================
// Function: adds two vectors and stores the result in the
//           third vector C = A + B
//    Input: vectors A and B 
//   Output: vector C
// ==================================================================
{
    // check for incompatible vectors
    int n = A.GetSize();
    if (n != B.GetSize() || n != C.GetSize())
        ErrorHandler (Error::VECERR_ADD);

    // add
    for (int i=1; i <= n; i++)
        C(i) = A(i) + B(i);
    m_dASOP += static_cast<double>(n);
}

template <class T>
void CMatToolBox<T>::Subtract (const CVector<T>& A,
                               const CVector<T>& B, CVector<T>& C)
// ==================================================================
// Function: subtracts one vector from another and stores the result
//           in the third vector C = A - B
//    Input: vectors A and B 
//   Output: vector C
// ==================================================================
{
    int n = A.GetSize();
    if (n != B.GetSize() || n != C.GetSize()) // check for incompatible vectors
        ErrorHandler(Error::VECERR_SUBTRACT);

    // Subtract
    for (int i = 1; i <= n; i++)
    {
        C(i) = A(i) - B(i);
    }
    m_dASOP += static_cast<double>(n);
}

template <class T>
T CMatToolBox<T>::DotProduct (const CVector<T>& A,
                              const CVector<T>& B)
// ==================================================================
// Function: computes the dot product of two vectors such that
//           product = A dot B
//    Input: vectors A and B 
//   Output: product 
// ==================================================================
{
    int n = A.GetSize();
    if (n != B.GetSize()) // check for incompatible vectors
        ErrorHandler(Error::VECERR_DOTPRODUCT);

    // dot product
    T product = 0;
    for (int i = 1; i <= n; i++)
        product += A(i) * B(i);
    m_dASOP += static_cast<double>(n);
    return product;
}

template <class T>
void CMatToolBox<T>::Normalize (CVector<T>& A)
// ==================================================================
// Function: normalizes a vector
//    Input: vector A 
//   Output: normalized vector A 
// ==================================================================
{
    int n = A.GetSize();
    T magnitude = 0;
    for (int i = 1; i <= n; i++) {
        magnitude += A(i) * A(i);
    }
    if(magnitude==0) //in case norm of the vector is zero
        ErrorHandler (Error::VECERR_NORMALIZE);
    for (int i = 1; i <= n; i++) {
        A(i)=A(i)/sqrt(magnitude); //unit vector
    }
}

template <class T>
void CMatToolBox<T>::Scale (CVector<T>& A, T c)
// ==================================================================
// Function: scales a vector by a constant c such that A = c A
//    Input: vector A and constant c 
//   Output: scaled vector A
// ==================================================================
{
    int n = A.GetSize();
    for(int i=1;i<=n;i++)
    {
        A(i) = c * A(i); //multiply each element of A by c
    }
}

template <class T>
T CMatToolBox<T>::MaxValue (const CVector<T>& A) const
// ==================================================================
// Function: finds the largest value among all the elements in A
//    Input: vector A 
//   Output: return value is the largest element in A
// ==================================================================
{
    T maxValue = FLT_MIN; //assigning a minimum value element to compare the 1st element of A
    int n = A.GetSize();
    for (int i = 1; i <= n; i++) {
        maxValue = std::max(maxValue, A(i)); // checking value of each element of A to be greater than the previous
    }
    return maxValue;
}

template <class T>
T CMatToolBox<T>::MinValue (const CVector<T>& A) const
// ==================================================================
// Function: finds the smallest value among all the elements in A
//    Input: vector A 
//   Output: return value is the smallest element in A
// ==================================================================
{
    T minValue = FLT_MAX;//assigning a maximum value element to compare the 1st element of A
    int n = A.GetSize();
    for (int i = 1; i <= n; i++) {
        minValue = std::min(minValue, A(i));// checking value of each element of A to be lesser than the previous
    }
    return minValue;
}

template <class T>
T CMatToolBox<T>::TwoNorm (const CVector<T>& A)
// ==================================================================
// Function: computes the two norm of vector A
//    Input: vector A 
//   Output: return value is the two-norm
// ==================================================================
{
    int n = A.GetSize();
    T magnitude = 0;
    for (int i = 1; i <= n; i++) {
        magnitude += A(i) * A(i); // sqaure of magnitude of the vector
    }
    return sqrt(magnitude);
}

template <class T>
T CMatToolBox<T>::MaxNorm (const CVector<T>& A) const
// ==================================================================
// Function: computes the max norm of vector A
//    Input: vector A 
//   Output: return value is the max norm
// ==================================================================
{
    int n = A.GetSize();
    T maxNorm = 0;
    for (int i = 1; i <= n; i++) {
        maxNorm = std::max(maxNorm, abs(A(i)));//absolute maximum value among all elements of the vector
    }
    return maxNorm;
}

template <class T>
void CMatToolBox<T>::CrossProduct (const CVector<T>& A,
                                   const CVector<T>& B,
                                   CVector<T>& C)
// ==================================================================
// Function: computes the cross-product of two vectors and stores the
//           result in the third vector such that C = A x B
//           (3-dimensional space)
//    Input: vectors A, B and C
//   Output: vector C
// ==================================================================
{
    int p = A.GetSize();
    int q = B.GetSize();
    int r = C.GetSize();
    if ((p != q) || (q != r) || (p != r) || p > 3 || q > 3) // check for incompatible vectors or vectors with 3directional+ space
    {
        ErrorHandler(Error::VECERR_CROSSPRODUCT);
    }
    // A × B = (bz – cy)i + (cx – az)j + (ay – bx)k
    C(1) = (A(2) * B(3)) - (A(3) * B(2));
    C(2) = (A(3) * B(1)) - (A(1) * B(3));
    C(3) = (A(1) * B(2)) - (A(2) * B(1));
}

// ---------------------------------------------------------------
// ------------------------ matrix functions ---------------------
// ---------------------------------------------------------------
template <class T>
void CMatToolBox<T>::Display (const std::string& strMessage,
                              const CMatrix<T>& A) const
// ==================================================================
// Function: displays a message and the elements of a matrix
//           rowwise
//    Input: message string and the matrix
//   Output: None
// ==================================================================
{
    std::cout << '\n' << strMessage << '\n';
    std::cout.setf(std::ios::left);
    for (int i=1; i <= A.GetRows(); i++)
    {
        int nC = 0;
        for (int j=1; j <= A.GetColumns(); j++)
        {
            ++nC;
            std::cout << "(" << i << "," << j << ") "
                      << std::setw(FW) << A(i,j) << " ";
            if ((nC % ELEMENTSPERLINE) == 0)
                std::cout << '\n';
        }
        std::cout << '\n';
    }
}

template <class T>
void CMatToolBox<T>::Multiply (const CMatrix<T>& A,
                               const CMatrix<T>& B, CMatrix<T>& C)
// ==================================================================
// Function: multiplies two matrices and stores the result
//           in the third matrix C = A * B
//    Input: matrices A and B 
//   Output: matrix C
// ==================================================================
{
    int ra = A.GetRows();
    int ca = A.GetColumns();
    int rb = B.GetRows();
    int cb = B.GetColumns();
    if ((ca != rb) || (C.GetRows()) != ra || C.GetColumns() != cb) // check for incompatible matrices
    {
        ErrorHandler(Error::MATERR_MULTIPLY);
    }
    for (int i = 1; i <= ra; i++) {
        for (int j = 1; j <= cb; j++) {
            C(i, j) = 0;//initializing C with each element equals zero
        }
    }
    for (int i = 1; i <= ra; i++) {
        for (int j = 1; j <= cb; j++) {
            for (int k = 1; k <= ca; k++) {
                C(i, j) += A(i, k) * B(k, j);
            }
        }
    }
}

template <class T>
T CMatToolBox<T>::MaxNorm (const CMatrix<T>& A) const
// ==================================================================
// Function: computes the max norm of matrix A
//    Input: matrix A 
//   Output: return value is the max norm
// ==================================================================
{
    //max norm of a matrix is the maximum row sum of a matrix taken absolute values
    int ra = A.GetRows();
    int ca = A.GetColumns();
    T maxNorm = 0;
    for (int i = 1; i <= ra; i++){
        T rowSum = 0;
        for (int j = 1; j <= ca; j++) {
            rowSum += abs(A(i, j));
        }
        maxNorm = std::max(maxNorm, rowSum); //
    }
    return maxNorm;
}

template <class T>
void CMatToolBox<T>::Transpose (const CMatrix<T>& A,
                                CMatrix<T>& B)
// ==================================================================
// Function: computes the transpose of a matrix and stores the result
//           in another matrix B = A(T)
//    Input: matrices A and B
//   Output: matrix B
// ==================================================================
{
    int ra = A.GetRows();
    int ca = A.GetColumns();
    int rb = B.GetRows();
    int cb = B.GetColumns();
    if ((ca != rb) || (cb!=ra)) //incompatibility check
    {
        ErrorHandler(Error::MATERR_TRANSPOSE);
    }
    for (int i = 1; i <= ra; i++) {
        for (int j = 1; j <= ca; j++) {
            B(j, i) = A(i, j); //interchange i,j position of the element
        }
    }
}

template <class T>
void CMatToolBox<T>::MatMultVec (const CMatrix<T>& A,
                                 const CVector<T>& x,
                                 CVector<T>& b)
// ==================================================================
// Function: multiplies a matrix and a vector and stores the result
//           in a vector b = A * x
//    Input: vectors A and x 
//   Output: vector b
// ==================================================================
{
    int ra = A.GetRows();
    int ca = A.GetColumns();
    int rx = x.GetSize(); //treating vector as a column matrix
    if (ca != rx || b.GetSize() != ra) //incompatibility check
    {
        ErrorHandler(Error::MATERR_MATMULTVEC);
    }
    for (int j = 1; j <= rx; j++) {
        b(j) = 0;//initializing resulting vector
    }
    for (int i = 1; i <= ra; i++) {
        for (int j = 1; j <= ca;j++) {
            b(i) += A(i, j) * x(j); //row i, column jth elements of A multiplied to jth row element of x and assigned to ith elemnt of b
        }
    }
}

template <class T>
void CMatToolBox<T>::AxEqb (CMatrix<T>& A,
                            CVector<T>& x,
                            CVector<T>& b,
                            T TOL)
// ==================================================================
// Function: solves A x = b using Gaussian Elimination Technique
//           (this version does NOT perform partial pivoting)
//    Input: Matrices A and b
//   Output: Matrix x
//           throws an exception for a singular A matrix.
// ==================================================================
{
    int i, j, k;   // loop indices
    T c;           // multiplier (Step 4)

    // number of equations to solve
    int n = A.GetRows();

    // x initially contains b
    x = b;

    // forward elimination
    for (k=1; k <= n-1; k++)              // Step 1
    {
        if (fabs(A(k,k)) <= TOL)          // Step 2
            ErrorHandler (Error::MATERR_GAUSSELIMINATION);
        for (i=k+1; i <= n; i++)          // Step 3
        {
            c = A(i,k)/A(k,k);            // Step 4
            for (j=k+1; j <= n; j++)      // Step 5
                A(i,j) -= c * A(k,j);     // Step 6
            x(i) -= c * x(k);             // Step 8
        }                                 // Step 9
        int nC = n-k;
        if (nC > 0)
        {
            m_dDOP += static_cast<double>(nC);
            m_dMOP += static_cast<double>(nC*nC);
            m_dASOP += static_cast<double>(nC+nC*nC);
        }
    }                                     // Step 10 

    // back substitution
    if (fabs(A(n,n)) <= TOL)              
        ErrorHandler (Error::MATERR_GAUSSELIMINATION);
    x(n) /= A(n,n);                       // Step 11

    for (i=n-1; i >= 1; i--)              // Step 12
    {
        T sum = T(0);
        for (j=i+1; j <= n; j++)
            sum += A(i,j) * x(j);         // Step 13
        if ((n-i) > 0)
        {
            m_dASOP += static_cast<double>(n-i);
            m_dMOP += static_cast<double>(n-i);
        }
        x(i) = (x(i) - sum)/A(i,i);       // Step 14
    }                                     // Step 15
    m_dASOP += static_cast<double>(n);
    m_dDOP += static_cast<double>(n+1);
}

template <class T>
void CMatToolBox<T>::LDLTFactorization (CMatrix<T>& A, T TOL)
// ==================================================================
// Function: carries out LDL(T) factorization of matrix A
//           A is replaced with L and D. A is a symmetric matrix.
//    Input: matrix A and tolerance value to detect singular A
//   Output: matrix A 
// ==================================================================
{
    int r = A.GetRows();
    int c = A.GetColumns();
    CMatrix<T> L(r, c); //L is the lower unit triangular matrix
    CMatrix<T> LT(r, c); //L is the UPPER unit triangular matrix
    CVector<T> D(r); //diagonal matrix(since it has elements equal to number of rows, vector is choosen for the sake of simplicity and memory allocation)
    if(r!=c)
        ErrorHandler (Error::MATERR_NTSQRMAT);
    for (int i = 1; i <= r; i++) //initializing L, D(vector) matrices
    {
        D(i) = 0.0;
        for (int j = 1; j <= c; j++)
        {
            LT(i, j) = A(i,j);
        }
    }
    //Obtaining LT matrix
    for (int i = 1; i <= r; i++) {
        for (int j = i + 1; j <= r; j++) {
            T fac = LT(j, i) / LT(i, i);
            for (int k = 1; k <= c; k++) {
                LT(j, k) -= fac * LT(i, k);
            }
        }
        D(i) = LT(i, i); //assigning diagonal elements to D matrix
        if (fabs(D(i)) <= TOL)
            ErrorHandler(Error::MATERR_LDLTFACTORIZATION);
    }
    //making LT upper UNIT triangular matrix
    for (int i = 1; i <= r; i++) {
        for (int j = 1; j <= c; j++) {
            LT(i, j) = LT(i, j) / D(i);
        }
    }
    Transpose(LT, L);
    // We have obtained L,D,L(T) matrices thus LDL(T) factorization of matrix A
    //now cholesky decomposition for LDL(T) solve assigining terms to matrix A
    for (int i = 1; i <= r; i++) {
        D(i) = pow(D(i), 0.5);
        L(i, i) *= D(i);
    }
    //Combining L and D^0.5 matrices to form CC(T) matrix A
    for (int i = 1; i <= r; i++) {
        for (int j = 1; j <= c; j++) {
            A(i,j)=L(i, j);
        }
    }
}

template <class T>
void CMatToolBox<T>::LDLTSolve (const CMatrix<T>& A,
                                CVector<T>& x,
                                const CVector<T>& b)
// ==================================================================
// Function: carries out forward and backward substitution so as to
//           solve A x = b. A contains L and D terms.
//    Input: matrix A, vectors x and b
//   Output: vector x 
// ==================================================================
{
    int r = A.GetRows();
    int c = A.GetColumns();
    if (r != x.GetSize() || b.GetSize() != c)
        ErrorHandler(Error::MATERR_LDLTSOLVE);
    //Forward substitution Cy=b, solve for y
    CVector<T> y(r);
    for (int i = 1; i <= r; i++) {
        float alpha = b(i);
        for (int j = 1; j <= i; j++) {
            alpha -= A(i, j) * y(j);
        }
        y(i) = alpha / A(i, i);  
    }

    //Backward substitution C(T)x=y, solve for x
    for (int i = r; i >= 1; i--) {
        float alpha = y(i);
        for (int j = 1; j <= i; j++) {
            alpha -= A(i, j) * x(j);
        }
        x(i) = alpha / A(i, i);
    }
}

template <class T>
void CMatToolBox<T>::ResidualVector (const CMatrix<T>& A, 
                                     const CVector<T>& x,
                                     const CVector<T>& b,
                                     CVector<T>& R,
                                     T& AbsError, T& RelError)
// ==================================================================
// Function: computes the residual vector arising from the solution
//           of A x = b, i.e., R = A x - b
//    Input: matrix A, vectors x and b, 
//   Output: vector R, abs. and rel. error values via TwoNorm function
// ==================================================================
{
    // check for incompatible sizes
    int n = A.GetRows();
    int nCols = A.GetColumns();
    int nx = x.GetSize();
    int nb = b.GetSize();

    if (n != nx || nCols != nx || nx != nb)
        ErrorHandler (Error::MATERR_RESIDUALVECTOR);

    MatMultVec (A, x, R);
    for (int i = 1; i <= n; i++)
    {
        R(i) -= b(i);
    }
    m_dASOP += static_cast<double>(n);
    AbsError = TwoNorm(R);
    RelError = TwoNorm(R)/TwoNorm(b);
}


template <class T>
void CMatToolBox<T>::GetFLOPStats (double& dAS, double& dM, 
                                   double& dD) const
// ==================================================================
// Function: retrieves floating point operations
//    Input: variables to store +-, * and / operations
//   Output: variables with their values
// ==================================================================
{
    dAS = m_dASOP;
    dM  = m_dMOP;
    dD  = m_dDOP;
}

template <class T>
bool CMatToolBox<T>::IsEqual (const T d1, const T d2, T TOL) const
// ==================================================================
// Function: checks if d1 and d2 are 'nearly' equal
//    Input: d1, d2, tolerance to use
//   Output: true if they are
// ==================================================================
{
    return (abs(d1-d2) <= TOL);
}

template <class T>
bool CMatToolBox<T>::IsEqual (const CMatrix<T>& dMA,
                              const CMatrix<T>& dMB, T TOL) const
// ==================================================================
// Function: checks if matrices A and B are 'nearly' equal
//    Input: A, B, tolerance to use
//   Output: true if they are
// ==================================================================
{
    for (int i=1; i <= dMA.GetRows(); i++)
    {
        for (int j=1; j <= dMA.GetColumns(); j++)
        {
            if (!IsEqual(dMA(i,j), dMB(i,j), TOL))
                return false;
        }
    }

    return true;
}

template <class T>
bool CMatToolBox<T>::IsEqual (const CVector<T>& dVA,
                              const CVector<T>& dVB, T TOL) const
// ==================================================================
// Function: checks if vectors A and B are 'nearly' equal
//    Input: A, B, tolerance to use
//   Output: true if they are
// ==================================================================
{
    for (int i=1; i <= dVA.GetSize(); i++)
    {
        if (!IsEqual(dVA(i), dVB(i), TOL))
            return false;
    }

    return true;
}


// ==================== Error Handler ========================
template <class T>
void CMatToolBox<T>::ErrorHandler (Error err) const
// ---------------------------------------------------------------------------
// Function: channels error message
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    if (err == CMatToolBox<T>::Error::VECERR_ADD)
    {
        throw std::exception ("Vector addition operation cannot take place due to "
                              "incompatible vectors.");
    }
    else if (err == CMatToolBox<T>::Error::VECERR_SUBTRACT)
    {
        throw std::exception("Vector subtraction operation cannot take place due to "
            "incompatible vectors.");
    }
    else if (err == CMatToolBox<T>::Error::VECERR_DOTPRODUCT)
    {
        throw std::exception("Vector dot product operation cannot take place due to "
            "incompatible vectors.");
    }
    else if (err == CMatToolBox<T>::Error::VECERR_NORMALIZE)
    {
        throw std::exception ("Cannot normalize vector since norm is zero.");
    }
    else if (err == CMatToolBox<T>::Error::VECERR_CROSSPRODUCT)
    {
        throw std::exception("Vector cross product operation cannot take place due to "
            "incompatible vectors.");
    }
    else if (err == CMatToolBox<T>::Error::MATERR_MULTIPLY)
    {
        throw std::exception("Matrix multiplication operation cannot take place due to incompatible matrices");
    }
    else if (err == CMatToolBox<T>::Error::MATERR_TRANSPOSE)
    {
        throw std::exception("Matrix transpose operation cannot take place due to incompatible matrices B(T)");
    }
    else if (err == CMatToolBox<T>::Error::MATERR_MATMULTVEC)
    {
        throw std::exception("Matrix vector multiplication operation cannot take place due to incompatible matrix and vector");
    }
    else if (err == CMatToolBox<T>::Error::MATERR_LDLTFACTORIZATION)
    {
        throw std::exception("The diagonal element is too small for LDLT factorization.");
    }
    else if (err == CMatToolBox<T>::Error::UNSUPPORTEDOPERATION)
    {
        throw std::exception ("This functionality is not supported.");
    }
    else if (err == CMatToolBox<T>::Error::MATERR_GAUSSELIMINATION)
    {
        throw std::exception ("Gaussian Elimination: Dependent equations. Diagonal element too small.");
    }
    else if (err == CMatToolBox<T>::Error::MATERR_NTSQRMAT)
    {
        throw std::exception("Cannot proceed with the operation as A is not a square matrix.");
    }
    else
        throw std::exception ("MatrixToolBox: Unknown error");
}
