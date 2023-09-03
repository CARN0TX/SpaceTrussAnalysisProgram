/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-22, S. D. Rajan
All rights reserved

Intermediate Structural Analysis and Design and
Object-Oriented Numerical Analysis via C++

Contains CLocalErrorHandler class implementation.

List of improvements made:

*********************************************/
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cfenv>
#include "LocalErrorHandler.h"

CLocalErrorHandler::CLocalErrorHandler()
// ---------------------------------------------------------------------------
// Function: ctor (initializes the member variables)
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nLineNumber = 0;
}

CLocalErrorHandler::~CLocalErrorHandler()
// ---------------------------------------------------------------------------
// Function: dtor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{

}

void CLocalErrorHandler::ErrorHandler (ERRORCODE& ErrorCode, int n)
// ---------------------------------------------------------------------------
// Function: handles the error messages that have a line number
// Input:    errorcode and line number
// Output:   
// ---------------------------------------------------------------------------
{
    m_nLineNumber = n;
    ErrorHandler (ErrorCode);
}

void CLocalErrorHandler::ErrorHandler (ERRORCODE& ErrorCode)
// ---------------------------------------------------------------------------
// Function: Displays the error message on the error stream
// Input:    Error #
// Output:   none
// ---------------------------------------------------------------------------
{
    std::cerr << '\n';
    std::ostringstream strMessage;

    if (ErrorCode == ERRORCODE::INVALIDNUMNODES)                                 // invalid number of nodes
        std::cerr << "Number of nodes must be >= 2.";
    else if (ErrorCode == ERRORCODE::INVALIDNUMELEMENTS)                         // invalid number of elements
        std::cerr << "Number of elements must be >= 1.";
    else if (ErrorCode == ERRORCODE::INVALIDDEBUGCODE)                           // invalid debug level
        std::cerr << "Debug level must be 0 or 1.";
    else if (ErrorCode == ERRORCODE::INVALIDNODENUM)                             // invalid node number
        std::cerr << "Invalid node number";
    else if (ErrorCode == ERRORCODE::INVALIDELEMENTNUM)                          // invalid element number
        std::cerr << "Invalid element number";
    else if (ErrorCode == ERRORCODE::INVALIDCSAREA)                              // invalid x/s area
        std::cerr << "Area must be positive.";
    else if (ErrorCode == ERRORCODE::INVALIDYM)                                  // invalid E
        std::cerr << "Modulus of elasticity must be positive.";
    else if (ErrorCode == ERRORCODE::UNSTABLETRUSS)                              // unstable structure
        std::cerr << "Unstable truss.";
    else if (ErrorCode == ERRORCODE::INVALIDINPUT)                               // invalid input
        std::cerr << "Invalid input detected in input file.";
    else if (ErrorCode == ERRORCODE::INVALIDNODALFIXITY)                         // invalid fixity code
        std::cerr << "Nodal fixity code must be 0 or 1.";
    else if (ErrorCode == ERRORCODE::MISSINGEND)                                 // missing end statement
        std::cerr << "Missing *END statement in input file.";
    else if (ErrorCode == ERRORCODE::CANNOTOPENIFILE)                            // cannot open input file
        std::cerr << "Cannot open specified input file.";
    else if (ErrorCode == ERRORCODE::CANNOTOPENOFILE)                            // cannot open output file
        std::cerr << "Cannot open specified output file.";
    else if (ErrorCode == ERRORCODE::INVALIDCOMMANDLINE)                         // need 1 or 3 command line arguments
        std::cerr << "Invalid number of command line arguments.";
    else if (ErrorCode == ERRORCODE::MTB_VECTORADDERROR)
        std::cerr << "Vector: addition error.";
    else if (ErrorCode == ERRORCODE::MTB_VECTORSUBTRACTERROR)
        std::cerr << "Vector: subtraction error.";
    else if (ErrorCode == ERRORCODE::MTB_VECTORNORMALIZEERROR)
        std::cerr << "Vector: normalize error.";
    else if (ErrorCode == ERRORCODE::MTB_VECTORTWONORMERROR)
        std::cerr << "Vector: two norm error.";
    else if (ErrorCode == ERRORCODE::MTB_VECTORCROSSPRODUCTERROR)
        std::cerr << "Vector: cross-product error.";
    else if (ErrorCode == ERRORCODE::MTB_MATRIXADDERROR)
        std::cerr << "Matrix: addition error.";
    else if (ErrorCode == ERRORCODE::MTB_MATRIXSUBTRACTERROR)
        std::cerr << "Matrix: subtraction error.";
    else if (ErrorCode == ERRORCODE::MTB_MATRIXMULTIPLYERROR)
        std::cerr << "Matrix: multiplication error.";
    else if (ErrorCode == ERRORCODE::MTB_MATRIXTRANSPOSEERROR)
        std::cerr << "Matrix: transpose error.";
    else if (ErrorCode == ERRORCODE::MTB_MATRIXMULTVECERROR)
        std::cerr << "Matrix: multiply vector error.";
    else if (ErrorCode == ERRORCODE::MTB_MATRIXLUFACTERROR)
        std::cerr << "Matrix: LU factorization error.";
    else if (ErrorCode == ERRORCODE::MTB_MATRIXLUSOLVEERROR)
        std::cerr << "Matrix: LU solve error.";
    else if (ErrorCode == ERRORCODE::MTB_MATRIXGAUSSELIMINATIONERROR)
        std::cerr << "Matrix: Gaussian elimination error.";
    else if (ErrorCode == ERRORCODE::MTB_MATRIXLDLTFACTERROR)
        std::cerr << "Matrix: LDLT factorization error.";
    else if (ErrorCode == ERRORCODE::MTB_MATRIXLDLTSOLVEERROR)
        std::cerr << "Matrix: LDLT solve error.";

    if (ErrorCode != ERRORCODE::UNSTABLETRUSS && ErrorCode != ERRORCODE::INVALIDCOMMANDLINE)
        std::cerr << '\n' << "Error in input file line : "
                  << m_nLineNumber;

    std::string strM = strMessage.str();
    std::cerr << strM << "\n";
}

