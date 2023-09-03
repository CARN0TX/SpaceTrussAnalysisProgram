/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-22, S. D. Rajan
All rights reserved

Intermediate Structural Analysis and Design and
Object-Oriented Numerical Analysis via C++

Contains CTruss class definition.

List of improvements made:

*********************************************/
#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "constants.h"
#include "..\LibraryEXH\arraycontainersEXH.h"
#include "..\LibraryEXH\parserEXH.h"
#include "..\LibraryEXH\GlobalErrorHandler.h"
#include "..\LibraryEXH\fileioEXH.h"
#include "node.h"
#include "element.h"
#include "nodalresponse.h"
#include "elementresponse.h"
#include "MatToolBox.h"
#include "LocalErrorHandler.h"

class CTruss 
{
    public:
        // ctor and dtor
        CTruss ();
        ~CTruss ();
        void Banner (std::ostream& OF) const;
        void PrepareIO (int argc, char* argv[]);
        void ReadTrussModel ();
        void Analyze ();
        void TerminateProgram ();
        void DisplayErrorMessage (CLocalErrorHandler::ERRORCODE);

    private:
        int m_nNodes;		// number of nodes
        int m_nElements;	// number of elements
        int m_nDOF;			// total degrees-of-freedom
        int m_nDebugLevel;	// debugging level

        CParser m_Parser;             // to parse input lines
        int     m_nLineNumber;	      // current line number in input file
        int     m_nTokens;            // number of tokens read
        std::vector<std::string> 
                m_strVTokens;         // stores the tokens read
        std::string m_strComment;     // comment characters in the input file
        std::string m_strDelimiters;  // delimiters in the input file
        CLocalErrorHandler m_LEH;     // for handling errors detected by program

        // data storage for 
        CVector<CNode> m_NodalData;                     // nodal data
        CVector<CElement> m_ElementData;                // element data
        CVector<CNodalResponse> m_NodalResponseData;    // nodal response
        CVector<CElementResponse> m_ElementResponseData;// element response

        std::ifstream m_FileInput;	// File Input
        std::ofstream m_FileOutput;	// File Output

        CMatrix<double>     m_dSSM;	 // structural stiffness matrix
        CVector<double>     m_dSND;	 // structural nodal displacements
        CVector<double>     m_dSNF;	 // structural nodal forces
        CMatToolBox<double> m_MTBDP; // double precision toolbox

        void ConstructK ();
        void CreateOutput ();
        void ImposeBC ();
        void Response ();
        void ReadProblemSize ();
        void SetSize ();
        void Solve ();
        void SuppressDOF (const int);

        void ErrorHandler (CLocalErrorHandler::ERRORCODE);        // gateway to local error handler
        void ErrorHandler (CGlobalErrorHandler::ERRORCODE) const; // gateway to global error handler
};
