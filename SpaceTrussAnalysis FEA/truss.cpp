/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-22, S. D. Rajan
All rights reserved

Intermediate Structural Analysis and Design and
Object-Oriented Numerical Analysis via C++

Contains CTruss class implementation.

List of improvements made:

*********************************************/
#include <iomanip>
#include <sstream>
#include <string>
#include "truss.h"
#include "MatToolBox.h"
#include "parserEXH.h"
#include "constants.h"

/* ==================================================================
   ======================= CTruss class =============================
   ================================================================== */

CTruss::CTruss ()
// ---------------------------------------------------------------------------
// Function: default ctor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nNodes = m_nElements = m_nDOF = 0;
    m_nDebugLevel = m_nLineNumber = 0;
    m_strComment = "**";                 // lines starts with **  
    m_strDelimiters = ", ";              // commas and blank spaces only
}

CTruss::~CTruss ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CTruss::Banner (std::ostream& OF) const
// ---------------------------------------------------------------------------
// Function: prints the program banner on the output stream
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    OF << '\n';
    OF << "\t\t--------------------------------------------" << '\n';
    OF << "\t\t        Planar Truss Analysis Program       " << '\n';
    OF << "\t\tIntroduction to Structural Analysis & Design" << '\n';
    OF << "\t\t          (c) 2000-23, S. D. Rajan          " << '\n';
    OF << "\t\t        Enhanced By: Ravi Chaudhary         " << '\n';
    OF << "\t\t--------------------------------------------" << '\n';
}

void CTruss::PrepareIO (int argc, char* argv[])
// ---------------------------------------------------------------------------
// Function: opens the input and output files
// Input:    command line arguments (currently unused)
// Output:   none
// ---------------------------------------------------------------------------
{
    CFileIO FIO;
    if (argc == 1)
    {
        // open the input file
        FIO.OpenInputFileByName ("Complete input file name: ", m_FileInput,
                                 std::ios::in);

        // open the output file
        FIO.OpenOutputFileByName ("Complete output file name: ", m_FileOutput,
                                  std::ios::out);
    }
    else if (argc == 3) // planartruss input_file output_file
    {
        m_FileInput.open (argv[1], std::ios::in);
        if (!m_FileInput)
            ErrorHandler (CGlobalErrorHandler::ERRORCODE::CANNOTOPENIFILE);
        m_FileOutput.open (argv[2], std::ios::out);
        if (!m_FileOutput)
            ErrorHandler (CGlobalErrorHandler::ERRORCODE::CANNOTOPENOFILE);
        std::cout << "\n";
        std::cout << argv[1] << " opened as input file.\n";
        std::cout << argv[2] << " opened as output file.\n";
    }
	else
    {
        ErrorHandler (CLocalErrorHandler::ERRORCODE::INVALIDCOMMANDLINE);
    }

    // print banner
    Banner (m_FileOutput);
}

void CTruss::ReadProblemSize ()
// ---------------------------------------------------------------------------
// Function: Reads the size of the problem being solved
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    bool bEOF = false;

    // initialize line number counter
    m_nLineNumber = 0;

    // header line (*spacetruss)
    m_Parser.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                        m_nTokens, m_strDelimiters, m_strComment,
                        bEOF);

    // header line (*heading)
    m_Parser.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                        m_nTokens, m_strDelimiters, m_strComment,
                        bEOF);

    // read the problem description
    m_Parser.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                        m_nTokens, m_strDelimiters, m_strComment,
                        bEOF);

    // header line (*control)
    m_Parser.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                        m_nTokens, m_strDelimiters, m_strComment,
                        bEOF);

    // read number of nodes, elements and debug level
    m_Parser.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                        m_nTokens, m_strDelimiters, m_strComment,
                        bEOF);

    m_Parser.GetIntValue (m_strVTokens[0], m_nNodes);
    m_Parser.GetIntValue (m_strVTokens[1], m_nElements);
    m_Parser.GetIntValue (m_strVTokens[2], m_nDebugLevel);

    // check data for validity
    if (m_nNodes <= 1)
        ErrorHandler (CLocalErrorHandler::ERRORCODE::INVALIDNUMNODES);
    if (m_nElements <= 0)
        ErrorHandler (CLocalErrorHandler::ERRORCODE::INVALIDNUMELEMENTS);
    if (m_nDebugLevel < 0 || m_nDebugLevel > 1)
        ErrorHandler (CLocalErrorHandler::ERRORCODE::INVALIDDEBUGCODE);

    // dynamic memory allocations for arrays
    SetSize ();
}

void CTruss::SetSize ()
// ---------------------------------------------------------------------------
// Function: Carries out memory allocation for all the major arrays
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // allocate space for nodal data
    m_NodalData.SetSize (m_nNodes);
    // allocate space for nodal response data
    m_NodalResponseData.SetSize (m_nNodes);
    // allocate space for element data
    m_ElementData.SetSize (m_nElements);
    // allocate space for element response data
    m_ElementResponseData.SetSize (m_nElements);

    // allocate and initialize major matrices
    m_nDOF = DOFPN*m_nNodes;
    m_dSSM.SetSize (m_nDOF, m_nDOF);
    m_dSND.SetSize (m_nDOF, 1);
    m_dSNF.SetSize (m_nDOF, 1);
    m_dSSM.Set (0.0);
    m_dSND.Set (0.0);
    m_dSNF.Set (0.0);
}

void CTruss::ReadTrussModel ()
// ---------------------------------------------------------------------------
// Function: Read the truss model data from the input file
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    bool bEOF = false;
    try
    {
        // read problem size
        ReadProblemSize ();

        // read nodal data
        int nN;
        float fX, fY, fZ;
        m_Parser.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                            m_nTokens, m_strDelimiters, m_strComment,
                            bEOF);

        for (int i = 1; i <= m_nNodes; i++)
        {
            m_Parser.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                                m_nTokens, m_strDelimiters, m_strComment,
                                bEOF);
            // node #
            m_Parser.GetIntValue (m_strVTokens[0], nN);
            // x-coordinate
            m_Parser.GetFloatValue (m_strVTokens[1], fX);
            // y-coordinate
            m_Parser.GetFloatValue (m_strVTokens[2], fY);
            // z-coordinate
            m_Parser.GetFloatValue(m_strVTokens[3], fZ);
            // valid node number?
            if (nN <= 0 || nN > m_nNodes)
                ErrorHandler (CLocalErrorHandler::ERRORCODE::INVALIDNODENUM);
            // store the data
            m_NodalData(nN).SetCoords (fX, fY, fZ);
        }

        // read nodal fixity conditions and dispalcement values
        float fXD, fYD, fZD;
        int nXFC, nYFC, nZFC;
        m_Parser.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                            m_nTokens, m_strDelimiters, m_strComment,
                            bEOF);

        for (int i = 1; i <= m_nNodes; i++)
        {
            m_Parser.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                                m_nTokens, m_strDelimiters, m_strComment,
                                bEOF);
            if (m_nTokens !=7 ) {
                break;
            }
            // node #
            m_Parser.GetIntValue (m_strVTokens[0], nN);
            // x-fixity code
            m_Parser.GetIntValue (m_strVTokens[1], nXFC);
            // y-fixity code
            m_Parser.GetIntValue (m_strVTokens[2], nYFC);
            // z-fixity code
            m_Parser.GetIntValue (m_strVTokens[3], nZFC);
            // x-displacement value
            m_Parser.GetFloatValue (m_strVTokens[4], fXD);
            // y-displacement value
            m_Parser.GetFloatValue (m_strVTokens[5], fYD);
            // z-displacement value
            m_Parser.GetFloatValue (m_strVTokens[6], fZD);
            // valid data?
            if (nN <= 0 || nN > m_nNodes)
                ErrorHandler (CLocalErrorHandler::ERRORCODE::INVALIDNODENUM);
            if (nXFC < 0 || nXFC > 1) 
                ErrorHandler (CLocalErrorHandler::ERRORCODE::INVALIDNODALFIXITY);
            if (nYFC < 0 || nYFC > 1) 
                ErrorHandler (CLocalErrorHandler::ERRORCODE::INVALIDNODALFIXITY);
            if (nZFC < 0 || nZFC > 1)
                ErrorHandler(CLocalErrorHandler::ERRORCODE::INVALIDNODALFIXITY);
            // store the data
            m_NodalData(nN).SetFixity (nXFC, nYFC, nZFC);
            m_NodalResponseData(nN).SetDisplacements(fXD, fYD, fZD);
        }


        // read nodal loads/forces
        float fXF, fYF, fZF, delT;
        m_Parser.GetTokens(m_FileInput, m_nLineNumber, m_strVTokens,
                           m_nTokens, m_strDelimiters, m_strComment,
                           bEOF);

        for (int i = 1; i <= m_nNodes; i++)
        {
            m_Parser.GetTokens(m_FileInput, m_nLineNumber, m_strVTokens,
                               m_nTokens, m_strDelimiters, m_strComment,
                               bEOF);
            if (m_nTokens != 5) {
                break;
            }
            // node #
            m_Parser.GetIntValue(m_strVTokens[0], nN);
            // x-force
            m_Parser.GetFloatValue(m_strVTokens[1], fXF);
            // y-force
            m_Parser.GetFloatValue(m_strVTokens[2], fYF);
            // z-force
            m_Parser.GetFloatValue(m_strVTokens[3], fZF);
            // delT-temprature change
            m_Parser.GetFloatValue(m_strVTokens[4], delT);

            // valid data?
            if (nN <= 0 || nN > m_nNodes)
                ErrorHandler(CLocalErrorHandler::ERRORCODE::INVALIDNODENUM);

            // store the data
            m_NodalData(nN).SetLoads(fXF, fYF, fZF, delT);
        }

        // read element data
        int nE, nSN, nEN;
        float fA, fE, falp;
        /*m_Parser.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                            m_nTokens, m_strDelimiters, m_strComment,
                            bEOF);*/
        for (int i = 1; i <= m_nElements; i++)
        {
            m_Parser.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                                m_nTokens, m_strDelimiters, m_strComment,
                                bEOF);
            // element #
            m_Parser.GetIntValue (m_strVTokens[0], nE);
            // start node #
            m_Parser.GetIntValue (m_strVTokens[1], nSN);
            // end node #
            m_Parser.GetIntValue (m_strVTokens[2], nEN);
            // x/s area
            m_Parser.GetFloatValue (m_strVTokens[3], fA);
            // young's modulus
            m_Parser.GetFloatValue (m_strVTokens[4], fE);
            // cofficient of expansion
            m_Parser.GetFloatValue(m_strVTokens[5], falp);
            // valid data?
            if (nE <= 0 || nE > m_nElements) 
                ErrorHandler (CLocalErrorHandler::ERRORCODE::INVALIDELEMENTNUM);
            if (nSN <= 0 || nSN > m_nNodes) 
                ErrorHandler (CLocalErrorHandler::ERRORCODE::INVALIDNODENUM);
            if (nEN <= 0 || nEN > m_nNodes) 
                ErrorHandler (CLocalErrorHandler::ERRORCODE::INVALIDNODENUM);
            if (fA <= 0.0f) 
                ErrorHandler (CLocalErrorHandler::ERRORCODE::INVALIDCSAREA);
            if (fE <= 0.0f) 
                ErrorHandler (CLocalErrorHandler::ERRORCODE::INVALIDYM);
            // store the data
            m_ElementData(nE).SetData (nSN, nEN, fA, fE, falp);
        }

        // end keyword?
        m_Parser.GetTokens (m_FileInput, m_nLineNumber, m_strVTokens,
                            m_nTokens, m_strDelimiters, m_strComment,
                            bEOF);
        if (m_strVTokens[0].substr(0,4) != "*end")
            ErrorHandler (CLocalErrorHandler::ERRORCODE::MISSINGEND);
    }

    // trap all input file errors here
    catch (CLocalErrorHandler::ERRORCODE &err)
    {
        m_LEH.ErrorHandler (err);
    }

    catch (CGlobalErrorHandler::ERRORCODE &err)
    {
        ErrorHandler (CLocalErrorHandler::ERRORCODE::INVALIDINPUT);
    }

    catch (...)
    {
        std::cout << "Sorry, could not catch the error whatever it is.\n";
    }
}

void CTruss::Analyze ()
// ---------------------------------------------------------------------------
// Function: Implements the FEA steps
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    ////////////////////////////////////////////////////////////////////////////
    //Thermal loading
    for (int i = 1; i <= m_nElements; i++) {
        int nSN, nEN;
        float fA, fE, falp;
        float fX1, fX2, fY1, fY2, fZ1, fZ2, fL;
        float sNfXF, sNfYF, sNfZF, sNfdelT;
        float eNfXF, eNfYF, eNfZF, eNfdelT;
        m_ElementData(i).GetData(nSN, nEN, fA, fE, falp);
        m_NodalData(nSN).GetLoads(sNfXF, sNfYF, sNfZF, sNfdelT);
        m_NodalData(nEN).GetLoads(eNfXF, eNfYF, eNfZF, eNfdelT);
        m_NodalData(nSN).GetCoords(fX1, fY1, fZ1);
        m_NodalData(nEN).GetCoords(fX2, fY2, fZ2);
        fL = sqrt((fX2 - fX1) * (fX2 - fX1) +
            (fY2 - fY1) * (fY2 - fY1) +
            (fZ2 - fZ1) * (fZ2 - fZ1));

        float fdelT = (sNfdelT + eNfdelT) / 2;
        double dl = double((fX2 - fX1) / fL);
        double dm = double((fY2 - fY1) / fL);
        double dn = double((fZ2 - fZ1) / fL);

        CVector<double> df1(2); // element load vector at both nodes
        df1(1) = double(fA * fE * falp * fdelT);
        df1(2) = double(fA * fE * falp * fdelT);

        CVector<double> dTLV(DOFPE);  //thermal load vector
        dTLV(1) = dl * df1(1);        //thermal load in x direction at start node of element i
        dTLV(2) = dm * df1(1);        //thermal load in y direction at start node of element i
        dTLV(3) = dn * df1(1);        //thermal load in z direction at start node of element i
        dTLV(4) = dl * df1(2);        //thermal load in x direction at end node of element i
        dTLV(5) = dm * df1(2);        //thermal load in y direction at end node of element i
        dTLV(6) = dn * df1(2);        //thermal load in z direction at end node of element i

        //adding thermal load value to the load variable for that node
        sNfXF += float(dTLV(1));  //x-force
        sNfYF += float(dTLV(2));  //y-force
        sNfZF += float(dTLV(3));  //z-force
        eNfXF += float(dTLV(4));  //x-force
        eNfYF += float(dTLV(5));  //y-force
        eNfZF += float(dTLV(6));  //z-force

        //assigning thermal load values to the nodes
        m_NodalData(nSN).SetLoads(sNfXF, sNfYF, sNfZF, sNfdelT);
        m_NodalData(nEN).SetLoads(eNfXF, eNfYF, eNfZF, eNfdelT);
    }

    ////////////////////////////////////////////////////////////////////////////

    // construct structural nodal load vector
    float fXF, fYF, fZF, fdelT;

    for (int i = 1; i <= m_nNodes; i++)
    {
        m_NodalData(i).GetLoads (fXF, fYF, fZF, fdelT);
        m_dSNF(DOFPN * i - 2) = fXF;
        m_dSNF(DOFPN * i- 1) = fYF;
        m_dSNF(DOFPN * i) = fZF;
    }

    // construct structural stiffness matrix
	ConstructK ();

	// impose boundary conditions
	ImposeBC ();

	// solve for the nodal displacements
	Solve ();

	// compute element response
	Response ();
	
    // create output file
	CreateOutput ();
}

void CTruss::ConstructK ()
// ---------------------------------------------------------------------------
// Function: Constructs the structural stiffness matrix
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    CVector<int> nE(DOFPE);                      // the nodal degrees-of-freedom vector
    CMatrix<double> dkl(2,2), dkg(DOFPE, DOFPE); // k' and k
    CMatrix<double> dT(2, DOFPE), dTT(DOFPE,2);  // T and T(t)
    CMatrix<double> dTemp(2, DOFPE);             // temporary storage
    float fX1, fX2, fY1, fY2, fZ1, fZ2, fL;      // (x1, y1, z1), (x2, y2, z2), L
    float fA, fE, falp;                          // A, E, alpha
    int nSN, nEN;                                // a, b
    
    // initialize
    dT.Set (0.0);
    
    // loop thro' all elements
    for (int i=1; i <= m_nElements; i++)
    {
        // construct k(6x6) in two steps
        m_ElementData(i).GetData (nSN, nEN, fA, fE, falp);
        m_NodalData(nSN).GetCoords (fX1, fY1, fZ1);
        m_NodalData(nEN).GetCoords (fX2, fY2, fZ2);
        fL = sqrt((fX2-fX1)*(fX2-fX1) +
                  (fY2-fY1)*(fY2-fY1) +
                  (fZ2-fZ1)*(fZ2-fZ1));
        // form AE/L 
        double dAEOL = double(fA*fE/fL);
        if (CGlobalErrorHandler::FloatException(dAEOL))
            throw (CGlobalErrorHandler::FLOATEXCEPTION);
        // construct klocal
        dkl(1,1) = dAEOL;
        dkl(1,2) = -dAEOL;
        dkl(2,1) = dkl(1,2);
        dkl(2,2) = dkl(1,1);
        // local-to-global transformation matrix
        double dl = double((fX2-fX1)/fL);
        double dm = double((fY2-fY1)/fL);
        double dn = double((fZ2-fZ1)/fL);
        dT(1, 1) = dl; dT(1, 2) = dm; dT(1, 3) = dn;
        dT(2, 4) = dl; dT(2, 5) = dm; dT(2, 6) = dn;
        // transpose of the T matrix
        m_MTBDP.Transpose (dT, dTT);
        // form k'*T
        m_MTBDP.Multiply (dkl, dT, dTemp);
        // construct kglobal=T(T)*k'*T
        m_MTBDP.Multiply (dTT, dTemp, dkg);

        // assemble into structural K
        nE(1) = DOFPN * nSN - 2; nE(2) = DOFPN * nSN - 1; nE(3) = DOFPN * nSN;
        nE(4) = DOFPN * nEN - 2; nE(5) = DOFPN * nEN - 1; nE(6) = DOFPN * nEN;
        for (int j=1; j <= DOFPE; j++)
        {
            int nRow = nE(j);
            for (int k=1; k <= DOFPE; k++)
            {
                int nCol = nE(k);
                m_dSSM(nRow, nCol) += dkg(j,k);
            }
        }

        
        

        // debug?
        if (m_nDebugLevel == 1)
        {
            std::ostringstream strPrompt;
            strPrompt << "Stiffness Matrix for element " << i;
            m_MTBDP.PrintMatrixRowWise (dkg, strPrompt.str(), m_FileOutput);
        }
    }

    // debug?
    if (m_nDebugLevel == 1)
    {
        m_MTBDP.PrintMatrixRowWise (m_dSSM, "Structural Stiffness (Before BCs)",
                                    m_FileOutput);
    }
}

void CTruss::ImposeBC ()
// ---------------------------------------------------------------------------
// Function: Imposes the essential boundary conditions
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    int nXFC, nYFC, nZFC;  // global degrees-of-freedom at node i
    
    // loop thro' all nodes
    for (int i=1; i <= m_nNodes; i++)
    {
        m_NodalData(i).GetFixity (nXFC, nYFC, nZFC);
        if (nXFC == 1) // fixed?
        {
            int nGDOF = DOFPN * i - 2;
            SuppressDOF (nGDOF);
        }
        if (nYFC == 1) // fixed?
        {
            int nGDOF = DOFPN * i - 1;
            SuppressDOF (nGDOF);
        }
        if (nZFC == 1) // fixed?
        {
            int nGDOF = DOFPN * i;
            SuppressDOF (nGDOF);
        }
    }

    // debug?
    if (m_nDebugLevel == 1)
    {
        m_MTBDP.PrintMatrixRowWise (m_dSSM, "Structural Stiffness (After BCs)",
                                    m_FileOutput);
        m_MTBDP.PrintVector (m_dSNF, "Structural Nodal Forces (After BCs)",
                             m_FileOutput);
    }
}

void CTruss::SuppressDOF (const int nEqn)
// ---------------------------------------------------------------------------
// Function: Imposes the essential boundary conditions
// Input:    Equation number
// Output:   none
// ---------------------------------------------------------------------------
{
    for (int j=1; j <= m_nDOF; j++)
    {
        // zero out the row
        m_dSSM(nEqn, j) = 0.0;
        // zero out the column
        m_dSSM(j, nEqn) = 0.0;
    }
    // set diagonal to 1
    m_dSSM(nEqn, nEqn) = 1.0;

    // set RHS to zero
    m_dSNF(nEqn) = 0.0;
}

void CTruss::Solve ()
// ---------------------------------------------------------------------------
// Function: Solves the system equations for the nodal displacements
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // implement Gauss Elimination Technique
    double TOL = 1.0e-6;

    /*m_MTBDP.LDLTFactorization(m_dSSM, TOL);
    m_MTBDP.LDLTSolve(m_dSSM, m_dSND, m_dSNF);*/
    m_MTBDP.GaussElimination (m_dSSM, m_dSND, m_dSNF, TOL);
    // ErrorHandler (CLocalErrorHandler::ERRORCODE::UNSTABLETRUSS);
    for (int i=1; i <= m_nNodes; i++)
    {
        float m_fXDisp = static_cast<float>(m_dSND(3*i-2));
        float m_fYDisp = static_cast<float>(m_dSND(3*i-1));
        float m_fZDisp = static_cast<float>(m_dSND(3*i));
        m_NodalResponseData(i).SetDisplacements (m_fXDisp, m_fYDisp, m_fZDisp);
    }
}

void CTruss::Response ()
// ---------------------------------------------------------------------------
// Function: Computes the element response
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    CVector<int> nE(DOFPE);                     // the nodal degrees-of-freedom vector
    CMatrix<double> dT(DOFPN,DOFPE),
                    dTT(DOFPE, DOFPN);          // T and T(t)
    CVector<double> dLD(3), dND(DOFPE);         // nodal displacements d', d
    float fX1, fX2, fY1, fY2, fZ1, fZ2, fL;     // (x1, y1, z1), (x2, y2, z2), L
    float fA, fE, falp;                         // A, E, alpha
    int nSN, nEN;                               // a, b
    
    // initialize
    dT.Set (0.0);
    
    // loop thro' all elements
    for (int i=1; i <= m_nElements; i++)
    {
        // form strain, stress and force in two steps
        m_ElementData(i).GetData (nSN, nEN, fA, fE, falp);
        m_NodalData(nSN).GetCoords (fX1, fY1, fZ1);
        m_NodalData(nEN).GetCoords (fX2, fY2, fZ2);
        fL = sqrt((fX2-fX1)*(fX2-fX1) +
                  (fY2-fY1)*(fY2-fY1) +
                  (fZ2-fZ1)*(fZ2-fZ1));
        // local-to-global transformation matrix
        double dl = double((fX2-fX1)/fL);
        double dm = double((fY2-fY1)/fL);
        double dn = double((fZ2 - fZ1) / fL);
        dT(1, 1) = dl; dT(1, 2) = dm; dT(1, 3) = dn;
        dT(2, 4) = dl; dT(2, 5) = dm; dT(2, 6) = dn;

        // get element nodal displacements
        nE(1) = DOFPN * nSN - 2; nE(2) = DOFPN * nSN - 1; nE(3) = DOFPN * nSN;
        nE(4) = DOFPN * nEN - 2; nE(5) = DOFPN * nEN - 1; nE(6) = DOFPN * nEN;
        for (int j=1; j <= DOFPE; j++)
        {
            dND(j) = m_dSND(nE(j));
        }

        // form d'=T*d
        m_MTBDP.MatMultVec (dT, dND, dLD);

        // strain
        float fStrain = float(dLD(2) - dLD(1))/fL;
        // stress
        float fStress = fE*fStrain;
        // force
        float fForce = fStress*fA;

        // store the element response
        m_ElementResponseData(i).SetData (fStrain,
                                          fStress, fForce);
    }
}

void CTruss::CreateOutput ()
// ---------------------------------------------------------------------------
// Function: Creates the output file containing the results.
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // print the problem size
    m_FileOutput << '\n';
    m_FileOutput << "PROBLEM SIZE" << '\n';
    m_FileOutput << "------------" << '\n';
    m_FileOutput << "   Number of nodes : " << m_nNodes << '\n';
    m_FileOutput << "Number of elements : " << m_nElements << '\n';

    // print the nodal coordinates
    m_FileOutput << '\n';
    float fX, fY, fZ;
    m_FileOutput << "NODAL COORDINATES" << '\n';
    m_FileOutput << "-----------------" << '\n';
    m_FileOutput << "Node" << "   " << "X-Coordinate"
                           << "    " << "Y-Coordinate"
                           << "    " << "Z-Coordinate" << '\n';
    m_FileOutput << "----" << "   " << "------------"
                 << "    " << "------------" 
                 << "    " << "------------" << '\n';
    for (int i = 1; i <= m_nNodes; i++)
    {
        m_NodalData(i).GetCoords(fX, fY, fZ);
        m_FileOutput << std::setw(4) << i << "   "
            << std::setw(12) << fX
            << "    " << std::setw(12) << fY
            << "    " << std::setw(12) << fZ<< '\n';
    }

    // print the nodal fixities
    m_FileOutput << '\n';
    int nXFC, nYFC, nZFC;
    m_FileOutput << "NODAL FIXITIES" << '\n';
    m_FileOutput << "--------------" << '\n';
    m_FileOutput << "Node" << "   " << "X-Fixity"
                           << "    " << "Y-Fixity"
                           << "    " << "Z-Fixity" << '\n';
    m_FileOutput << "----" << "   " << "--------"
                 << "    " << "--------"
                 << "    " << "--------" << '\n';
    for (int i = 1; i <= m_nNodes; i++)
    {
        m_NodalData(i).GetFixity(nXFC, nYFC, nZFC);
        m_FileOutput << std::setw(4) << i << "   "
            << std::setw(8) << nXFC
            << "    " << std::setw(8) << nYFC
            << "    " << std::setw(8) << nZFC<< '\n';
    }

    // print the nodal forces
    m_FileOutput << '\n';
    float fXF, fYF, fZF, fdelT;
    m_FileOutput << "NODAL FORCES" << '\n';
    m_FileOutput << "------------" << '\n';
    m_FileOutput << "Node" << "    " << "X-Force" 
                           << "     " << "Y-Force"
                           << "     " << "Z-Force"
                           << "     " << "Delta T" << '\n';
    m_FileOutput << "----" 
                 << "    " << "-------" 
                 << "     " << "-------"
                 << "     " << "-------"
                 << "     " << "-------" << '\n';
    for (int i = 1; i <= m_nNodes; i++)
    {
        m_NodalData(i).GetLoads(fXF, fYF, fZF, fdelT);
        m_FileOutput << std::setw(4) << i << "   "
            << std::setw(8) << fXF
            << "    " << std::setw(8) << fYF
            << "    " << std::setw(8) << fZF
            << "    " << std::setw(8) << fdelT << '\n';
    }
    m_FileOutput << "\n===================== FE RESULTS ========================" << '\n';
    // print the nodal displacements
    m_FileOutput << '\n';
    float fXDisp, fYDisp, fZDisp;
    m_FileOutput << "NODAL DISPLACEMENTS" << '\n';
    m_FileOutput << "-------------------" << '\n';
    m_FileOutput << "Node" << "   " << "X-Displacement"
                 << "    " << "Y-Displacement"
                 << "    " << "Z-Displacement" << '\n';
    m_FileOutput << "----" << "   " << "--------------"
                 << "    " << "--------------" 
                 << "    " << "--------------" << '\n';
    for (int i=1; i <= m_nNodes; i++)
    {
        m_NodalResponseData(i).GetDisplacements (fXDisp, fYDisp, fZDisp);
        m_FileOutput << std::setw(4) << i << "   "
                     << std::setw(14) << fXDisp 
                     << "    " << std::setw(14) << fYDisp
                     << "    " << std::setw(14) << fZDisp << '\n';
    }
    
    // print the element strain, stress and force
    m_FileOutput << '\n';
    float fStrain, fStress, fForce;
    m_FileOutput << "ELEMENT RESPONSE  (Tension is positive)" << '\n';
    m_FileOutput << "----------------" << '\n';
    m_FileOutput << "Element"
                 << "         " << "Strain"
                 << "         " << "Stress"
                 << "         " << "Force" << '\n';
    m_FileOutput << "-------"
                 << "        " << "-------"
                 << "        " << "-------"
                 << "        " << "------" << '\n';
    for (int i = 1; i <= m_nElements; i++)
    {
        m_ElementResponseData(i).GetData(fStrain, fStress, fForce);
        m_FileOutput << std::setw(7) << i << "   "
            << std::setw(12) << fStrain
            << "    " << std::setw(11) << fStress
            << "    " << std::setw(10) << fForce << '\n';
    }

    // print the nodal reactions
    m_FileOutput << '\n';
    m_FileOutput << "NODAL Reactions" << '\n';
    m_FileOutput << "---------------" << '\n';
    m_FileOutput << "Node" << "   " << "X-Reaction"
        << "    " << "Y-Reaction"
        << "    " << "Z-Reaction" << '\n';
    m_FileOutput << "----" << "   " << "----------"
        << "    " << "----------"
        << "    " << "----------" << '\n';
    for (int i = 1; i <= m_nNodes; i++)
    {
        m_NodalData(i).GetLoads(fXF, fYF, fZF, fdelT);
        m_FileOutput << std::setw(4) << i << "   ";
        m_NodalData(i).GetFixity(nXFC, nYFC, nZFC);
        if (nXFC == 1) // fixed?
        {
            m_FileOutput << std::setw(10) << fXF << "   ";
        }
        if (nYFC == 1) // fixed?
        {
            m_FileOutput << std::setw(11) << fYF << "   ";
        }
        if (nZFC == 1) // fixed?
        {
            m_FileOutput << std::setw(11) << fZF << "   ";
        }
        m_FileOutput << "\n";
    }
}

void CTruss::TerminateProgram ()
// ---------------------------------------------------------------------------
// Function: Closes input and output files
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    // close the input and output files
    m_FileInput.close ();
    m_FileOutput.close ();

    std::cout << "\nExecution completed successfully." << std::endl;
}

void CTruss::ErrorHandler (CLocalErrorHandler::ERRORCODE ErrorCode) 
{
    throw ErrorCode;
}

void CTruss::ErrorHandler (CGlobalErrorHandler::ERRORCODE ErrorCode) const
{
    throw ErrorCode;
}

void CTruss::DisplayErrorMessage (CLocalErrorHandler::ERRORCODE err)
{
    m_LEH.ErrorHandler (err, m_nLineNumber);
}
