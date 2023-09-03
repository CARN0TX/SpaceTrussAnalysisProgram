/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-22, S. D. Rajan
All rights reserved

Intermediate Structural Analysis and Design and
Object-Oriented Numerical Analysis via C++

Contains CElement class implementation.

List of improvements made:

*********************************************/
#include "Element.h"

CElement::CElement ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nSN = m_nEN = 0;
    m_fArea = m_fE = m_falp = 0.0f;
}

CElement::~CElement ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CElement::GetData (int& nSN, int& nEN, float& fA,
                        float& fE, float& falp) const
// ---------------------------------------------------------------------------
// Function: gets element-related values
// Input:    variables to hold start node#, end node #, x/s area, E
// Output:   start node#, end node #, x/s area, E
// ---------------------------------------------------------------------------
{
    nSN = m_nSN;
    nEN = m_nEN;
    fA = m_fArea;
    fE = m_fE;
    falp = m_falp;
}

void CElement::SetData (const int nSN, const int nEN,
                        const float fA, const float fE, const float falp)
// ---------------------------------------------------------------------------
// Function: sets element-related values
// Input:    variables holding values for start node#, end node #, 
//           x/s area, E
// Output:   modified value for start node#, end node #, x/s area, E
// ---------------------------------------------------------------------------
{
    m_nSN = nSN;
    m_nEN = nEN;
    m_fArea = fA;
    m_fE = fE;
    m_falp = falp;
}
