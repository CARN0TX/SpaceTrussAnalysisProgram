/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-22, S. D. Rajan
All rights reserved

Intermediate Structural Analysis and Design and
Object-Oriented Numerical Analysis via C++

Contains CNode class implementation.

List of improvements made:

*********************************************/
#include "Node.h"

CNode::CNode ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_fXCoor = m_fYCoor = m_fZCoor = m_fXForce = m_fYForce = m_fZForce = m_fdelT = 0.0f;
    m_nXFC = m_nYFC = m_nZFC = 0;
}

CNode::~CNode ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CNode::GetCoords (float& fX, float& fY, float& fZ) const
// ---------------------------------------------------------------------------
// Function: gets nodal coordinates
// Input:    variables to hold x-coordinate, y-coordinate values
// Output:   x-coordinate, y-coordinate values
// ---------------------------------------------------------------------------
{
    fX = m_fXCoor;
    fY = m_fYCoor;
    fZ = m_fZCoor;
}

void CNode::GetFixity (int& nXFC, int& nYFC, int& nZFC) const
// ---------------------------------------------------------------------------
// Function: gets nodal fixity conditions
// Input:    variables to hold x-fixity, y-fixity values
// Output:   x-fixity, y-fixity values
// ---------------------------------------------------------------------------
{
    nXFC = m_nXFC;
    nYFC = m_nYFC;
    nZFC = m_nZFC;
}

//void CNode::GetCoords(float&, float&, float&) const
//{
//}

void CNode::GetLoads (float& fXF, float& fYF, float& fZF, float& fdelT) const
// ---------------------------------------------------------------------------
// Function: gets nodal loads
// Input:    variables to hold x-force, y-force values
// Output:   x-force, y-force values
// ---------------------------------------------------------------------------
{
    fXF = m_fXForce;
    fYF = m_fYForce;
    fZF = m_fZForce;
    fdelT = m_fdelT;
}

void CNode::SetCoords (const float fX, const float fY, const float fZ)
// ---------------------------------------------------------------------------
// Function: sets nodal coordinates
// Input:    variables holding x-coordinate, y-coordinate values
// Output:   modified x-coordinate, y-coordinate values
// ---------------------------------------------------------------------------
{
    m_fXCoor = fX;
    m_fYCoor = fY;
    m_fZCoor = fZ;
}

void CNode::SetFixity (const int nXFC, const int nYFC, const int nZFC)
// ---------------------------------------------------------------------------
// Function: sets nodal fixities
// Input:    variables x-fixity, y-fixity values
// Output:   modified x-fixity, y-fixity values
// ---------------------------------------------------------------------------
{
    m_nXFC = nXFC;
    m_nYFC = nYFC;
    m_nZFC = nZFC;
}

void CNode::SetLoads (const float fXF, const float fYF, const float fZF, const float fdelT)
// ---------------------------------------------------------------------------
// Function: sets nodal loads
// Input:    variables holding x-force, y-force values
// Output:   modified x-force, y-force values
// ---------------------------------------------------------------------------
{
    m_fXForce = fXF;
    m_fYForce = fYF;
    m_fZForce = fZF;
    m_fdelT = fdelT;
}

