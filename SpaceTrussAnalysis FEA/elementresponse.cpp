/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-22, S. D. Rajan
All rights reserved

Intermediate Structural Analysis and Design and
Object-Oriented Numerical Analysis via C++

Contains CElementResponse class implementation.

List of improvements made:

*********************************************/
#include "ElementResponse.h"

CElementResponse::CElementResponse ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_fStrain = m_fStress = m_fForce = 0.0f;
}

CElementResponse::~CElementResponse ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CElementResponse::GetData (float& fE,
                                float& fS, float& fF) const
// ---------------------------------------------------------------------------
// Function: gets the element strain, stress and force
// Input:    variables to hold the values
// Output:   element strain, stress and force values
// ---------------------------------------------------------------------------
{
    fE = m_fStrain;
    fS = m_fStress;
    fF = m_fForce;
}

void CElementResponse::SetData (const float fE,
                                const float fS, const float fF)
// ---------------------------------------------------------------------------
// Function: sgets the element strain, stress and force
// Input:    variables holding the values
// Output:   modified element strain, stress and force values
// ---------------------------------------------------------------------------
{
    m_fStrain = fE;
    m_fStress = fS;
    m_fForce = fF;
}