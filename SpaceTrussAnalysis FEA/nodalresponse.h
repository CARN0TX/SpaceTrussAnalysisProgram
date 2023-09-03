/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-22, S. D. Rajan
All rights reserved

Intermediate Structural Analysis and Design and
Object-Oriented Numerical Analysis via C++

Contains CNodalResponse class definition.

List of improvements made:

*********************************************/
#pragma once

class CNodalResponse
{
    public:
        // ctor and dtor
        CNodalResponse ();
        ~CNodalResponse ();

        // accessor function
        void GetDisplacements (float& fXDisp, float& fYDisp, float& fZDisp) const;

        // modifier function
        void SetDisplacements (const float fXDisp, const float fYDisp, const float fZDisp);

    private:
        float m_fXDisp;	// x-displacement
        float m_fYDisp;	// y-displacement
        float m_fZDisp;	// z-displacement
};
