/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-22, S. D. Rajan
All rights reserved

Intermediate Structural Analysis and Design and
Object-Oriented Numerical Analysis via C++

Contains CElementResponse class definition.

List of improvements made:

*********************************************/
#pragma once

class CElementResponse
{
    public:
        // ctor and dtor
        CElementResponse ();
        ~CElementResponse ();

        // accessor function
        void GetData (float& fStrain, float& fStress, 
                      float& fForce) const;

        // modifier function
        void SetData (const float fStrain, const float fStress,
                      const float fForce);

    private:
        float m_fStrain; // strain
        float m_fStress; // stress
        float m_fForce;	 // force
};
