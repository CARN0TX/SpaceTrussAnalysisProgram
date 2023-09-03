/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-22, S. D. Rajan
All rights reserved

Intermediate Structural Analysis and Design and
Object-Oriented Numerical Analysis via C++

Contains CElement class definition.

List of improvements made:

*********************************************/
#pragma once

class CElement
{
    public:
        // ctor and dtor
        CElement ();
        ~CElement ();

        // accessor function
        void GetData (int& nSN, int& nEN, float& fArea, 
                      float& fE, float& falp) const;

        // modifier functions
        void SetData (const int nSN, const int nEN, const float fArea,
                      const float fE, const float falp);

    private:
        int m_nSN;		// start node
        int m_nEN;		// end node
        float m_fArea;	// x/s area
        float m_fE;		// modulus of elasticity
        float m_falp;   // cofficient of thermal expansion
};
