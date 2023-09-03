/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-22, S. D. Rajan
All rights reserved

Intermediate Structural Analysis and Design and
Object-Oriented Numerical Analysis via C++

Contains CNode class definition.

List of improvements made:

*********************************************/
#pragma once

class CNode
{
    public:
        // ctor and dtor
        CNode ();
        ~CNode ();

        // accessor functions
        void GetCoords (float&, float&, float&) const;
        void GetFixity (int&, int&, int&) const;
        void GetLoads (float&, float&, float&, float&) const;

        // modifier functions
        void SetCoords (const float, const float, const float);
        void SetFixity (const int, const int, const int);
        void SetLoads (const float, const float, const float, const float);

    private:
        float m_fXCoor;		// x-coordinate
        float m_fYCoor;		// y-coordinate
        float m_fZCoor;		// z-coordinate
        int   m_nXFC;		// x-fixity code
        int   m_nYFC;		// y-fixity code
        int   m_nZFC;		// z-fixity code
        float m_fXForce;	// x-force
        float m_fYForce;	// y-force
        float m_fZForce;	// z-force
        float m_fdelT;      // temprature change
};
