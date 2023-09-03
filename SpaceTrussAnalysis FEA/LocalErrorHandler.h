/*********************************************
Copyright(c) 2000-22, S. D. Rajan
All rights reserved

Intermediate Structural Analysis and Design and
Object-Oriented Numerical Analysis via C++

Contains ClocalErrorHandler class definition.

List of improvements made:

*********************************************/
#pragma once

#include <fstream>
#include <string>

class CLocalErrorHandler
{
    public:
        CLocalErrorHandler ();
        ~CLocalErrorHandler ();

        enum class ERRORCODE
        {
            // model data/input errors
            INVALIDNUMNODES, INVALIDNUMELEMENTS, INVALIDDEBUGCODE,
            INVALIDNODENUM, INVALIDNODALFIXITY, INVALIDELEMENTNUM,
            INVALIDCSAREA, INVALIDYM, MISSINGEND,
            UNSTABLETRUSS, INVALIDINPUT, INVALIDCOMMANDLINE,
            CANNOTOPENIFILE, CANNOTOPENOFILE,
            // matrix toolbox errors
            MTB_VECTORADDERROR, MTB_VECTORSUBTRACTERROR, 
            MTB_VECTORDOTPRODUCTTERROR,
            MTB_VECTORNORMALIZEERROR, MTB_VECTORTWONORMERROR,
            MTB_VECTORCROSSPRODUCTERROR,
            MTB_MATRIXADDERROR, MTB_MATRIXSUBTRACTERROR,
            MTB_MATRIXMULTIPLYERROR, MTB_MATRIXTRANSPOSEERROR,
            MTB_MATRIXMULTVECERROR, MTB_MATRIXLUFACTERROR, 
            MTB_MATRIXLUSOLVEERROR, MTB_MATRIXGAUSSELIMINATIONERROR,
            MTB_MATRIXLDLTFACTERROR, MTB_MATRIXLDLTSOLVEERROR
        };

        void ErrorHandler (ERRORCODE& err);
        void ErrorHandler (ERRORCODE& err, int n);

    protected:
        int m_nLineNumber; // line number in the input file
};


