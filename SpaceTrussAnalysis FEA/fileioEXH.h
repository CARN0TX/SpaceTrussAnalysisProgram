/*********************************************
Utility Library Function
Copyright(c) 2000-22, S. D. Rajan
All rights reserved

(non-class) Functions to handle simple file operations

Object-Oriented Numerical Analysis
*********************************************/
#pragma once

#include <iostream>
#include <fstream>	
#include <string>
#include <vector>
#include "getinteractiveEXH.h"

class CFileIO
{
    public:
        void OpenInputFileByName (const std::string& strPrompt, 
                                  std::ifstream& IFile, 
                                  const std::ios::openmode&);
        void OpenOutputFileByName (const std::string& strPrompt,
                                   std::ofstream& OFile,
                                   const std::ios::openmode&);
        void Rewind (std::ifstream& IOFile);
        bool OpenInputFile (std::ifstream& IFile,
                            const std::string& strFileName);
        bool OpenOutputFile (std::ofstream& OFile,
                             const std::string& strFileName);
    private:
        CGetInteractive m_GI;
};
//#endif
