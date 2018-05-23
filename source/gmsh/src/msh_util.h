#ifndef __MSH_UTIL_H__
#define __MSH_UTIL_H__

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include "mex.h"

namespace msh_n {
  
  const std::string WhitespaceString = " \t\f\v\n\r";

  inline void open_msh_file(const mxArray *cpkMshFileName, std::ifstream &rkMshFileStream)
  {
    // read file name
    char *pcFileName = mxArrayToString(cpkMshFileName);
    if (!pcFileName)
      throw std::string("open_msh_file: mxArrayToString failed");
    
    std::string kFileName(pcFileName);
    mxFree(pcFileName);
    
    // check for *.msh
    if (kFileName.find(".msh") == std::string::npos)
        kFileName.append(".msh");
    
    // open file
    rkMshFileStream.open(kFileName.c_str());
    
    if (!rkMshFileStream.is_open())
      throw std::string("open_msh_file: ifstream::open failed");
  }

  inline bool getline_checked(std::ifstream &rkInpFileStream, std::string &rkLine)
  {
    getline(rkInpFileStream, rkLine);
    
    // check stream state
    if (!rkInpFileStream.good())
    {
      if (rkInpFileStream.eof())
        return false;
      else
        throw std::string("getline_checked: Bad stream state");
    }

    // remove trailing whitespace
    std::size_t nPos = rkLine.find_last_not_of(WhitespaceString);    
    if (nPos != std::string::npos)
      rkLine.erase(nPos+1);
    else
      rkLine.clear();
    
    return true;
  }
  
  inline size_t getSizeNum(const std::string kLine)
  {
    std::stringstream kStream(kLine);
    size_t value;
    kStream >> value;
    if (kStream.fail())
      return 0;
    else
      return value;
  }

  /**
   * \brief Creates a matrix with crnLength rows and 1 column
   * \param[in] crnLength Length of the array
   * \param[out] rpkBuffer Data storage location
   */
  template<typename T>
  void
  create_vector_array_1d(const size_t  &crnLength,
                         T            *&rpkBuffer)
  {
    rpkBuffer = NULL;
    rpkBuffer = (T*) mxCalloc(crnLength, sizeof(T));
    
    if (!rpkBuffer)
      throw std::string("create_vector_array_1d: mxCalloc failed");
  }

  /**
   * \brief Creates a matrix with crnLength rows and 2 columns
   * \param[in] crnLength Length of the array
   * \param[out] rpkBuffer Data storage location
   * \param[out] rpkX column of \f$x\f$-values
   * \param[out] rpkY column of \f$y\f$-values
   */
  template<typename T>
  void
  create_vector_array_2d(const size_t  &crnLength,
                         T            *&rpkBuffer ,
                         T            *&rpkX      ,
                         T            *&rpkY      )
  {
    rpkBuffer = NULL;
    rpkBuffer = (T*) mxCalloc(2*crnLength, sizeof(T));

    if (!rpkBuffer)
      throw std::string("create_vector_array_2d: mxCalloc failed");

    rpkX = rpkBuffer;
    rpkY = rpkBuffer+crnLength;
  }

  /**
   * \brief Creates a matrix with crnLength rows and 3 columns
   * \param[in] crnLength Length of the array
   * \param[out] rpkBuffer Data storage location
   * \param[out] rpkX column of \f$x\f$-values
   * \param[out] rpkY column of \f$y\f$-values
   * \param[out] rpkZ column of \f$z\f$-values
   */
  template<typename T>
  void
  create_vector_array_3d(const size_t  &crnLength,
                         T            *&rpkBuffer ,
                         T            *&rpkX      ,
                         T            *&rpkY      ,
                         T            *&rpkZ      )
  {
    rpkBuffer = NULL;
    rpkBuffer = (T*) mxCalloc(3*crnLength, sizeof(T));
    
    if (!rpkBuffer)
      throw std::string("create_vector_array_3d: mxCalloc failed");
    
    rpkX = rpkBuffer;
    rpkY = rpkBuffer+  crnLength;
    rpkZ = rpkBuffer+2*crnLength;
  }

  /**
   * \brief Creates a matrix with crnLength rows and 4 columns
   * \param[in] crnLength Length of the array
   * \param[out] rpkBuffer Data storage location
   * \param[out] rpkX column of \f$x\f$-values
   * \param[out] rpkY column of \f$y\f$-values
   * \param[out] rpkZ column of \f$z\f$-values
   * \param[out] rpkW column of \f$w\f$-values
   */
  template<typename T>
  void
  create_vector_array_4d(const size_t  &crnLength,
                         T            *&rpkBuffer,
                         T            *&rpkX     ,
                         T            *&rpkY     ,
                         T            *&rpkZ     ,
                         T            *&rpkW     )
  {
    rpkBuffer = NULL;
    rpkBuffer = (T*) mxCalloc(4*crnLength, sizeof(T));
    
    if (!rpkBuffer)
      throw std::string("create_vector_array_4d: mxCalloc failed");
    
    rpkX = rpkBuffer;
    rpkY = rpkBuffer+  crnLength;
    rpkZ = rpkBuffer+2*crnLength;
    rpkW = rpkBuffer+3*crnLength;
  }
  
  /**
   * \brief Creates a matrix with crnLength rows and 5 columns
   * \param[in] crnLength Length of the array
   * \param[out] rpkBuffer Data storage location
   * \param[out] rpkX column of \f$x\f$-values
   * \param[out] rpkY column of \f$y\f$-values
   * \param[out] rpkZ column of \f$z\f$-values
   * \param[out] rpkW column of \f$w\f$-values
   */
  template<typename T>
  void
  create_vector_array_5d(const size_t  &crnLength,
                         T            *&rpkBuffer,
                         T            *&rpkX     ,
                         T            *&rpkY     ,
                         T            *&rpkZ     ,
                         T            *&rpkW     ,
                         T            *&rpkV     )
  {
    rpkBuffer = NULL;
    rpkBuffer = (T*) mxCalloc(5*crnLength, sizeof(T));
        
    if (!rpkBuffer)
      throw std::string("create_vector_array_5d: mxCalloc failed");
        
    rpkX = rpkBuffer;
    rpkY = rpkBuffer+  crnLength;
    rpkZ = rpkBuffer+2*crnLength;
    rpkW = rpkBuffer+3*crnLength;
    rpkV = rpkBuffer+4*crnLength;
  }

  /**
   * \brief Free space occupied by rpkBuffer
   * \param[inout] pkBuffer Pointer to the sapce to be freed
   */
  template<typename T>
  void
  safe_free(T *&rpkBuffer)
  {
    if (rpkBuffer)
    {
      mxFree(rpkBuffer);
      rpkBuffer=NULL;
    }
  }

}; // namespace msh_n

#endif