/**
 * inp_utils.cpp
 *
 * Copyright: Michael Dudzinski,
 *  Department of the Theory of Electrical Engineering,
 *  University of the Federal Armed Forces Hamburg,
 *  Hamburg,
 *  Germany
 */

/**
 * \file
 * \author Michael Dudzinski
 * \date 10.09.2014
 *
 * Any of the functions defined here will throw an exception string in case of
 * any error.
 */

#include "inp_utils.hpp"

#include <iostream>
#include <fstream>

#include <sstream>
#include <string>
#include <vector>

#include "mex.h"

namespace inp_n {

  /**
   * \brief Trim trailing whitespaces in rkString
   * \param[inout] crkString string to be trimmed
   */
  void
  trimr(std::string &rkString)
  {
    std::size_t nPos = rkString.find_last_not_of(WhitespaceString);
    
    if (nPos!=std::string::npos)
      rkString.erase(nPos+1);
    else
      rkString.clear();
  }
  
  
  /**
   * \brief Trim leading whitespaces in rkString
   * \param[inout] crkString string to be trimmed
   */
  void
  triml(std::string &rkString)
  {
    std::size_t nPos = rkString.find_first_not_of(WhitespaceString);
    
    if (nPos!=std::string::npos)
      rkString.erase(0,nPos);
    else
      rkString.clear();
  }
  
  
  /**
   * \brief Trim leading and trailing whitespaces in rkString
   * \param[inout] crkString string to be trimmed
   */
  void
  trim(std::string &rkString)
  {
    std::stringstream trimmer;
    
    trimmer << rkString;
    rkString.clear();
    trimmer >> rkString;
  }


  /**
   * \brief Get line from stream. Check state of the stream and return a right
   * trimmed line provided stream state is good.
   * \param[inout] rkInpFileStream The file stream
   * \param[out] rkLine The read line
   */
  bool
  getline_checked(std::ifstream &rkInpFileStream, std::string &rkLine)
  {
    getline(rkInpFileStream, rkLine);
    
    /* check stream state */
    if (!rkInpFileStream.good())
    {
      if (rkInpFileStream.eof())
        return false;
      else
        throw std::string("getline_checked: Bad stream state");
    }

    /* remove trailing whitespace */
    trimr(rkLine);
//    /* remove possible carriage return */
//    if (*(rkLine.end()-1)=='\r') kLine.erase(kLine.end()-1,kLine.end());

    return true;
  }

  /**
   * \brief Read all numeric values present in crkLine
   * \param[in] crkLine string to be parsed
   * \param[out] rafValues vector of numeric value read from crkLine
   */
  void
  read_numeric_values(const std::string &crkLine, std::vector<uint32_t> &rkValues)
  {
    std::istringstream kLineStream(crkLine);
    
    unit32_t nValue=0.;
    
    // parse line
    while (!kLineStream.eof())
    {
      // read until numerical value is encountered
      while (true)
      {
        /* try to read a numeric value */
        kLineStream>>nValue;
        
        /* stream is corrupted */
        if (kLineStream.bad())
          throw std::string("Bad stream state.");
        
        if (kLineStream.fail())
        {
          /* last read operation failed */
          kLineStream.clear();
          kLineStream.ignore();
          
          /* either eof ... */
          if (kLineStream.eof()) break;
          
          /* ... or not a numeric value */
          continue;
        }
        
        /* numeric value successfully read from stream */
        rkValues.push_back(nValue);
        break;
      }
    }
  }

  /**
   * \brief Read all numeric values present in crkLine
   * \param[in] crkLine string to be parsed
   * \param[out] rafValues vector of numeric value read from crkLine
   */
  void
  read_numeric_values(const std::string &crkLine, std::vector<double> &rafValues)
  {
    std::istringstream kLineStream(crkLine);
    
    double fValue=0.;
    
    // parse line
    while (!kLineStream.eof())
    {
      // read until numerical value is encountered
      while (true)
      {
        /* try to read a numeric value */
        kLineStream>>fValue;

        /* stream is corrupted */
        if (kLineStream.bad())
          throw std::string("Bad stream state.");

        if (kLineStream.fail())
        {
          /* last read operation failed */
          kLineStream.clear();
          kLineStream.ignore();

          /* either eof ... */
          if (kLineStream.eof()) break;

          /* ... or not a numeric value */
          continue;
        }

        /* numeric value successfully read from stream */
        rafValues.push_back(fValue);
        break;
      }
    }
  }


  /**
   * \brief Open inp file given by cpkInpFileName.
   *
   * It is assumed that cpkInpFileName is of class string. If the name given in
   * cpkInpFileName does not contain the *.inp ending it will automatically be
   * appended.
   *
   * \param[in] crkInpFileName Name of the inp-file to be openend
   * \param[out] rkInpFileStream Stream containing file contents
   */
  void
  open_inp_file(const mxArray *cpkInpFileName,
                std::ifstream &rkInpFileStream)
  {
    /* read file name */
    char *pcFileName = mxArrayToString(cpkInpFileName);
    if (!pcFileName)
      throw std::string("open_inp_file: mxArrayToString failed");
    
    std::string kFileName(pcFileName);
    mxFree(pcFileName);
    
    /* check for *.inp */
    if (kFileName.find(".inp")==std::string::npos)
      kFileName.append(".inp");
    
    /* open file */
    rkInpFileStream.open(kFileName.c_str());
    
    if (!rkInpFileStream.is_open())
      throw std::string("open_inp_file: ifstream::open failed");
  }


  /**
   * \brief Creates a matrix with crnLength rows and 1 columns
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
    
    rpkX = rpkBuffer;
    rpkY = rpkBuffer+crnLength;
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

}; /* namespace inp_n */
