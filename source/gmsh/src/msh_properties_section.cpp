/**
 * \file
 * \author Michael Dudzinski
 * \date 10.09.2014
 */

#include "inp_properties_section.hpp"
#include "inp_utils.hpp"

#include <algorithm>

using namespace std;

namespace inp_n {

  const string inp_properties_section::ms_SectionHeader = "********************************** P R O P E R T I E S ************************";

  /****************************************************************************/
  inp_properties_section::inp_properties_section()
  :
  m_akProperties(0)
  {
    /**/
  }
  /****************************************************************************/
  mxArray*
  inp_properties_section::create_and_fill_matlab_structure()
  {
    const mwSize nSetsQuantity = m_akProperties.size();

    if (nSetsQuantity==0)
      return mxCreateCellArray(0, NULL);

    mwSize anDims[2] = { 2, nSetsQuantity };
    mxArray *pkCell  = mxCreateCellArray(2,anDims);

    for (size_t i=0; i<nSetsQuantity; ++i)
    {
      /* name of property */
      mxArray *pkNameMatrix = mxCreateString(m_akProperties.back().m_kPropertyName.c_str());
      mxSetCell(pkCell,i*2,pkNameMatrix);

      /* content of property */
      mwSize nPropertyDim = 2;
      mxArray *pkPropertyCell  = mxCreateCellArray(1,&nPropertyDim);
      {
        mxArray *pkNameMatrix = mxCreateString(m_akProperties[i].m_kSetName.c_str());
        mxSetCell(pkPropertyCell,0,pkNameMatrix);

        pkNameMatrix = mxCreateString(m_akProperties[i].m_kMaterialName.c_str());
        mxSetCell(pkPropertyCell,1,pkNameMatrix);
      }
      mxSetCell(pkCell,i*2+1,pkPropertyCell);
    }
    
    return pkCell;
  }
  /****************************************************************************/
  void
  inp_properties_section::printState(const eScanState &eState,
                                   const std::string &crkLine)
  const
  {
    switch (eState)
    {
      case SS_in_properties_section:
        mexPrintf("current state: SS_in_properties_section\n");
        break;
        
      case SS_read_solid_section_header:
        mexPrintf("current state: SS_read_solid_section_header\n");
        break;
        
      case SS_read_shell_section_header:
        mexPrintf("current state: SS_read_shell_section_header\n");
        break;
        
      case SS_out_properties_section:
        mexPrintf("current state: SS_out_properties_section\n");
        break;
    }
    mexPrintf("current line: \"%s\"\n",crkLine.c_str());
  }
  /****************************************************************************/
  void
  inp_properties_section::read_solid_section_header(const string &crkLine       ,
                                                    string       &rkSetName     ,
                                                    string       &rkMaterialName)
  const
  {
    size_t nPos    = crkLine.find("ELSET=")+6;
    size_t nLength = crkLine.substr(nPos).find_first_of(WhitespaceString+",");
    
    rkSetName = crkLine.substr(nPos,nLength);
    trim(rkSetName);

    nPos    = crkLine.find("MATERIAL=")+9;
    nLength = crkLine.substr(nPos).find_first_of(WhitespaceString+",");

    rkMaterialName = crkLine.substr(nPos,nLength);
    trim(rkMaterialName);
  }
  /****************************************************************************/
  void
  inp_properties_section::read_shell_section_header(const string &crkLine    ,
                                                    string &rkSetName        ,
                                                    string &rkMaterialName   ,
                                                    string &rkIntegrationName)
  const
  {
    size_t nPos    = crkLine.find("ELSET=")+6;
    size_t nLength = crkLine.substr(nPos).find_first_of(WhitespaceString+",");
    
    rkSetName = crkLine.substr(nPos,nLength);
    trim(rkSetName);

    nPos    = crkLine.find("INTEGRATION=")+12;
    nLength = crkLine.substr(nPos).find_first_of(WhitespaceString+",");
    
    rkIntegrationName = crkLine.substr(nPos,nLength);
    trim(rkIntegrationName);
    
    nPos    = crkLine.find("MATERIAL=")+9;
    nLength = crkLine.substr(nPos).find_first_of(WhitespaceString+",");
    
    rkMaterialName = crkLine.substr(nPos,nLength);
    trim(rkMaterialName);
  }
  /****************************************************************************/
  bool
  inp_properties_section::read_contents(ifstream &rkInpFileStream)
  {
//    mexPrintf("Reading structure\n");

    eScanState eState=SS_in_properties_section;
    std::streampos nPrevStreamPos=rkInpFileStream.tellg();

    size_t nCurrSetQuantity=0;

    string kLine;
    
    while ((eState!=SS_out_properties_section) && getline_checked(rkInpFileStream,kLine))
    {
//      printState(eState,kLine);
      
      /* check for scanner state */
      switch (eState)
      {
        case SS_in_properties_section:
        {
          if (!kLine.find("*SOLID SECTION"))
          {
            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_solid_section_header;
            break;
          }
          else if (!kLine.find("*SHELL SECTION"))
          {
            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_shell_section_header;
            break;
          }
          else if (!kLine.find("**"))
            eState = SS_out_properties_section;
          else
            throw string("inp_elements_section::read_structure: corrupted syntax. I expected *SOLID, *SHELL or **");


          break;
        }

        case SS_read_solid_section_header:
        {
          if (!kLine.find("*SOLID SECTION"))
          {
            m_akProperties.resize(m_akProperties.size()+1);
            
            m_akProperties.back().m_kPropertyName="solid";
            
            read_solid_section_header(kLine                                ,
                                      m_akProperties.back().m_kSetName     ,
                                      m_akProperties.back().m_kMaterialName);
          }
          else if (!kLine.compare("**"))
          {
            eState = SS_out_properties_section;
          }
          else if (!kLine.find("*"))
          {
            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_in_properties_section;
          }
          else
          {
            /* just ignore */
          }

          break;
        }

        case SS_read_shell_section_header:
        {
          if (!kLine.find("*SHELL SECTION"))
          {
            m_akProperties.resize(m_akProperties.size()+1);
            
            m_akProperties.back().m_kPropertyName="shell";
            
            string kIntegrationName;
            read_shell_section_header(kLine                                ,
                                      m_akProperties.back().m_kSetName     ,
                                      m_akProperties.back().m_kMaterialName,
                                      kIntegrationName);
          }
          else if (!kLine.compare("**"))
          {
            eState = SS_out_properties_section;
          }
          else if (!kLine.find("*"))
          {
            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_in_properties_section;
          }
          else
          {
            /* just ignore */
          }
          
          break;
        }
          
        default:
          throw string("inp_elements_section::read_structure: eState contains an undefined state");

      } /* switch (eState) */

      nPrevStreamPos=rkInpFileStream.tellg();
    } /* while ((eState!=SS_out_properties_section) && getline_checked(rkInpFileStream,kLine)) */
    
    if (eState != SS_out_properties_section)
      throw std::string("inp_properties_section::read_structure: the end of the ***...* P R O P E R T I E S *...*** section is missing");

//    for (size_t i=0; i<m_akProperties.size(); ++i)
//    {
//      mexPrintf("Property solid: elset=%s, material=%s\n",
//                m_akProperties[i].m_kSetName     .c_str(),
//                m_akProperties[i].m_kMaterialName.c_str());
//    }

    return true;
  }
  /****************************************************************************/
  mxArray*
  inp_properties_section::scan(ifstream &rkInpFileStream)
  {
    if (!read_contents(rkInpFileStream)) return mxCreateCellArray(0, NULL);

    return create_and_fill_matlab_structure();
//    return mxCreateCellArray(0, NULL);
  }
  /****************************************************************************/
  
}; /* namespace inp_n */
