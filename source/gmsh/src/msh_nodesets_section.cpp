/**
 * \file
 * \author Michael Dudzinski
 * \date 10.09.2014
 */

#include "inp_nodesets_section.hpp"
#include "inp_utils.hpp"

#include <algorithm>

using namespace std;

namespace inp_n {

  const string inp_nodesets_section::ms_SectionHeader = "********************************** N O D E S E T S **********************************";

  /****************************************************************************/
  inp_nodesets_section::inp_nodesets_section()
  :
  m_akNodesets(0)
  {
    /**/
  }
  /****************************************************************************/
  mxArray*
  inp_nodesets_section::create_and_fill_matlab_structure()
  {
    const mwSize nSetsQuantity = m_akNodesets.size();

    if (nSetsQuantity==0)
      return mxCreateCellArray(0, NULL);

//    mexPrintf("Number of nodesets: %u\n",nSetsQuantity);

    mwSize anDims[2] = { 2, nSetsQuantity };
    mxArray *pkCell  = mxCreateCellArray(2,anDims);

    for (size_t i=0; i<nSetsQuantity; ++i)
    {
      /* name of nodeset */
      mxArray *pkNameMatrix = mxCreateString(m_akNodesets[i].m_kSectionName.c_str());
      mxSetCell(pkCell,i*2,pkNameMatrix);

//      mexPrintf("Number of id in nodeset(%u): %u\n",i,m_akNodesets[i].m_kNodeIDs.size());
//
//      mexPrintf("Values:\n");
//      for (size_t j=0; j<m_akNodesets[i].m_kNodeIDs.size(); ++j)
//      {
//        mexPrintf("\t%d\n",m_akNodesets[i].m_kNodeIDs[j]);
//      }

      size_t M = m_akNodesets[i].m_kNodeIDs.size();

      /* create buffer and copy content */
      uint32_t *pkData = NULL;
      create_vector_array_1d(M,pkData);
      copy(m_akNodesets[i].m_kNodeIDs.begin(),m_akNodesets[i].m_kNodeIDs.end(),pkData);
        
      mxArray *pkElemsMatrix = mxCreateNumericMatrix(M,1,mxUINT32_CLASS,mxREAL);
      mxSetData(pkElemsMatrix, pkData);

      mxSetCell(pkCell,i*2+1,pkElemsMatrix);
    }
    
    return pkCell;
  }
  /****************************************************************************/
  void
  inp_nodesets_section::printState(const eScanState &eState,
                                   const std::string &crkLine)
  const
  {
    switch (eState)
    {
      case SS_in_nodesets_section:
        mexPrintf("current state: SS_in_nodesets_section\n");
        break;
        
      case SS_read_set_header:
        mexPrintf("current state: SS_read_set_header\n");
        break;
        
      case SS_read_node_ids:
        mexPrintf("current state: SS_read_node_ids\n");
        break;
        
      case SS_out_nodesets_section:
        mexPrintf("current state: SS_out_nodesets_section\n");
        break;
    }
    mexPrintf("current line: \"%s\"\n",crkLine.c_str());
  }
  /****************************************************************************/
  void
  inp_nodesets_section::read_set_header(const string &crkLine      ,
                                        string       &rkSectionName)
  const
  {
    size_t nPos    = crkLine.find("NSET=")+5;
    size_t nLength = crkLine.substr(nPos).find_first_of(WhitespaceString+",");
    
    rkSectionName = crkLine.substr(nPos,nLength);
    trim(rkSectionName);
  }
  /****************************************************************************/
  void
  inp_nodesets_section::read_node_ids(const string     &crkLine ,
                                      vector<uint32_t> &rkNodeIDs)
  {
    read_numeric_values(crkLine,rkNodeIDs);
  }
  /****************************************************************************/
  void
  inp_nodesets_section::read_contents(ifstream &rkInpFileStream)
  {
    eScanState eState=SS_in_nodesets_section;
    std::streampos nPrevStreamPos=rkInpFileStream.tellg();
    
    string kLine;
    
    while ((eState!=SS_out_nodesets_section) && getline_checked(rkInpFileStream,kLine))
    {
//      printState(eState,kLine);
      
      /* check for scanner state */
      switch (eState)
      {
        case SS_in_nodesets_section:
        {
          if (!kLine.find("*NSET"))
          {
            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
          }
          else if (!kLine.find("**"))
          {
            eState = SS_out_nodesets_section;
          }
          
          /* already checked above */
//          if (!kLine.compare("**"))
//          {
//            return false; /* signal empty section */
//          }
          
          break;
        }
          
        case SS_read_set_header:
        {
          if (!kLine.find("*NSET"))
          {
            /* new nodeset found */
            std::string kSectionName;
            read_set_header(kLine,kSectionName);
            
            m_akNodesets.resize(m_akNodesets.size()+1);
            m_akNodesets.back().m_kSectionName = kSectionName;

            eState=SS_read_node_ids;

          }
          else
          {
            throw std::string("inp_nodesets_section::read_content: NSET header expected");
          }
          
          break;
        }
          
        case SS_read_node_ids:
        {
          if (!kLine.find("*NSET"))
          {
            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
            
            break;
          }
          else if (!kLine.find("**"))
          {
            eState = SS_out_nodesets_section;
          }
          else if (!kLine.find("*"))
          {
            throw std::string("inp_nodesets_section::read_content: syntax corrupted. I expected either numerical values or *NSET");
          }

          read_node_ids(kLine,m_akNodesets.back().m_kNodeIDs);

          break;
        }
          
        default:
          throw std::string("inp_nodesets_section::read_content: eState contains an undefined state");
          
      } /* switch (eState) */
      
      nPrevStreamPos=rkInpFileStream.tellg();
    } /* while ((eState!=SS_out_nodesets_section) && getline_checked(rkInpFileStream,kLine)) */
  }
  /****************************************************************************/
  mxArray*
  inp_nodesets_section::scan(ifstream &rkInpFileStream)
  {
    streampos nSectionStart = rkInpFileStream.tellg();

    read_contents(rkInpFileStream);

    return create_and_fill_matlab_structure();
  }
  /****************************************************************************/
  
}; /* namespace inp_n */
