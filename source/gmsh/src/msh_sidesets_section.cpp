/**
 * \file
 * \author Michael Dudzinski
 * \date 10.09.2014
 */

#include "inp_sidesets_section.hpp"
#include "inp_utils.hpp"

#include <algorithm>
#include <sstream>

using namespace std;

namespace inp_n {

  const string inp_sidesets_section::ms_SectionHeader = "********************************** S I D E S E T S **********************************";

  /****************************************************************************/
  inp_sidesets_section::inp_sidesets_section()
  :
  m_akSurfaces(0)
  {
    /**/
  }
  /****************************************************************************/
  mxArray*
  inp_sidesets_section::create_and_fill_matlab_structure()
  {
    const mwSize nSetsQuantity = m_akSurfaces.size();

    /* first row: surface name, second row: surface cell structure */
    mwSize anDims[2] = { 2, nSetsQuantity };
    mxArray *pkCell  = mxCreateCellArray(2,anDims);
//    return pkCell;

    for (size_t i=0; i<nSetsQuantity; ++i)
    {
      /* first row, name of surface */
      mxArray *pkNameMatrix = mxCreateString(m_akSurfaces[i].m_kSectionName.c_str());
      mxSetCell(pkCell,i*2,pkNameMatrix);

      /* second row, surface cell structure */

      /* number of ELSETs at most expected */
      mwSize nNumELSET  = 0;

      /* we assume that all ELSETs are of same type */
      const sideset::type_e eType = m_akSurfaces[i].m_akSidesets[0].m_eType;
      /* char corresponding to type of ELSETs expected */
      string kELSETType = "i";

      switch (eType)
      {
        case sideset::edge:
        {
          /* an edge has three ELSETs */
          nNumELSET = 3;
          /* an edge has type E */
          kELSETType = "E";

          break;
        }

        case sideset::surface:
        {
          /* a surface has four ELSETs */
          nNumELSET = 4;
          /* a surface has type S */
          kELSETType = "S";

          break;
        }
      }

      /* first row: ELSET name, second row: matrix containing element ids */
      mwSize anSurfDims[2] = { 2, nNumELSET };
      mxArray *pkSurfCell  = mxCreateCellArray(2,anSurfDims);
      {
        /* uniqueness of types */
        bool *abELSET = new bool[nNumELSET];
        for (size_t j=0; j<nNumELSET; ++j) abELSET[j]=false;

        for (size_t j=0; j<m_akSurfaces[i].m_akSidesets.size(); ++j)
        {
          size_t nCurrID = m_akSurfaces[i].m_akSidesets[j].m_nID-1;
//          mexPrintf("nCurrID  : %d\n",nCurrID);
          
          if (abELSET[nCurrID])
            throw string("inp_sidesets_section::create_and_fill_matlab_structure: multiple "+kELSETType+to_string(nCurrID)+" element sets.");
          
          abELSET[nCurrID]=true;
          
          {
//			mexPrintf("ELSET  : name=%s\n", m_akSurfaces[i].m_akSidesets[j].m_kSectionName.c_str());
            /* first row, name of element set */
            mxArray *pkNameMatrix = mxCreateString(m_akSurfaces[i].m_akSidesets[j].m_kSectionName.c_str());
            mxSetCell(pkSurfCell,2*nCurrID,pkNameMatrix);
            
            /* second row, matrix of element ids */
//			mexPrintf("ELSET  : size=%d\n", m_akSurfaces[i].m_akSidesets[j].m_nLength);
            mxArray *pkElemsMatrix = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
            
            mxSetData(pkElemsMatrix, m_akSurfaces[i].m_akSidesets[j].m_pfBuffer);
            mxSetM   (pkElemsMatrix, m_akSurfaces[i].m_akSidesets[j].m_nLength );
            mxSetN   (pkElemsMatrix, 1                                         );
            
            mxSetCell(pkSurfCell,2*nCurrID+1,pkElemsMatrix);
          }
        }
        delete abELSET;
        abELSET = NULL;
      }

      /* add surface cell */
      mxSetCell(pkCell,i*2+1,pkSurfCell);
    }
//    
//      
//      
//      
//
//      bool S1=false, S2=false, S3=false, S4=false;
//      for (size_t j=0; j<m_akSurfaces[i].m_akSidesets.size(); ++j)
//      {
//        int nCellIdx = m_akSurfaces[i].m_akSidesets[j].m_kSectionType[1]-'1'+1;
//
//        switch (nCellIdx)
//        {
//          case 1:
//          {
//            if (S1)
//              throw string("inp_sidesets_section::create_and_fill_matlab_structure: multiple S1 element sets encountered.");
//
//            S1=true;
//
//            break;
//          }
//
//          case 2:
//          {
//            if (S2)
//              throw string("inp_sidesets_section::create_and_fill_matlab_structure: multiple S2 element sets encountered.");
//            
//            S2=true;
//            
//            break;
//          }
//
//          case 3:
//          {
//            if (S3)
//              throw string("inp_sidesets_section::create_and_fill_matlab_structure: multiple S3 element sets encountered.");
//            
//            S3=true;
//            
//            break;
//          }
//
//          case 4:
//          {
//            if (S4)
//              throw string("inp_sidesets_section::create_and_fill_matlab_structure: multiple S4 element sets encountered.");
//            
//            S4=true;
//            
//            break;
//          }
//
//          default:
//            throw string("inp_sidesets_section::create_and_fill_matlab_structure: unknown element sets type found.");
//        }
//
//        mwSize anElementDims[2] = { 2, 1 };
//        mxArray *pkElementSetCell  = mxCreateCellArray(2,anElementDims);
//        {
//          /* name of element set */
//          mxArray *pkNameMatrix = mxCreateString(m_akSurfaces[i].m_akSidesets[j].m_kSectionName.c_str());
//          mxSetCell(pkElementSetCell,0,pkNameMatrix);
//
//          mxArray *pkElemsMatrix = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
//        
//          mxSetData(pkElemsMatrix, m_akSurfaces[i].m_akSidesets[j].m_pfBuffer);
//          mxSetM   (pkElemsMatrix, m_akSurfaces[i].m_akSidesets[j].m_nLength );
//          mxSetN   (pkElemsMatrix, 1                                         );
//        
//          mxSetCell(pkElementSetCell,1,pkElemsMatrix);
//        }
//
//        mxSetCell(pkCell,i*5+nCellIdx,pkElementSetCell);
//      }
//    }
    
    return pkCell;
  }
  /****************************************************************************/
  void
  inp_sidesets_section::read_elementset_header(const string &crkLine      ,
                                               string       &rkSectionName)
  const
  {
    size_t nPos    = crkLine.find("ELSET=")+6;
    size_t nLength = crkLine.substr(nPos).find_first_of(WhitespaceString+",");
    
    rkSectionName = crkLine.substr(nPos,nLength);
    trim(rkSectionName);
  }
  /****************************************************************************/
  void
  inp_sidesets_section::read_element_ids(const string &crkLine     ,
                                         double       *pfElementIDs,
                                         size_t       &rnLastIdx   )
  {
    vector<size_t> anElementIDs;
    read_numeric_values(crkLine,anElementIDs);

    copy(anElementIDs.begin(),anElementIDs.end(),pfElementIDs+rnLastIdx);
    rnLastIdx+=anElementIDs.size();
  }
  /****************************************************************************/
  size_t
  inp_sidesets_section::count_element_ids(const std::string &crkLine)
  {
    vector<size_t> anElementIDs;
    read_numeric_values(crkLine,anElementIDs);

    return anElementIDs.size();
  }
  /****************************************************************************/
  void
  inp_sidesets_section::read_surface_header(const string &crkLine,
                                            string       &rkName )
  const
  {
    size_t nPos    = crkLine.find("NAME=")+5;
    size_t nLength = crkLine.substr(nPos).find_first_of(WhitespaceString+",");
    
    rkName = crkLine.substr(nPos,nLength);
    trim(rkName);
  }
  /****************************************************************************/
  void
  inp_sidesets_section::read_elementset_name_type_and_id(const string    &crkLine,
                                                         string          &rkName ,
                                                         sideset::type_e &reType ,
                                                         size_t          &rkID   )
  const
  {
    rkName = crkLine.substr(0,crkLine.find_first_of(WhitespaceString+","));
    trim(rkName);
//    mexPrintf("read_elementset_name_type_and_id: name=%s\n", rkName.c_str());

    std::string kTypeStr = crkLine.substr(crkLine.find_first_of(WhitespaceString+",")+1);
    trim(kTypeStr);
//    mexPrintf("read_elementset_name_type_and_id: type=%n\n", kTypeStr.c_str());

    switch (kTypeStr[0])
    {
      case 'C': reType = sideset::curve  ; break;
      case 'S': reType = sideset::surface; break;
      case 'E': reType = sideset::edge   ; break;
      case 'F': reType = sideset::face   ; break;
      case 'T': reType = sideset::tri    ; break;
      default:
        throw string("inp_elements_section::read_elementset_name_type_and_id: file format corrupted. Sideset of invalid type found in SURFACE section");
    }

    std::stringstream stream(kTypeStr.substr(1));
    stream>>rkID;
  }
  /****************************************************************************/
  bool
  inp_sidesets_section::read_structure(ifstream &rkInpFileStream)
  {
//    mexPrintf("Reading structure\n");

    eScanState eState=SS_in_sidesets_section;
    std::streampos nPrevStreamPos=rkInpFileStream.tellg();

    size_t nCurrSetQuantity=0;

    string kLine;
    
    while ((eState!=SS_out_sidesets_section) && getline_checked(rkInpFileStream,kLine))
    {
//      mexPrintf("read structure line  : %s.\n",kLine.c_str());
      
      /* check for scanner state */
      switch (eState)
      {
        case SS_in_sidesets_section:
        {
//		  mexPrintf("SIDESETS: entering SS_in_sidesets_section\n");
          if (!kLine.find("*ELSET"))
          {
            /* just entered sidesets section, thus a new surface found */
            m_akSurfaces.resize(1);

            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_elementset_header;
            break;
          }
          if (!kLine.compare("**"))
            return false; /* signal empty section */

          break;
        }

        case SS_read_elementset_header:
        {
//		  mexPrintf("SIDESETS: entering SS_read_elementset_header\n");
          nCurrSetQuantity = 0;

          /*
           * we came here since a *ELSET at the beginning of the line was
           * encountered => new sideset
           */
          m_akSurfaces.back().m_akSidesets.resize(m_akSurfaces.back().m_akSidesets.size()+1);
          read_elementset_header(kLine,m_akSurfaces.back().m_akSidesets.back().m_kSectionName);

          eState=SS_read_element_ids;

          break;
        }

        case SS_read_element_ids:
        {
//			mexPrintf("SIDESETS: entering SS_read_element_ids\n");
//          mexPrintf("reading tri: %s\n",kLine.c_str());

          if (!kLine.find("*ELSET"))
          {
            m_akSurfaces.back().m_akSidesets.back().m_nLength=nCurrSetQuantity;

            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_elementset_header;
            break;
          }
          if (!kLine.find("*SURFACE"))
          {
            m_akSurfaces.back().m_akSidesets.back().m_nLength=nCurrSetQuantity;

            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_surface_header;
            break;
          }
          if (!kLine.compare("**"))
            throw string("inp_elements_section::read_structure: file format corrupted. Sideset section must be terminated by a surface section");


          nCurrSetQuantity+=count_element_ids(kLine);

          break;
        }


        case SS_read_surface_header:
        {
//		  mexPrintf("SIDESETS: entering SS_read_surface_header\n");
          read_surface_header(kLine,m_akSurfaces.back().m_kSectionName);

          eState = SS_read_elementset_name_and_type;

          break;
        }

        case SS_read_elementset_name_and_type:
        {
//		  mexPrintf("SIDESETS: entering SS_read_elementset_name_and_type\n");
          if (!kLine.find("*ELSET"))
          {
            if (m_akSurfaces.back().m_akSidesets.size()==0)
              throw string("inp_elements_section::read_structure: file format corrupted. SURFACE with empty element sets found");
              
            /* new surface */
            m_akSurfaces.resize(m_akSurfaces.size()+1);

            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_elementset_header;
            break;
          }
          if (!kLine.compare("**"))
          {
            if (m_akSurfaces.back().m_akSidesets.size()==0)
              throw string("inp_elements_section::read_structure: file format corrupted. SURFACE with empty element sets found");

            eState = SS_out_sidesets_section;
            break;
          }

          string kName;
          sideset::type_e eType;
          size_t nID;
		  //mexPrintf("SIDESETS: name=%s, type=%d, ID=%d\n", "bla", 1, 2);
          read_elementset_name_type_and_id(kLine, kName, eType, nID);
//		  mexPrintf("SIDESETS: name=%s, type=%d, ID=%d\n",kName.c_str(),eType,nID);

          bool bFound=false;
          for (size_t i=0; i<m_akSurfaces.back().m_akSidesets.size(); ++i)
          {
            if (!m_akSurfaces.back().m_akSidesets[i].m_kSectionName.compare(kName))
            {
              m_akSurfaces.back().m_akSidesets[i].m_eType = eType;
              m_akSurfaces.back().m_akSidesets[i].m_nID   = nID  ;
              bFound=true;
            }
          }

          if (!bFound)
            throw string("inp_elements_section::read_structure: file format corrupted. Element set named in SURFACE section not found");
          

          break;
        }
          
        default:
          throw string("inp_elements_section::read_structure: eState contains an undefined state");

      } /* switch (eState) */

      nPrevStreamPos=rkInpFileStream.tellg();
    } /* while ((eState!=SS_out_sidesets_section) && getline_checked(rkInpFileStream,kLine)) */
    
    if (eState != SS_out_sidesets_section)
      throw std::string("inp_sidesets_section::read_structure: the end of the ***...* S I D E S E T S *...*** section is missing");

//    for (size_t i=0; i<m_akSurfaces.size(); ++i)
//    {
//      mexPrintf("Surface %s\n",m_akSurfaces[i].m_kSectionName.c_str());
//      for (size_t j=0; j<m_akSurfaces[i].m_akSidesets.size(); ++j)
//      {
//        mexPrintf("\tSideset %s: quantity=%d, type=%s.\n",
//                  m_akSurfaces[i].m_akSidesets[j].m_kSectionName.c_str(),
//                  m_akSurfaces[i].m_akSidesets[j].m_nLength             ,
//                  m_akSurfaces[i].m_akSidesets[j].m_kSectionType.c_str());
//      }
//    }

    return true;
  }
  /****************************************************************************/
  void
  inp_sidesets_section::read_contents(ifstream &rkInpFileStream)
  {
    eScanState eState=SS_in_sidesets_section;
    std::streampos nPrevStreamPos=rkInpFileStream.tellg();
    
    size_t nCurrSetQuantity=0;
    int nCurrSet=-1;
    size_t nCurrSurface=0;
    
    string kLine;
    
    while ((eState!=SS_out_sidesets_section) && getline_checked(rkInpFileStream,kLine))
    {
//      mexPrintf("read contents line  : %s.\n",kLine.c_str());
      
      /* check for scanner state */
      switch (eState)
      {
        case SS_in_sidesets_section:
        {
          if (!kLine.find("*ELSET"))
          {
            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_elementset_header;

            break;
          }
//          if (!kLine.compare("**"))
//            return false; /* signal empty section */
          
          
          break;
        }
          
        case SS_read_elementset_header:
        {
          //          mexPrintf("reading header: %s\n",kLine.c_str());
          nCurrSetQuantity = 0;
          ++nCurrSet;

          create_vector_array_1d(m_akSurfaces[nCurrSurface].m_akSidesets[nCurrSet].m_nLength ,
                                 m_akSurfaces[nCurrSurface].m_akSidesets[nCurrSet].m_pfBuffer);
          
          eState=SS_read_element_ids;
          //          mexPrintf("Switching from SS_read_set_header to SS_read_node_ids\n");
          
          break;
        }
          
        case SS_read_element_ids:
        {
          //          mexPrintf("reading tri: %s\n",kLine.c_str());
          
          if (!kLine.find("*ELSET"))
          {
            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_elementset_header;
            break;
          }
          if (!kLine.find("*SURFACE"))
          {
            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_surface_header;
            break;
          }

          read_element_ids(kLine,
                           m_akSurfaces[nCurrSurface].m_akSidesets[nCurrSet].m_pfBuffer,
                           nCurrSetQuantity);
          
          break;
        }
          
          
        case SS_read_surface_header:
        {
          eState = SS_read_elementset_name_and_type;
          
          break;
        }
          
        case SS_read_elementset_name_and_type:
        {
          if (!kLine.find("*ELSET"))
          {
            /* new surface */
            ++nCurrSurface;
            nCurrSet=-1;
            
            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_elementset_header;
            break;
          }
          if (!kLine.compare("**"))
          {
            eState = SS_out_sidesets_section;
            break;
          }
          
          break;
        }
          
        default:
          throw string("inp_elements_section::read_structure: eState contains an undefined state");
          
      } /* switch (eState) */
      
      nPrevStreamPos=rkInpFileStream.tellg();
    } /* while ((eState!=SS_out_sidesets_section) && getline_checked(rkInpFileStream,kLine)) */
  }
  /****************************************************************************/
  mxArray*
  inp_sidesets_section::scan(ifstream &rkInpFileStream)
  {
//	mexPrintf("Entering scan\n");
    streampos nSectionStart = rkInpFileStream.tellg();

    if (!read_structure(rkInpFileStream)) return mxCreateCellArray(0, NULL);

    rkInpFileStream.clear();
    rkInpFileStream.seekg(nSectionStart);

    read_contents(rkInpFileStream);

    return create_and_fill_matlab_structure();
//    return mxCreateCellArray(0, NULL);
  }
  /****************************************************************************/
  
}; /* namespace inp_n */
