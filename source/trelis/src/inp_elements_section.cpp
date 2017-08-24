/**
 * \file
 * \author Michael Dudzinski
 * \date 10.09.2014
 */

#include "inp_elements_section.hpp"
#include "inp_utils.hpp"

using namespace std;

namespace inp_n {

  const string inp_elements_section::ms_SectionHeader = "********************************** E L E M E N T S ****************************";

  /****************************************************************************/
  inp_elements_section::inp_elements_section()
  :
  m_akSetsTri  (0),
  m_akSetsTetra(0)
  {
    /**/
  }
  /****************************************************************************/
  mxArray*
  inp_elements_section::create_and_fill_matlab_structure()
  {
    const mwSize nSetsQuantity = m_akSetsTri.size()+m_akSetsTetra.size();

    mwSize anDims[2] = { 2, nSetsQuantity };
    mxArray *pkCell  = mxCreateCellArray(2,anDims);

    size_t nCurrSetTri=0,nCurrSetTetra=0;

    for (size_t i=0; i<nSetsQuantity; ++i)
    {
      if (m_akSetsTri.size() && m_akSetsTri[nCurrSetTri].m_nSectionIdx==i)
      {
        /* name of nodeset */
        mxArray *pkNameMatrix = mxCreateString(m_akSetsTri[nCurrSetTri].m_kSectionName.c_str());
        mxSetCell(pkCell,i*2,pkNameMatrix);

        /* content of nodeset */
        mxArray *pkElemsMatrix = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
        
        mxSetData(pkElemsMatrix, m_akSetsTri[nCurrSetTri].m_pfBuffer);
        mxSetM   (pkElemsMatrix, m_akSetsTri[nCurrSetTri].m_nLength );
        mxSetN   (pkElemsMatrix, 4                                  );

        mxSetCell(pkCell,i*2+1,pkElemsMatrix);

        ++nCurrSetTri;
      }

      if (m_akSetsTetra.size() && m_akSetsTetra[nCurrSetTetra].m_nSectionIdx==i)
      {
        /* name of nodeset */
        mxArray *pkNameMatrix = mxCreateString(m_akSetsTetra[nCurrSetTetra].m_kSectionName.c_str());
        mxSetCell(pkCell,i*2,pkNameMatrix);

        /* content of nodeset */
        mxArray *pkElemsMatrix = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
        
        mxSetData(pkElemsMatrix, m_akSetsTetra[nCurrSetTetra].m_pfBuffer);
        mxSetM   (pkElemsMatrix, m_akSetsTetra[nCurrSetTetra].m_nLength );
        mxSetN   (pkElemsMatrix, 5                                      );

        mxSetCell(pkCell,i*2+1,pkElemsMatrix);

        ++nCurrSetTetra;
      }
    }
    
    return pkCell;
  }
  /****************************************************************************/
  void
  inp_elements_section::read_set_header(const string &crkLine      ,
                                        string       &rkSectionName,
                                        string       &rkSectionType)
  const
  {
    size_t nPos    = crkLine.find("TYPE=")+5;
    size_t nLength = crkLine.substr(nPos).find_first_of(WhitespaceString+",");
//    if (nLength!=string::npos) ++nLength;
    
    rkSectionType = crkLine.substr(nPos,nLength);
    trim(rkSectionType);

//    mexPrintf("section type: %s at %d, %d\n",rkSectionType.c_str(),nPos,nLength);

    nPos    = crkLine.find("ELSET=")+6;
    nLength = crkLine.substr(nPos).find_first_of(WhitespaceString+",");

    rkSectionName = crkLine.substr(nPos,nLength);
    trim(rkSectionName);
  }
  /****************************************************************************/
  void
  inp_elements_section::read_element_tri(const string &crkLine,
                                         double       &rfN1   ,
                                         double       &rfN2   ,
                                         double       &rfN3   ,
                                         double       &rfID
                                         )
  const
  {
    size_t nNodeID=0;
    stringstream kStream(crkLine);
    
    char c=0;
      
    kStream>>nNodeID;
    rfID=nNodeID;
      
    do { kStream.get(c); } while ( c!=',');
    kStream>>nNodeID;
    rfN1=nNodeID;
    
    do { kStream.get(c); } while ( c!=',');
    kStream>>nNodeID;
    rfN2=nNodeID;
    
    do { kStream.get(c); } while ( c!=',');
    kStream>>nNodeID;
    rfN3=nNodeID;
  }
  /****************************************************************************/
  void
  inp_elements_section::read_element_tetra(const string &crkLine,
                                           double       &rfN1   ,
                                           double       &rfN2   ,
                                           double       &rfN3   ,
                                           double       &rfN4   ,
                                           double       &rfID    )
  const
  {
    size_t nNodeID=0;
    stringstream kStream(crkLine);

    char c=0;
      
    kStream>>nNodeID;
    rfID=nNodeID;

    do { kStream.get(c); } while ( c!=',');
    kStream>>nNodeID;
    rfN1=nNodeID;

    do { kStream.get(c); } while ( c!=',');
    kStream>>nNodeID;
    rfN2=nNodeID;

    do { kStream.get(c); } while ( c!=',');
    kStream>>nNodeID;
    rfN3=nNodeID;

    do { kStream.get(c); } while ( c!=',');
    kStream>>nNodeID;
    rfN4=nNodeID;
  }
  /****************************************************************************/
  bool
  inp_elements_section::read_structure(ifstream &rkInpFileStream)
  {
//    mexPrintf("Reading structure\n");

    eScanState eState=SS_in_elements_section;
    std::streampos nPrevStreamPos=rkInpFileStream.tellg();

    size_t nCurrSetQuantity=0;
    size_t nCurrSet=0;

    string kLine;
    
    while ((eState!=SS_out_elements_section) && getline_checked(rkInpFileStream,kLine))
    {
//      mexPrintf("line  : %s.\n",kLine.c_str());
      
      /* check for scanner state */
      switch (eState)
      {
        case SS_in_elements_section:
        {
          if (!kLine.find("*ELEMENT"))
          {
            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
//            mexPrintf("Switching from SS_in_elements_section to SS_read_set_header\n");
            break;
          }

          if (!kLine.compare("**"))
          {
            return false; /* signal empty section */
          }

          break;
        }

        case SS_read_set_header:
        {
//          mexPrintf("reading header: %s\n",kLine.c_str());
          nCurrSetQuantity = 0;
          nCurrSet         = m_akSetsTri.size()+m_akSetsTetra.size();

          std::string kSectionName, kSectionType;
          read_set_header(kLine,kSectionName, kSectionType);

          /* dependent on element switch state to read tri or tetra */
          if (!kSectionType.compare("STRI3"))
          {
            m_akSetsTri.resize(m_akSetsTri.size()+1);
            m_akSetsTri.back().m_kSectionName = kSectionName;
            m_akSetsTri.back().m_nSectionIdx  = nCurrSet;

            eState=SS_read_element_tri;
//            mexPrintf("Switching from SS_read_set_header to SS_read_element_tri\n");
            break;
          }
          else if (!kSectionType.compare("C3D4"))
          {
            m_akSetsTetra.resize(m_akSetsTetra.size()+1);
            m_akSetsTetra.back().m_kSectionName = kSectionName;
            m_akSetsTetra.back().m_nSectionIdx  = nCurrSet;
            
            eState=SS_read_element_tetra;
//            mexPrintf("Switching from SS_read_set_header to SS_read_element_tetra\n");
            break;
          }
          else
            throw std::string("inp_elements_section::read_structure: element line must contain exactly 4 or exactly 5 values");

          break;
        }

        case SS_read_element_tri:
        {
//          mexPrintf("reading tri: %s\n",kLine.c_str());

          if (!kLine.find("*ELEMENT"))
          {
            m_akSetsTri.back().m_nLength=nCurrSetQuantity;

            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
//            mexPrintf("Switching from SS_read_element_tri to SS_read_set_header\n");
            break;
          }

          if (!kLine.compare("**"))
          {
            m_akSetsTri.back().m_nLength=nCurrSetQuantity;

            eState = SS_out_elements_section;
//            mexPrintf("Switching from SS_read_element_tri to SS_out_elements_section\n");
            break;
          }

          ++nCurrSetQuantity;

          break;
        }

        case SS_read_element_tetra:
        {
//          mexPrintf("reading tetra: %s\n",kLine.c_str());

          if (!kLine.find("*ELEMENT"))
          {
            m_akSetsTetra.back().m_nLength=nCurrSetQuantity;
            
            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
//            mexPrintf("Switching from SS_read_element_tetra to SS_read_set_header\n");
            break;
          }
          
          if (!kLine.compare("**"))
          {
            m_akSetsTetra.back().m_nLength=nCurrSetQuantity;
            
            eState = SS_out_elements_section;
//            mexPrintf("Switching from SS_read_element_tetra to SS_out_elements_section\n");
            break;
          }
          
          ++nCurrSetQuantity;

          break;
        }
          
        default:
          throw std::string("inp_elements_section::read_structure: eState contains an undefined state");

      } /* switch (eState) */

      nPrevStreamPos=rkInpFileStream.tellg();
    } /* while ((eState!=SS_out_elements_section) && getline_checked(rkInpFileStream,kLine)) */
    
    if (eState != SS_out_elements_section)
      throw std::string("inp_elements_section::read_structure: the end of the ***...* E L E M E N T S *...*** section is missing");

//    for (size_t i=0; i<m_akSetsTri.size(); ++i)
//    {
//        mexPrintf("Elementset tri %s: idx=%d, quantity=%d.\n",
//                  m_akSetsTri[i].m_kSectionName.c_str(),
//                  m_akSetsTri[i].m_nSectionIdx,
//                  m_akSetsTri[i].m_nLength);
//    }
//
//    for (size_t i=0; i<m_akSetsTetra.size(); ++i)
//    {
//      mexPrintf("Elementset tetra %s: idx=%d, quantity=%d.\n",
//                m_akSetsTetra[i].m_kSectionName.c_str(),
//                m_akSetsTetra[i].m_nSectionIdx,
//                m_akSetsTetra[i].m_nLength);
//    }

    return true;
  }
  /****************************************************************************/
  void
  inp_elements_section::read_contents(ifstream &rkInpFileStream)
  {
//    mexPrintf("Reading content\n");

    eScanState eState=SS_in_elements_section;
    std::streampos nPrevStreamPos=rkInpFileStream.tellg();
    
    size_t nCurrSetQuantity=0;
    int    nCurrSetTri=-1;
    int    nCurrSetTetra=-1;

    double n1,n2,n3,n4, ID;
    
    string kLine;

//    getline_checked(rkInpFileStream,kLine);
//    rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
//    mexPrintf("line: %s\n",kLine.c_str());
    
    while ((eState!=SS_out_elements_section) && getline_checked(rkInpFileStream,kLine))
    {
      
      /* check for scanner state */
      switch (eState)
      {
        case SS_in_elements_section:
        {
          if (!kLine.find("*ELEMENT"))
          {
            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
//            mexPrintf("Switching from SS_in_elements_section to SS_read_set_header\n");
            break;
          }

          // already checked in read_structure
//          if (!kLine.compare("**"))
//          {
//            return false; /* signal empty section */
//          }

          break;
        }
          
        case SS_read_set_header:
        {
          nCurrSetQuantity = 0;
//          mexPrintf("reading header 2: %s\n",kLine.c_str());

//          std::string kSectionName;
//          read_set_header(kLine,kSectionName);
          
          
          /* check type of coordinate to read */
//          std::vector<double> afValues;
          
//          nPrevStreamPos=rkInpFileStream.tellg();
//          if (!getline_checked(rkInpFileStream,kLine)) break;
//          read_numeric_values(kLine,afValues);
//          rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
//
////          mexPrintf("line contains %d values.\n",afValues.size());
//          
//          /* dependent on coordinate switch state to read 2D or 3D */
//          switch (afValues.size())
//          {
//            case 4:
//            {
//              ++nCurrSetTri;
//
//              create_vector_array_3d(m_akSetsTri[nCurrSetTri].m_nLength ,
//                                     m_akSetsTri[nCurrSetTri].m_pfBuffer,
//                                     m_akSetsTri[nCurrSetTri].m_pfN1    ,
//                                     m_akSetsTri[nCurrSetTri].m_pfN2    ,
//                                     m_akSetsTri[nCurrSetTri].m_pfN3    );
//
//              eState=SS_read_element_tri;
//              mexPrintf("Switching from SS_read_set_header to SS_read_element_tri\n");
//              break;
//            }
//              
//            case 5:
//            {
//              ++nCurrSetTetra;
//
//              create_vector_array_4d(m_akSetsTetra[nCurrSetTetra].m_nLength ,
//                                     m_akSetsTetra[nCurrSetTetra].m_pfBuffer,
//                                     m_akSetsTetra[nCurrSetTetra].m_pfN1    ,
//                                     m_akSetsTetra[nCurrSetTetra].m_pfN2    ,
//                                     m_akSetsTetra[nCurrSetTetra].m_pfN3    ,
//                                     m_akSetsTetra[nCurrSetTetra].m_pfN4    );
//              
//              eState=SS_read_element_tetra;
//              mexPrintf("Switching from SS_read_set_header to SS_read_element_tetra\n");
//              break;
//            }
//
//            default: /* already checked in read_structure */
//              throw std::string("inp_elements_section::read_structure: element line must contain exactly 4 or exactly 5 values");
//          }

          std::string kSectionName, kSectionType;
          read_set_header(kLine,kSectionName, kSectionType);

          /* dependent on element switch state to read tri or tetra */
          if (!kSectionType.compare("STRI3"))
          {
            ++nCurrSetTri;
            // vector for indices
              
         
              
\
             /*
            create_vector_array_3d(m_akSetsTri[nCurrSetTri].m_nLength ,
                                   m_akSetsTri[nCurrSetTri].m_pfBuffer,
                                   m_akSetsTri[nCurrSetTri].m_pfN1    ,
                                   m_akSetsTri[nCurrSetTri].m_pfN2    ,
                                   m_akSetsTri[nCurrSetTri].m_pfN3    );
              */
              create_vector_array_4d(m_akSetsTri[nCurrSetTri].m_nLength ,
                                     m_akSetsTri[nCurrSetTri].m_pfBuffer,
                                     m_akSetsTri[nCurrSetTri].m_pfId    ,
                                     m_akSetsTri[nCurrSetTri].m_pfN1    ,
                                     m_akSetsTri[nCurrSetTri].m_pfN2    ,
                                     m_akSetsTri[nCurrSetTri].m_pfN3    );
            
            eState=SS_read_element_tri;
//            mexPrintf("Switching from SS_read_set_header to SS_read_element_tri\n");
            break;

          }
          else if (!kSectionType.compare("C3D4"))
          {
            ++nCurrSetTetra;
            
            create_vector_array_5d(m_akSetsTetra[nCurrSetTetra].m_nLength ,
                                   m_akSetsTetra[nCurrSetTetra].m_pfBuffer,
                                   m_akSetsTetra[nCurrSetTetra].m_pfId      ,
                                   m_akSetsTetra[nCurrSetTetra].m_pfN1    ,
                                   m_akSetsTetra[nCurrSetTetra].m_pfN2    ,
                                   m_akSetsTetra[nCurrSetTetra].m_pfN3    ,
                                   m_akSetsTetra[nCurrSetTetra].m_pfN4    );
            
            eState=SS_read_element_tetra;
//            mexPrintf("Switching from SS_read_set_header to SS_read_element_tetra\n");
            break;

          }
          else
            throw std::string("inp_elements_section::read_structure: element line must contain exactly 4 or exactly 5 values");

          break;
        }
          
        case SS_read_element_tri:
        {
          if (!kLine.find("*ELEMENT"))
          {
            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
//            mexPrintf("Switching from SS_read_element_tri to SS_read_set_header\n");
            break;
          }
          
          if (!kLine.compare("**"))
          {
            eState = SS_out_elements_section;
//            mexPrintf("Switching from SS_read_element_tri to SS_out_elements_section\n");
            break;
          }

          read_element_tri(kLine,n1,n2,n3, ID);
          m_akSetsTri[nCurrSetTri].m_pfN1[nCurrSetQuantity]=n1;
          m_akSetsTri[nCurrSetTri].m_pfN2[nCurrSetQuantity]=n2;
          m_akSetsTri[nCurrSetTri].m_pfN3[nCurrSetQuantity]=n3;
            //
          m_akSetsTri[nCurrSetTri].m_pfId[nCurrSetQuantity]=ID;
          ++nCurrSetQuantity;

          break;
        }
          
        case SS_read_element_tetra:
        {
          if (!kLine.find("*ELEMENT"))
          {
            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
//            mexPrintf("Switching from SS_read_element_tetra to SS_read_set_header\n");
            break;
          }
          
          if (!kLine.compare("**"))
          {
            eState = SS_out_elements_section;
//            mexPrintf("Switching from SS_read_element_tetra to SS_out_elements_section\n");
            break;
          }

          read_element_tetra(kLine,n1,n2,n3,n4,ID);
//          mexPrintf("coordinate: %f, %f, %f,\n",x,y,z);
          m_akSetsTetra[nCurrSetTetra].m_pfN1[nCurrSetQuantity]=n1;
          m_akSetsTetra[nCurrSetTetra].m_pfN2[nCurrSetQuantity]=n2;
          m_akSetsTetra[nCurrSetTetra].m_pfN3[nCurrSetQuantity]=n3;
          m_akSetsTetra[nCurrSetTetra].m_pfN4[nCurrSetQuantity]=n4;
          m_akSetsTetra[nCurrSetTetra].m_pfId[nCurrSetQuantity]=ID;
          
          ++nCurrSetQuantity;

          break;
        }
          
        default:
          throw std::string("inp_elements_section::read_structure: eState contains an undefined state");
          
      } /* switch (eState) */
      
      nPrevStreamPos=rkInpFileStream.tellg();
    } /* while ((eState!=SS_out_elements_section) && getline_checked(rkInpFileStream,kLine)) */
  }
  /****************************************************************************/
  mxArray*
  inp_elements_section::scan(ifstream &rkInpFileStream)
  {
    streampos nSectionStart = rkInpFileStream.tellg();

    if (!read_structure(rkInpFileStream)) return mxCreateCellArray(0, NULL);

    rkInpFileStream.clear();
    rkInpFileStream.seekg(nSectionStart);

//    mexPrintf("Reading content\n");

    read_contents(rkInpFileStream);

    return create_and_fill_matlab_structure();
//    return mxCreateCellArray(0, NULL);
  }
  /****************************************************************************/
  
}; /* namespace inp_n */
