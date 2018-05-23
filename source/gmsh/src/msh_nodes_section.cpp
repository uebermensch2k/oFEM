/**
 * \file
 * \author Michael Dudzinski
 * \date 10.09.2014
 */

#include "msh_nodes_section.hpp"
#include "msh_utils.hpp"

using namespace std;

namespace msh_n {

  const string msh_nodes_section::ms_SectionHeader = "$Nodes";

  /****************************************************************************/
  msh_nodes_section::msh_nodes_section()
  :
  m_akSets2D(0),
  m_akSets3D(0)//,
//  m_nCoordinateSectionStart(0)
  {
    /**/
  }
  /****************************************************************************/
  mxArray*
  msh_nodes_section::create_and_fill_matlab_structure()
  {
    const mwSize nSetsQuantity = m_akSets2D.size()+m_akSets3D.size();

    mwSize anDims[2] = { 2, nSetsQuantity };
    mxArray *pkCell  = mxCreateCellArray(2,anDims);

    size_t nCurrSet2D=0,nCurrSet3D=0;

    for (size_t i=0; i<nSetsQuantity; ++i)
    {
      if (m_akSets2D.size() && m_akSets2D[nCurrSet2D].m_nSectionIdx==i)
      {
        /* name of nodeset */
        mxArray *pkNameMatrix = mxCreateString(m_akSets2D[nCurrSet2D].m_kSectionName.c_str());
        mxSetCell(pkCell,i*2,pkNameMatrix);

        /* content of nodeset */
        mxArray *pkElemsMatrix = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

        mxSetData(pkElemsMatrix, m_akSets2D[nCurrSet2D].m_pfBuffer);
        mxSetM   (pkElemsMatrix, m_akSets2D[nCurrSet2D].m_nLength );
        mxSetN   (pkElemsMatrix, 2                                );

        mxSetCell(pkCell,i*2+1,pkElemsMatrix);

        ++nCurrSet2D;
      }

      if (m_akSets3D.size() && m_akSets3D[nCurrSet3D].m_nSectionIdx==i)
      {
        /* name of nodeset */
        mxArray *pkNameMatrix = mxCreateString(m_akSets3D[nCurrSet3D].m_kSectionName.c_str());
        mxSetCell(pkCell,i*2,pkNameMatrix);

        /* content of nodeset */
        mxArray *pkElemsMatrix = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

        mxSetData(pkElemsMatrix, m_akSets3D[nCurrSet3D].m_pfBuffer);
        mxSetM   (pkElemsMatrix, m_akSets3D[nCurrSet3D].m_nLength );
        mxSetN   (pkElemsMatrix, 3                                );

        mxSetCell(pkCell,i*2+1,pkElemsMatrix);

        ++nCurrSet3D;
      }
    }

    return pkCell;
  }
  /****************************************************************************/
  void
  msh_nodes_section::printState(const eScanState &eState,
                                const std::string &crkLine)
  const
  {
    switch (eState)
    {
      case SS_in_nodes_section:
        mexPrintf("current state: SS_in_nodes_section\n");
        break;

      case SS_read_set_header:
        mexPrintf("current state: SS_read_set_header\n");
        break;

      case SS_read_node_2d:
        mexPrintf("current state: SS_read_node_2d\n");
        break;

      case SS_read_node_3d:
        mexPrintf("current state: SS_read_node_3d\n");
        break;

      case SS_out_nodes_section:
        mexPrintf("current state: SS_out_nodes_section\n");
        break;
    }
    mexPrintf("current line: \"%s\"\n",crkLine.c_str());
  }
  /****************************************************************************/
  void
  msh_nodes_section::read_set_header(const string &crkLine      ,
                                     string       &rkSectionName)
  const
  {
    //size_t nPos    = crkLine.find("NSET=")+5;
//    size_t nLength = crkLine.find_last_not_of(WhitespaceString+",",nPos);
//    if (nLength!=string::npos) ++nLength;
    //size_t nLength = crkLine.substr(nPos).find_first_of(WhitespaceString+",",nPos);

    //rkSectionName = crkLine.substr(nPos,nLength);
    //trim(rkSectionName);
	rkSectionName = "All";
  }
  /****************************************************************************/
  void
  msh_nodes_section::read_node_2d(const string &crkLine,
                                  double       &rfX    ,
                                  double       &rfY    )
  const
  {
//    std::vector<double> afValues;
//    read_numeric_values(crkLine,afValues);
//    afValues.erase(afValues.begin()); /* first value is node id */
//      if (afValues.size()!=2)
//        throw std::string("msh_nodes_section::read_node_2d: exactly 2 values expected on the line");
//
//    rfX=afValues[0];
//    rfY=afValues[1];

    double fValue=0.;
    stringstream kStream(crkLine);

    char c=0;

    do { kStream.get(c); } while ( c!=' ');
    kStream>>fValue;
    rfX=fValue;

    do { kStream.get(c); } while ( c!=' ');
    kStream>>fValue;
    rfY=fValue;
  }
  /****************************************************************************/
  void
  msh_nodes_section::read_node_3d(const string &crkLine,
                                  double       &rfX    ,
                                  double       &rfY    ,
                                  double       &rfZ    )
  const
  {
//    mexPrintf("line: %s\n",crkLine.c_str());
//    vector<double> afValues;
//    read_numeric_values(crkLine,afValues);
//    afValues.erase(afValues.begin()); /* first value is node id */
//    if (afValues.size()!=3)
//      throw std::string("msh_nodes_section::read_node_3d: exactly 3 values expected on the line");
//
//    rfX=afValues[0];
//    rfY=afValues[1];
//    rfZ=afValues[2];

    double fValue=0.;
    stringstream kStream(crkLine);

    char c=0;

    do { kStream.get(c); } while ( c!=' ');
    kStream>>fValue;
    rfX=fValue;

    do { kStream.get(c); } while ( c!=' ');
    kStream>>fValue;
    rfY=fValue;

    do { kStream.get(c); } while ( c!=' ');
    kStream>>fValue;
    rfZ=fValue;

//    mexPrintf("coordinate: %f, %f, %f\n",rfX,rfY,rfZ);
  }
  /****************************************************************************/
  bool
  msh_nodes_section::read_structure(ifstream &rkMshFileStream)
  {
    eScanState eState=SS_read_set_header;
    std::streampos nPrevStreamPos=rkMshFileStream.tellg();

    size_t nCurrSetQuantity=0;
    size_t nCurrSet=0;

    string kLine;
	int i = 0;

    while ((eState!=SS_out_nodes_section) && getline_checked(rkMshFileStream,kLine))
    {
		// if (i<10){
	  	// 	printState(eState,kLine);
		// 	++i;
		// }

      /* check for scanner state */
      switch (eState)
      {
        case SS_in_nodes_section:
        {
          if (!kLine.find("$Nodes"))
          {
            rkMshFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
            break;
          }

          if (!kLine.compare("$EndNodes"))
            return false; /* signal empty section */

          break;
        }

        case SS_read_set_header:
        {
          nCurrSetQuantity = 0;
          nCurrSet         = m_akSets2D.size()+m_akSets3D.size();

          std::string kSectionName;
          read_set_header(kLine,kSectionName);


          /* check type of coordinate to read */
          std::vector<double> afValues;

          nPrevStreamPos=rkMshFileStream.tellg();
          if (!getline_checked(rkMshFileStream,kLine)) break;
//          mexPrintf("line: %s.\n",kLine.c_str());
          read_numeric_values(kLine,afValues);
          afValues.erase(afValues.begin()); /* first value is node id */
          rkMshFileStream.seekg(nPrevStreamPos); /* line to be reread */

//          mexPrintf("line contains %d values.\n",afValues.size());

          /* dependent on coordinate switch state to read 2D or 3D */
          switch (afValues.size())
          {
            case 2:
            {
              m_akSets2D.resize(m_akSets2D.size()+1);
              m_akSets2D.back().m_kSectionName = kSectionName;
              m_akSets2D.back().m_nSectionIdx  = nCurrSet;

              eState=SS_read_node_2d;
              break;
            }

            case 3:
            {
              m_akSets3D.resize(m_akSets3D.size()+1);
              m_akSets3D.back().m_kSectionName = kSectionName;
              m_akSets3D.back().m_nSectionIdx  = nCurrSet;

              eState=SS_read_node_3d;
              break;
            }

            default:
              throw std::string("msh_nodes_section::read_structure: coordinate line must contain exactly 2 or exactly 3 values");
          }

          break;
        }

        case SS_read_node_2d:
        {
          if (!kLine.find("$Nodes"))
          {
            m_akSets2D.back().m_nLength=nCurrSetQuantity;

            rkMshFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
            break;
          }

          if (!kLine.compare("$EndNodes"))
          {
            m_akSets2D.back().m_nLength=nCurrSetQuantity;

            eState = SS_out_nodes_section;
            break;
          }

          ++nCurrSetQuantity;

          break;
        }

        case SS_read_node_3d:
        {
          if (!kLine.find("$Nodes"))
          {
            m_akSets3D.back().m_nLength=nCurrSetQuantity;

            rkMshFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
            break;
          }

          if (!kLine.compare("$EndNodes"))
          {
            m_akSets3D.back().m_nLength=nCurrSetQuantity;

            eState = SS_out_nodes_section;
            break;
          }

          ++nCurrSetQuantity;

          break;
        }

        default:
          throw std::string("msh_nodes_section::read_structure: eState contains an undefined state");

      } /* switch (eState) */

      nPrevStreamPos=rkMshFileStream.tellg();
    } /* while ((eState!=SS_out_nodes_section) && getline_checked(rkMshFileStream,kLine)) */

    if (eState != SS_out_nodes_section)
      throw std::string("msh_nodes_section::read_structure: the end of the ***...* N O D E S *...*** section is missing");

//    for (size_t i=0; i<m_akSets2D.size(); ++i)
//    {
//        mexPrintf("Nodeset 2D %s: idx=%d, quantity=%d.\n",
//                  m_akSets2D[i].m_kSectionName.c_str(),
//                  m_akSets2D[i].m_nSectionIdx,
//                  m_akSets2D[i].m_nLength);
//    }
//
//    for (size_t i=0; i<m_akSets3D.size(); ++i)
//    {
//      mexPrintf("Nodeset 3D %s: idx=%d, quantity=%d.\n",
//                m_akSets3D[i].m_kSectionName.c_str(),
//                m_akSets3D[i].m_nSectionIdx,
//                m_akSets3D[i].m_nLength);
//    }

    return true;
  }
  /****************************************************************************/
  void
  msh_nodes_section::read_contents(ifstream &rkMshFileStream)
  {
    eScanState eState=SS_read_set_header;
    std::streampos nPrevStreamPos=rkMshFileStream.tellg();

    size_t nCurrSetQuantity=0;
    int    nCurrSet2D=-1;
    int    nCurrSet3D=-1;

    double x,y,z;

    string kLine;

    while ((eState!=SS_out_nodes_section) && getline_checked(rkMshFileStream,kLine))
    {
//      printState(eState,kLine);

      /* check for scanner state */
      switch (eState)
      {
        case SS_in_nodes_section:
        {
          if (!kLine.find("$Nodes"))
          {
            rkMshFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
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
//          mexPrintf("line: %s\n",kLine.c_str());

//          std::string kSectionName;
//          read_node_header(kLine,kSectionName);


          /* check type of coordinate to read */
          std::vector<double> afValues;

          nPrevStreamPos=rkMshFileStream.tellg();
          if (!getline_checked(rkMshFileStream,kLine)) break;
          read_numeric_values(kLine,afValues);
          afValues.erase(afValues.begin()); /* first value is node id */
          rkMshFileStream.seekg(nPrevStreamPos); /* line to be reread */

//          mexPrintf("line contains %d values.\n",afValues.size());

          /* dependent on coordinate switch state to read 2D or 3D */
          switch (afValues.size())
          {
            case 2:
            {
              ++nCurrSet2D;

              create_vector_array_2d(m_akSets2D[nCurrSet2D].m_nLength ,
                                     m_akSets2D[nCurrSet2D].m_pfBuffer,
                                     m_akSets2D[nCurrSet2D].m_pfX     ,
                                     m_akSets2D[nCurrSet2D].m_pfY     );

              eState=SS_read_node_2d;
              break;
            }

            case 3:
            {
              ++nCurrSet3D;

              create_vector_array_3d(m_akSets3D[nCurrSet3D].m_nLength ,
                                     m_akSets3D[nCurrSet3D].m_pfBuffer,
                                     m_akSets3D[nCurrSet3D].m_pfX     ,
                                     m_akSets3D[nCurrSet3D].m_pfY     ,
                                     m_akSets3D[nCurrSet3D].m_pfZ     );

              eState=SS_read_node_3d;
              break;
            }

            default: /* already checked in read_structure */
              throw std::string("msh_nodes_section::read_structure: coordinate line must contain exactly 2 or exactly 3 values");
          }

          break;
        }

        case SS_read_node_2d:
        {
          if (!kLine.find("$Nodes"))
          {
            rkMshFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
            break;
          }

          if (!kLine.compare("$EndNodes"))
          {
            eState = SS_out_nodes_section;
            break;
          }

          read_node_2d(kLine,x,y);
          m_akSets2D[nCurrSet2D].m_pfX[nCurrSetQuantity]=x;
          m_akSets2D[nCurrSet2D].m_pfY[nCurrSetQuantity]=y;

          ++nCurrSetQuantity;

          break;
        }

        case SS_read_node_3d:
        {
          if (!kLine.find("$Nodes"))
          {
            rkMshFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
            break;
          }

          if (!kLine.compare("$EndNodes"))
          {
            eState = SS_out_nodes_section;
            break;
          }

          read_node_3d(kLine,x,y,z);
//          mexPrintf("coordinate: %f, %f, %f,\n",x,y,z);
          m_akSets3D[nCurrSet3D].m_pfX[nCurrSetQuantity]=x;
          m_akSets3D[nCurrSet3D].m_pfY[nCurrSetQuantity]=y;
          m_akSets3D[nCurrSet3D].m_pfZ[nCurrSetQuantity]=z;

          ++nCurrSetQuantity;

          break;
        }

        default:
          throw std::string("msh_nodes_section::read_structure: eState contains an undefined state");

      } /* switch (eState) */

      nPrevStreamPos=rkMshFileStream.tellg();
    } /* while ((eState!=SS_out_nodes_section) && getline_checked(rkMshFileStream,kLine)) */
  }
  /****************************************************************************/
  mxArray*
  msh_nodes_section::scan(ifstream &rkMshFileStream)
  {
    streampos nSectionStart = rkMshFileStream.tellg();

    if (!read_structure(rkMshFileStream)) return mxCreateCellArray(0, NULL);

//    return mxCreateCellArray(0, NULL);

    rkMshFileStream.clear();
    rkMshFileStream.seekg(nSectionStart);

    read_contents(rkMshFileStream);

    return create_and_fill_matlab_structure();
  }
  /****************************************************************************/

}; /* namespace msh_n */
