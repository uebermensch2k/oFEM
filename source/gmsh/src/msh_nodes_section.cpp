/**
 * \file
 * \author Harald Scharf
 * \date 13.12.2017
 */

#include "msh_nodes_section.h"
#include "msh_util.h"

using namespace std;

namespace msh_n {

  const string msh_nodes_section::ms_DimHeader = "$Dimension";
  const string msh_nodes_section::ms_NodesHeader = "$Nodes";
  const string msh_nodes_section::ms_NodesEnd = "$EndNodes";

  /****************************************************************************/
  msh_nodes_section::msh_nodes_section() :
    m_akSets2D(0),
    m_akSets3D(0)
  {
    /**/
  }
    
  /****************************************************************************/
  mxArray* msh_nodes_section::create_and_fill_matlab_structure()
  {
    const mwSize nSetsQuantity = m_akSets2D.size() + m_akSets3D.size();

    mwSize anDims[2] = { 2, nSetsQuantity };
    mxArray *pkCell  = mxCreateCellArray(2, anDims);

    size_t nCurrSet2D = 0, nCurrSet3D = 0;

    for (size_t i = 0; i < nSetsQuantity; ++i)
    {
      if (m_akSets2D.size() && m_akSets2D[nCurrSet2D].m_nSectionIdx == i)
      {
        // name of nodeset
        mxArray *pkNameMatrix = mxCreateString(m_akSets2D[nCurrSet2D].m_kSectionName.c_str());
        mxSetCell(pkCell, i*2, pkNameMatrix);

        // content of nodeset
        mxArray *pkElemsMatrix = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
        
        mxSetData(pkElemsMatrix, m_akSets2D[nCurrSet2D].m_pfBuffer);
        mxSetM   (pkElemsMatrix, m_akSets2D[nCurrSet2D].m_nLength );
        mxSetN   (pkElemsMatrix, 2                                );

        mxSetCell(pkCell, i*2+1, pkElemsMatrix);

        ++nCurrSet2D;
      }

      if (m_akSets3D.size() && m_akSets3D[nCurrSet3D].m_nSectionIdx == i)
      {
        // name of nodeset
        mxArray *pkNameMatrix = mxCreateString(m_akSets3D[nCurrSet3D].m_kSectionName.c_str());
        mxSetCell(pkCell, i*2, pkNameMatrix);

        // content of nodeset
        mxArray *pkElemsMatrix = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
        
        mxSetData(pkElemsMatrix, m_akSets3D[nCurrSet3D].m_pfBuffer);
        mxSetM   (pkElemsMatrix, m_akSets3D[nCurrSet3D].m_nLength );
        mxSetN   (pkElemsMatrix, 3                                );

        mxSetCell(pkCell, i*2+1, pkElemsMatrix);

        ++nCurrSet3D;
      }
    }
    return pkCell;
  }  
  
  /****************************************************************************/
  void msh_nodes_section::read_node_2d(const string &crkLine,
                                       double &rfX, double &rfY) const
  {
    double fValue;
    std::vector<double> vfValues;
    stringstream kStream(crkLine);
    
    while (!kStream.eof()) {
      kStream >> fValue;
      vfValues.push_back(fValue);
    }
    if (kStream.fail() || vfValues.size() < 3)
      throw std::string("msh_nodes_section::read_node_2d: at least 1+2 values expected on the line");
    
    // discard first value (node id)
    rfX = vfValues[1];
    rfY = vfValues[2];
  }

  /****************************************************************************/
  void msh_nodes_section::read_node_3d(const string &crkLine,
                                       double &rfX, double &rfY, double &rfZ) const
  {
    double fValue;
    std::vector<double> vfValues;
    stringstream kStream(crkLine);
    
    while (!kStream.eof()) {
      kStream >> fValue;
      vfValues.push_back(fValue);
    }
    if (kStream.fail() || vfValues.size() != 4)
      throw std::string("msh_nodes_section::read_node_3d: exactly 1+3 values expected on the line");
    
    // discard first value (node id)
    rfX = vfValues[1];
    rfY = vfValues[2];    
    rfZ = vfValues[3];    
  }
  
  /****************************************************************************/
  void msh_nodes_section::read_contents(ifstream &rkInpFileStream)
  {
    eScanState eState = SS_read_header_dim;
    size_t nCurrSetQuantity = 0;
    int nCurrSet2D = -1;
    int nCurrSet3D = -1;
    double x,y,z;
    string kLine;
    
    while ((eState != SS_out_nodes_section) && getline_checked(rkInpFileStream, kLine))
    {
      switch (eState)
      {          
        case SS_read_header_dim:
          if (!kLine.compare(ms_DimHeader))
            eState = SS_read_nb_dim;
          break;
        
        case SS_read_nb_dim:
          nGeometryDim = getSizeNum(kLine);
          if (nGeometryDim < 2 || nGeometryDim > 3)
            throw std::string("msh_nodes_section::read_contents: dimension must be 2 or 3");
          eState = SS_read_header_nodes;
          break;
          
        case SS_read_header_nodes:
          if (!kLine.compare(ms_NodesHeader))
            eState = SS_read_nb_nodes;
          break;
          
        case SS_read_nb_nodes:
          if (nGeometryDim == 2) {
            ++nCurrSet2D;
            m_akSets2D.resize(m_akSets2D.size()+1);
            m_akSets2D.back().m_kSectionName = "all_nodes";
            m_akSets2D.back().m_nSectionIdx  = 0;
            m_akSets2D.back().m_nLength = getSizeNum(kLine);
            create_vector_array_2d(m_akSets2D[nCurrSet2D].m_nLength,
                                   m_akSets2D[nCurrSet2D].m_pfBuffer,
                                   m_akSets2D[nCurrSet2D].m_pfX,
                                   m_akSets2D[nCurrSet2D].m_pfY);
          } else {
            ++nCurrSet3D;
            m_akSets3D.resize(m_akSets3D.size()+1);
            m_akSets3D.back().m_kSectionName = "all_nodes";
            m_akSets3D.back().m_nSectionIdx  = 0;
            m_akSets3D.back().m_nLength = getSizeNum(kLine);
            create_vector_array_3d(m_akSets3D[nCurrSet3D].m_nLength,
                                   m_akSets3D[nCurrSet3D].m_pfBuffer,
                                   m_akSets3D[nCurrSet3D].m_pfX,
                                   m_akSets3D[nCurrSet3D].m_pfY,
                                   m_akSets3D[nCurrSet3D].m_pfZ );
          }
          eState = SS_read_node;
          break;
                  
        case SS_read_node:
          if (!kLine.compare(ms_NodesEnd))
          {
            eState = SS_out_nodes_section;
            break;
          }
          if (nGeometryDim == 2) {
            read_node_2d(kLine, x, y);
            m_akSets2D[nCurrSet2D].m_pfX[nCurrSetQuantity]=x;
            m_akSets2D[nCurrSet2D].m_pfY[nCurrSetQuantity]=y;
          } else {
            read_node_3d(kLine, x, y, z);
            m_akSets3D[nCurrSet3D].m_pfX[nCurrSetQuantity]=x;
            m_akSets3D[nCurrSet3D].m_pfY[nCurrSetQuantity]=y;
            m_akSets3D[nCurrSet3D].m_pfZ[nCurrSetQuantity]=z;
          }
          ++nCurrSetQuantity;
          break;
          
        default:          
          throw std::string("msh_nodes_section::read_contents: eState contains an undefined state");        
      } // switch (eState)      
    } // while ((eState != SS_out_nodes_section) && getline_checked(rkInpFileStream, kLine))
  
    if (eState != SS_out_nodes_section) {
      if (eState == SS_read_header_dim)
        throw std::string("msh_nodes_section::read_contents: header of DIMENSION section not found!");
      else
        throw std::string("msh_nodes_section::read_contents: header/end of NODES section not found!");
    }
  }
  
  /****************************************************************************/
  mxArray* msh_nodes_section::scan(ifstream &rkMshFileStream)
  {
    read_contents(rkMshFileStream);
    return create_and_fill_matlab_structure();
  }
  
}; // namespace msh_n  