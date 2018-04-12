/**
 * \file
 * \author Harald Scharf
 * \date 14.12.2017
 */

#include "msh_elements_section.h"

using namespace std;

namespace msh_n {

  const string msh_elements_section::ms_DimHeader = "$Dimension";
  const string msh_elements_section::ms_PhysicalEntitiesHeader = "$PhysicalNames";
  const string msh_elements_section::ms_PhysicalEntitiesEnd = "$EndPhysicalNames";
  const string msh_elements_section::ms_ElementsHeader = "$Elements";
  const string msh_elements_section::ms_ElementsEnd = "$EndElements";

  /****************************************************************************/
  msh_elements_section::msh_elements_section()
  :
  m_akSetsLine(0),
  m_akSetsTri(0),
  m_akSetsTetra(0),
  m_akSurfaces(0)
  {
    /**/
  }

  /****************************************************************************/
  int msh_elements_section::findIndex(const std::vector<size_t> &v, const size_t nId)
  {
    for (size_t i = 0; i < v.size(); ++i)
    {
      if (v[i] == nId)
        return i;
    }
    return -1;
  }

  /****************************************************************************/
  void msh_elements_section::findElementByLine(const size_t nLineN1, const size_t nLineN2,
                                               std::vector<double> &rvnE1,
                                               std::vector<double> &rvnE2,
                                               std::vector<double> &rvnE3)
  {
    size_t N1, N2, N3, nId;
    bool bFound1, bFound2, bFound3;
    std::unordered_set<size_t> elSet = {nLineN1, nLineN2};
    for (size_t nCurrSet = 0; nCurrSet < m_akSetsTri.size(); ++nCurrSet)
    {
      for (size_t i = 0; i < m_akSetsTri[nCurrSet].m_nLength; ++i)
      {
        N1 = m_akSetsTri[nCurrSet].m_pfN1[i];
        N2 = m_akSetsTri[nCurrSet].m_pfN2[i];
        N3 = m_akSetsTri[nCurrSet].m_pfN3[i];
        bFound1 = (elSet.find(N1) != elSet.end());
        bFound2 = (elSet.find(N2) != elSet.end());
        bFound3 = (elSet.find(N3) != elSet.end());
        if ( (bFound1 && (bFound2 || bFound3)) || (bFound2 && bFound3) ) {
          nId = (size_t)m_akSetsTri[nCurrSet].m_pfId[i];
          if (bFound1) {
            if (bFound2)
              rvnE1.push_back((double)nId);
            else
              rvnE3.push_back((double)nId);
          } else
            rvnE2.push_back((double)nId);
        }
      }
    }
  }
  
  /****************************************************************************/
  mxArray* msh_elements_section::create_and_fill_matlab_structure()
  {
    mwSize nElSetsQuantity = 0;
    if (nGeometryDim == 2) {
      nElSetsQuantity = m_akSetsTri.size();
    }
    if (nGeometryDim == 3) {
      nElSetsQuantity = m_akSetsTetra.size();
    }
    if (nElSetsQuantity == 0)
      throw std::string("msh_elements_section::create_and_fill_matlab_structure: number of elements is zero");
    
    mwSize anDims[2] = { 2, nElSetsQuantity };
    mxArray *pkCell  = mxCreateCellArray(2, anDims);

    size_t nCurrSetTri = 0, nCurrSetTetra = 0;

    for (size_t i = 0; i < nElSetsQuantity; ++i)
    {
      if (nGeometryDim == 2 && m_akSetsTri[nCurrSetTri].m_nSectionIdx == i)
      {
        // name of nodeset
        mxArray *pkNameMatrix = mxCreateString(m_akSetsTri[nCurrSetTri].m_kSectionName.c_str());
        mxSetCell(pkCell, i*2, pkNameMatrix);

        // content of nodeset
        mxArray *pkElemsMatrix = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
        
        mxSetData(pkElemsMatrix, m_akSetsTri[nCurrSetTri].m_pfBuffer);
        mxSetM   (pkElemsMatrix, m_akSetsTri[nCurrSetTri].m_nLength );
        mxSetN   (pkElemsMatrix, 4                                  );

        mxSetCell(pkCell, i*2+1, pkElemsMatrix);

        ++nCurrSetTri;
      }

      if (nGeometryDim == 3 && m_akSetsTetra[nCurrSetTetra].m_nSectionIdx == i)
      {
        // name of nodeset
        mxArray *pkNameMatrix = mxCreateString(m_akSetsTetra[nCurrSetTetra].m_kSectionName.c_str());
        mxSetCell(pkCell, i*2, pkNameMatrix);

        // content of nodeset
        mxArray *pkElemsMatrix = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
        
        mxSetData(pkElemsMatrix, m_akSetsTetra[nCurrSetTetra].m_pfBuffer);
        mxSetM   (pkElemsMatrix, m_akSetsTetra[nCurrSetTetra].m_nLength );
        mxSetN   (pkElemsMatrix, 5                                      );

        mxSetCell(pkCell, i*2+1, pkElemsMatrix);

        ++nCurrSetTetra;
      }
    }
    
    return pkCell;
  }

  /****************************************************************************/
  mxArray* msh_elements_section::create_and_fill_matlab_structure_bd()
  {    
    const mwSize nSetsQuantity = m_akSurfaces.size();

    // first row: surface name, second row: surface cell structure
    mwSize anDims[2] = { 2, nSetsQuantity };
    mxArray *pkCell  = mxCreateCellArray(2, anDims);

    for (size_t nCurrSet = 0; nCurrSet < nSetsQuantity; ++nCurrSet)
    {
      // first row, name of surface
      mxArray *pkNameMatrix = mxCreateString(m_akSurfaces[nCurrSet].m_kSectionName.c_str());
      mxSetCell(pkCell, nCurrSet*2, pkNameMatrix);

      // second row, surface cell structure
      
      // number of ELSETs at most expected
      mwSize nNumELSET  = 0;

      // we assume that all ELSETs are of same type
      const sideset::type_e eType = m_akSurfaces[nCurrSet].m_akSidesets[0].m_eType;
      switch (eType)
      {
        case sideset::edge:
          nNumELSET = 3;  // an edge has three ELSETs
          break;

        case sideset::surface:
          nNumELSET = 4; // a surface has four ELSETs
          break;
        
        default:
          throw "msh_elements_section::create_and_fill_matlab_structure_bd: invalid eType";
      }

      // first row: ELSET name, second row: matrix containing element ids
      mwSize anSurfDims[2] = { 2, nNumELSET };
      mxArray *pkSurfCell  = mxCreateCellArray(2, anSurfDims);
      
      if (nGeometryDim == 2 && m_akSurfaces[nCurrSet].m_akSidesets.size() != nNumELSET)
        throw "msh_elements_section::create_and_fill_matlab_structure_bd: wrong number of Sidesets";

      if (nGeometryDim == 3 && m_akSurfaces[nCurrSet].m_akSidesets.size() != 1)
        throw "msh_elements_section::create_and_fill_matlab_structure_bd: wrong number of Sidesets";
      
      if (nGeometryDim == 2)
      {
        for (size_t i = 0; i < nNumELSET; ++i)
        {
          // first row, name of element set
          //mexPrintf("ELSET  : name=%s\n", m_akSurfaces[nCurrSet].m_akSidesets[i].m_kSectionName.c_str());
          mxArray *pkElsetNameMatrix = mxCreateString(m_akSurfaces[nCurrSet].m_akSidesets[i].m_kSectionName.c_str());
          mxSetCell(pkSurfCell, 2*i, pkElsetNameMatrix);
            
          // second row, matrix of element ids
          //mexPrintf("ELSET  : size=%d\n", m_akSurfaces[nCurrSet].m_akSidesets[i].m_nLength);
          mxArray *pkElemsMatrix = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);         
          mxSetData(pkElemsMatrix, m_akSurfaces[nCurrSet].m_akSidesets[i].m_pfBuffer);
          mxSetM   (pkElemsMatrix, m_akSurfaces[nCurrSet].m_akSidesets[i].m_nLength );
          mxSetN   (pkElemsMatrix, 1                                         );
          mxSetCell(pkSurfCell, 2*i+1, pkElemsMatrix);
        }
      } else {  // surface ELSET has dimension 1 (list of element IDs)
        // first row, name of element set
        mxArray *pkElsetNameMatrix = mxCreateString(m_akSurfaces[nCurrSet].m_akSidesets[0].m_kSectionName.c_str());
        mxSetCell(pkSurfCell, 0, pkElsetNameMatrix);

        // second row, matrix of element ids
        mxArray *pkElemsMatrix = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);         
        mxSetData(pkElemsMatrix, m_akSurfaces[nCurrSet].m_akSidesets[0].m_pfBuffer);
        mxSetM   (pkElemsMatrix, m_akSurfaces[nCurrSet].m_akSidesets[0].m_nLength );
        mxSetN   (pkElemsMatrix, 1                                         );
        mxSetCell(pkSurfCell, 1, pkElemsMatrix);
      }

      // add surface cell
      mxSetCell(pkCell, nCurrSet*2+1, pkSurfCell);
    }
    return pkCell;
  }
  
  /****************************************************************************/
  void msh_elements_section::read_entity(const string &crkLine, size_t &rnId,
                                         size_t &rnDim, std::string &rkName) const
  {
    stringstream kStream(crkLine);
    char c = 0;
    rkName.clear();
    
    kStream >> rnDim;
    kStream >> rnId;
    do { kStream.get(c); } while (c != '"');
    do
    {
      kStream.get(c);
      if (c != '"')
        rkName += c;
    } while (c != '"');
    if (kStream.fail())
      throw std::string("msh_elements_section::read_entity: exactly 3 values expected on the line");    
  }

  /****************************************************************************/
  void msh_elements_section::read_element(const string &crkLine,
                                          size_t &rType,
                                          size_t &rN1, size_t &rN2, size_t &rN3, size_t &rN4,
                                          size_t &rEntity, size_t &rID) const
  {
    size_t nTags = 0;
    size_t nElementaryEntity = 0;
    size_t nUnused;
    stringstream kStream(crkLine);
    
    kStream >> rID;
    kStream >> rType;
    kStream >> nTags;
    for (size_t i = 0; i < nTags; ++i)
    {
      if (i == 0)
        kStream >> rEntity;
      else
        kStream >> nUnused;
    }
    
    kStream >> rN1;
    kStream >> rN2;
    if (rType > 1)
      kStream >> rN3;
    else
      rN3 = 0;
    if (rType > 2)
      kStream >> rN4;
    else
      rN4 = 0;
  }
    
  /****************************************************************************/
  bool msh_elements_section::read_structure(ifstream &rkMshFileStream)
  {
    eScanState eState = SS_read_header_dim;
    size_t nCurrPartId = 0;
    size_t nCurrBdId = 0;
    size_t nPhysicalEntities = 0;
    size_t nId = 0, nDim = 0;
    size_t nType, n1, n2, n3, n4, nEntity;
    size_t nIndex;
    string kName;
    string kLine;
    elements_line *p_line = NULL;
    elements_tri *p_tri = NULL;
    elements_tetra *p_tetra = NULL;
    
    while ((eState != SS_out_elements_section) && getline_checked(rkMshFileStream, kLine))
    {
      switch (eState)
      {
        case SS_read_header_dim:
           if (!kLine.compare(ms_DimHeader))
            eState = SS_read_nb_dim;
          break;
          
        case SS_read_nb_dim:
          nGeometryDim = getSizeNum(kLine);
          eState = SS_read_header_entities;
          break;
          
        case SS_read_header_entities:
          if (!kLine.compare(ms_PhysicalEntitiesHeader))
            eState = SS_read_nb_entities;
          break;
          
        case SS_read_nb_entities:
          eState = SS_read_entity;
          break;
                  
        case SS_read_entity:
          if (!kLine.compare(ms_PhysicalEntitiesEnd))
          {
            eState = SS_read_header_elements;
            break;
          }
          
          read_entity(kLine, nId, nDim, kName);
          if (nDim > 3)
            throw std::string("msh_elements_section::read_structure: physical entity has invalid dimension");
          if (nDim == nGeometryDim) {
            if (nDim == 2) {  // surface
              p_tri = new elements_tri;
              p_tri->m_nSectionIdx = nCurrPartId;
              p_tri->m_kSectionName = kName;
              p_tri->m_nLength = 0;
              m_akSetsTri.push_back(*p_tri);
              delete p_tri;
            } else { // volume
              p_tetra = new elements_tetra;
              p_tetra->m_nSectionIdx = nCurrPartId;
              p_tetra->m_kSectionName = kName;
              p_tetra->m_nLength = 0;
              m_akSetsTetra.push_back(*p_tetra);
              delete p_tetra;
            }
            m_anPartMapping.push_back(nId);  // Part <nId> is assigned to new Id <nCurrPartId>
            nCurrPartId++;
          } else {  // lower dimension (for BDs)
            if (nDim == 1) {  // line
              p_line = new elements_line;
              p_line->m_nSectionIdx = nCurrBdId;
              p_line->m_kSectionName = kName;
              p_line->m_nLength = 0;
              m_akSetsLine.push_back(*p_line);
              delete p_line;
            } else { // surface
              p_tri = new elements_tri;
              p_tri->m_nSectionIdx = nCurrBdId;
              p_tri->m_kSectionName = kName;
              p_tri->m_nLength = 0;
              m_akSetsTri.push_back(*p_tri);
              delete p_tri;
            }
            m_anBdPartMapping.push_back(nId);  // Part <nId> is assigned to new Id <nCurrBdId>
            nCurrBdId++;
          }
          break;

        case SS_read_header_elements:
          if (!kLine.compare(ms_ElementsHeader))
            eState = SS_read_nb_elements;
          break;
          
        case SS_read_nb_elements:
          eState = SS_read_element;
          break;

        case SS_read_element:
          if (!kLine.compare(ms_ElementsEnd))
          {
            eState = SS_out_elements_section;
            break;
          }

          read_element(kLine, nType, n1, n2, n3, n4, nEntity, nId);
          if ((nType == 2 && nGeometryDim == 2) || (nType == 4 && nGeometryDim == 3)) {
            nIndex = findIndex(m_anPartMapping, nEntity);
            if (nType == 2) // 3-node triangle
              m_akSetsTri[nIndex].m_nLength++;
            else  // 4-node tetrahedron
              m_akSetsTetra[nIndex].m_nLength++;
          } else {  // lower dimension (for BDs)
            nIndex = findIndex(m_anBdPartMapping, nEntity);
            if (nType == 1) // 2-node line
              m_akSetsLine[nIndex].m_nLength++;
            else  // 3-node triangle
              m_akSetsTri[nIndex].m_nLength++;
          }
          break;
          
        default:
          throw std::string("msh_elements_section::read_structure: eState contains an undefined state");
      } // switch (eState)
    } // while ((eState != SS_out_elements_section) && getline_checked(rkMshFileStream, kLine))
    
    if (eState != SS_out_elements_section) {
      if (eState == SS_read_header_dim)
        throw std::string("msh_elements_section::read_structure: header of DIMENSION section not found!");
      else
        throw std::string("msh_elements_section::read_structure: header/end of ELEMENTS section not found!");
    }
    
    for (size_t i = 0; i < m_akSetsLine.size(); ++i)
    {
      //mexPrintf("Elementset line '%s': idx=%d, quantity=%d.\n",
      //          m_akSetsLine[i].m_kSectionName.c_str(),
      //          m_akSetsLine[i].m_nSectionIdx,
      //          m_akSetsLine[i].m_nLength);
      create_vector_array_3d(m_akSetsLine[i].m_nLength,
                             m_akSetsLine[i].m_pfBuffer,
                             m_akSetsLine[i].m_pfId,
                             m_akSetsLine[i].m_pfN1,
                             m_akSetsLine[i].m_pfN2);
    }
    
    for (size_t i = 0; i < m_akSetsTri.size(); ++i)
    {
      //mexPrintf("Elementset tri '%s': idx=%d, quantity=%d.\n",
      //          m_akSetsTri[i].m_kSectionName.c_str(),
      //          m_akSetsTri[i].m_nSectionIdx,
      //          m_akSetsTri[i].m_nLength);
      create_vector_array_3d(m_akSetsTri[i].m_nLength,
                             m_akSetsTri[i].m_pfBuffer,
                             m_akSetsTri[i].m_pfN1,
                             m_akSetsTri[i].m_pfN2,
                             m_akSetsTri[i].m_pfN3);
      create_vector_array_4d(m_akSetsTri[i].m_nLength,
                             m_akSetsTri[i].m_pfBuffer,
                             m_akSetsTri[i].m_pfId,
                             m_akSetsTri[i].m_pfN1,
                             m_akSetsTri[i].m_pfN2,
                             m_akSetsTri[i].m_pfN3);
    }

    for (size_t i=0; i<m_akSetsTetra.size(); ++i)
    {
      //mexPrintf("Elementset tetra '%s': idx=%d, quantity=%d.\n",
      //          m_akSetsTetra[i].m_kSectionName.c_str(),
      //          m_akSetsTetra[i].m_nSectionIdx,
      //          m_akSetsTetra[i].m_nLength);
      create_vector_array_5d(m_akSetsTetra[i].m_nLength,
                             m_akSetsTetra[i].m_pfBuffer,
                             m_akSetsTetra[i].m_pfId,
                             m_akSetsTetra[i].m_pfN1,
                             m_akSetsTetra[i].m_pfN2,
                             m_akSetsTetra[i].m_pfN3,
                             m_akSetsTetra[i].m_pfN4);
    }
    
    return true;
  }
  
  /****************************************************************************/
  void msh_elements_section::read_contents(ifstream &rkMshFileStream)
  {
    eScanState eState = SS_read_header_elements;
    size_t nType, n1, n2, n3, n4, nEntity, nId;
    size_t nIndex;
    int nIdOffset = -1;
    int nBdOffset = -1;
    string kLine;
    std::vector<size_t> anCurrEl;
    std::vector<size_t> anCurrBd;
    if (nGeometryDim == 2) {
      anCurrEl = std::vector<size_t>(m_akSetsTri.size(), 0);
      anCurrBd = std::vector<size_t>(m_akSetsLine.size(), 0);
    } else {
      anCurrEl = std::vector<size_t>(m_akSetsTetra.size(), 0);
      anCurrBd = std::vector<size_t>(m_akSetsTri.size(), 0);
    }
      
    while ((eState != SS_out_elements_section) && getline_checked(rkMshFileStream, kLine))
    {
      switch (eState)
      {
        case SS_read_header_elements:
          if (!kLine.compare(ms_ElementsHeader))
            eState = SS_read_nb_elements;
          break;

        case SS_read_nb_elements:
          eState = SS_read_element;
          break;
          
        case SS_read_element:
          if (!kLine.compare(ms_ElementsEnd))
          {
            eState = SS_out_elements_section;
            break;
          }

          read_element(kLine, nType, n1, n2, n3, n4, nEntity, nId);
          if ((nType == 2 && nGeometryDim == 2) || (nType == 4 && nGeometryDim == 3)) {
            if (nIdOffset == -1)
              nIdOffset = nId - 1;
            nIndex = findIndex(m_anPartMapping, nEntity);
            if (nType == 2) { // 3-node triangle
              m_akSetsTri[nIndex].m_pfN1[anCurrEl[nIndex]] = n1;
              m_akSetsTri[nIndex].m_pfN2[anCurrEl[nIndex]] = n2;
              m_akSetsTri[nIndex].m_pfN3[anCurrEl[nIndex]] = n3;
              m_akSetsTri[nIndex].m_pfId[anCurrEl[nIndex]] = nId - nIdOffset;
            } else {  // 4-node tetrahedron
              m_akSetsTetra[nIndex].m_pfN1[anCurrEl[nIndex]] = n1;
              m_akSetsTetra[nIndex].m_pfN2[anCurrEl[nIndex]] = n2;
              m_akSetsTetra[nIndex].m_pfN3[anCurrEl[nIndex]] = n3;
              m_akSetsTetra[nIndex].m_pfN4[anCurrEl[nIndex]] = n4;
              m_akSetsTetra[nIndex].m_pfId[anCurrEl[nIndex]] = nId - nIdOffset;
            }
            anCurrEl[nIndex]++;
          } else {  // lower dimension (for BDs)
            if (nBdOffset == -1)
              nBdOffset = nId - 1;
            nIndex = findIndex(m_anBdPartMapping, nEntity);
            if (nType == 1) { // 2-node line
              m_akSetsLine[nIndex].m_pfN1[anCurrBd[nIndex]] = n1;
              m_akSetsLine[nIndex].m_pfN2[anCurrBd[nIndex]] = n2;
              m_akSetsLine[nIndex].m_pfId[anCurrBd[nIndex]] = nId - nBdOffset;
            } else {  // 3-node triangle
              m_akSetsTri[nIndex].m_pfN1[anCurrBd[nIndex]] = n1;
              m_akSetsTri[nIndex].m_pfN2[anCurrBd[nIndex]] = n2;
              m_akSetsTri[nIndex].m_pfN3[anCurrBd[nIndex]] = n3;
              m_akSetsTri[nIndex].m_pfId[anCurrBd[nIndex]] = nId - nBdOffset;
            }
            anCurrBd[nIndex]++;
          }
          break;
                   
        default:
          throw std::string("msh_elements_section::read_contents: eState contains an undefined state");  
      } // switch (eState)
    }
  }
  
  /****************************************************************************/
  void msh_elements_section::find_sidesets()
  {
    size_t N1, N2;
    size_t nIndex;
    std::vector<double> vnE1;
    std::vector<double> vnE2;
    std::vector<double> vnE3;
    surface* pSurface = NULL;
    sideset* pSideset = NULL;
    if (nGeometryDim == 2) {
      for (size_t i = 0; i < m_akSetsLine.size(); ++i) {
        pSurface = new surface;
        pSurface->m_kSectionName = m_akSetsLine[i].m_kSectionName;
        vnE1.clear();
        vnE2.clear();
        vnE3.clear();
        for (size_t j = 0; j < m_akSetsLine[i].m_nLength; ++j) {
          N1 = (size_t)m_akSetsLine[i].m_pfN1[j];
          N2 = (size_t)m_akSetsLine[i].m_pfN2[j];
          findElementByLine(N1, N2, vnE1, vnE2, vnE3);
        }
        
        pSideset = new sideset;
        pSideset->m_eType = sideset::edge;
        pSideset->m_kSectionName = pSurface->m_kSectionName + "_E1";
        pSideset->m_nLength = vnE1.size();
        if (vnE1.size() > 0) {
          create_vector_array_1d(vnE1.size(), pSideset->m_pfBuffer);
          memcpy(pSideset->m_pfBuffer, vnE1.data(), vnE1.size() * sizeof(double));
        }
        pSurface->m_akSidesets.push_back(*pSideset);
        delete pSideset;

        pSideset = new sideset;
        pSideset->m_eType = sideset::edge;
        pSideset->m_kSectionName = pSurface->m_kSectionName + "_E2";
        pSideset->m_nLength = vnE2.size();
        if (vnE2.size() > 0) {
          create_vector_array_1d(vnE2.size(), pSideset->m_pfBuffer);
          memcpy(pSideset->m_pfBuffer, vnE2.data(), vnE2.size() * sizeof(double));
        }
        pSurface->m_akSidesets.push_back(*pSideset);
        delete pSideset;

        pSideset = new sideset;
        pSideset->m_eType = sideset::edge;
        pSideset->m_kSectionName = pSurface->m_kSectionName + "_E3";
        pSideset->m_nLength = vnE3.size();
        if (vnE3.size() > 0) {
          create_vector_array_1d(vnE3.size(), pSideset->m_pfBuffer);
          memcpy(pSideset->m_pfBuffer, vnE3.data(), vnE3.size() * sizeof(double));
        }
        pSurface->m_akSidesets.push_back(*pSideset);
        delete pSideset;

        m_akSurfaces.push_back(*pSurface);
        delete pSurface;
      }
    } else {    // 3D
      for (size_t i = 0; i < m_akSetsTri.size(); ++i) {
        pSurface = new surface;
        pSurface->m_kSectionName = m_akSetsTri[i].m_kSectionName;       
        pSideset = new sideset;
        pSideset->m_eType = sideset::surface;
        pSideset->m_kSectionName = pSurface->m_kSectionName + "_S1";
        pSideset->m_nLength = m_akSetsTri[i].m_nLength;
        create_vector_array_1d(m_akSetsTri[i].m_nLength, pSideset->m_pfBuffer);
        for (size_t j = 0; j < m_akSetsTri[i].m_nLength; ++j)
          pSideset->m_pfBuffer[j] = m_akSetsTri[i].m_pfId[j];
        pSurface->m_akSidesets.push_back(*pSideset);
        delete pSideset;  
        m_akSurfaces.push_back(*pSurface);
        delete pSurface;
      }
    }
  }
  
  /****************************************************************************/
  mxArray* msh_elements_section::scan(ifstream &rkMshFileStream)
  {
    if (!read_structure(rkMshFileStream)) return mxCreateCellArray(0, NULL);
    rkMshFileStream.clear();
    rkMshFileStream.seekg(0);
    read_contents(rkMshFileStream);
    find_sidesets();
    return create_and_fill_matlab_structure();
  }
  
  /****************************************************************************/
  mxArray* msh_elements_section::scan_bd()
  {
    if (m_akSurfaces.size() == 0) return mxCreateCellArray(0, NULL);
/*    for (size_t nCurrSurface = 0; nCurrSurface < m_akSurfaces.size(); ++nCurrSurface)
    {
      mexPrintf("surface name = '%s', #sidesets = %d\n", m_akSurfaces[nCurrSurface].m_kSectionName.c_str(), m_akSurfaces[nCurrSurface].m_akSidesets.size());      
      for (size_t nSideset = 0; nSideset < m_akSurfaces[nCurrSurface].m_akSidesets.size(); ++nSideset)
      {
        sideset kSideset = m_akSurfaces[nCurrSurface].m_akSidesets[nSideset];
        mexPrintf("  sideset name = '%s', #elements = %d\n", kSideset.m_kSectionName.c_str(), kSideset.m_nLength);
      }
    } */
    return create_and_fill_matlab_structure_bd();
  }

  /****************************************************************************/
  
}; // namespace msh_n