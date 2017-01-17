/**
 * \file
 * \author Michael Dudzinski
 * \date 10.09.2014
 */

#include "inp_assembly_section.hpp"
#include "inp_utils.hpp"

using namespace std;

namespace inp_n {

  const string inp_assembly_section::ms_SectionHeader = "********************************** A S S E M B L Y ************************************";

  /****************************************************************************/
  inp_assembly_section::inp_assembly_section()
  :
  m_akMaterials(0)
  {
    /**/
  }
  /****************************************************************************/
  mxArray*
  inp_assembly_section::create_and_fill_matlab_structure()
  {
    if (m_akMaterials.size()==0)
      return mxCreateCellArray(0, NULL);

    const mwSize nSetsQuantity = m_akMaterials.size();

    mwSize anDims[2] = { 2, nSetsQuantity };
    mxArray *pkCell  = mxCreateCellArray(2,anDims);

    size_t i=0, j=0;
    double *pfBuffer = NULL;

    vector<material>::iterator pMatIter;
    vector<material::property>::iterator pPropIter;

    for (pMatIter=m_akMaterials.begin(), i=0; pMatIter!=m_akMaterials.end(); ++pMatIter, ++i)
    {
      /* name of material */
      mxArray *pkMatName = mxCreateString(pMatIter->m_kName.c_str());
      mxSetCell(pkCell,i*2,pkMatName);

      /*  */
      const mwSize nPropsQuantity = pMatIter->m_akProperties.size();
      mwSize anPropDims[2]        = { 2, nPropsQuantity };
      mxArray *pkMaterialCell     = mxCreateCellArray(2,anPropDims);

      for (pPropIter=pMatIter->m_akProperties.begin(), j=0; pPropIter!=pMatIter->m_akProperties.end(); ++pPropIter, ++j)
      {
        /* name of material property */
        mxArray *pkPropName = mxCreateString(pPropIter->m_kName.c_str());
        mxSetCell(pkMaterialCell,j*2,pkPropName);

        /* contents of property */
        mxArray *pkElemsMatrix = mxCreateNumericMatrix(pPropIter->m_akValues.size(),
                                                       1                           ,
                                                       mxDOUBLE_CLASS              ,
                                                       mxREAL                      );
        copy(pPropIter->m_akValues.begin(),
             pPropIter->m_akValues.end()  ,
             mxGetPr(pkElemsMatrix)       );

        mxSetCell(pkMaterialCell,j*2+1,pkElemsMatrix);
      }

      mxSetCell(pkCell,i*2+1,pkMaterialCell);
    }
    
    return pkCell;
  }
  /****************************************************************************/
  void
  inp_assembly_section::printState(const eScanState &eState,const std::string &crkLine)
  const
  {
    switch (eState)
    {
      case SS_in_assembly_section:
        mexPrintf("current state: SS_in_assembly_section\n");
        break;

      case SS_read_material_header:
        mexPrintf("current state: SS_read_material_header\n");
        break;

      case SS_read_property_header:
        mexPrintf("current state: SS_read_property_header\n");
        break;

      case SS_read_property:
        mexPrintf("current state: SS_read_property\n");
        break;

      case SS_out_assembly_section:
        mexPrintf("current state: SS_out_assembly_section\n");
        break;
    }
    mexPrintf("current line: \"%s\"\n",crkLine.c_str());
  }
  /****************************************************************************/
  void
  inp_assembly_section::read_material_header(const string &crkLine      ,
                                             string       &rkSectionName)
  const
  {
    size_t nPos    = crkLine.find("NAME =")+6;
//    size_t nLength = crkLine.find_last_not_of(WhitespaceString+",",nPos);
//    if (nLength!=string::npos) ++nLength;
    size_t nLength = crkLine.substr(nPos).find_first_of(WhitespaceString+",",nPos);

    rkSectionName = crkLine.substr(nPos,nLength);
    trim(rkSectionName);

//    mexPrintf("Material found: %s\n",rkSectionName.c_str());
  }
  /****************************************************************************/
  void
  inp_assembly_section::read_property_header(const std::string &crkLine,
                                             std::string &rkSectionName)
  const
  {
//    mexPrintf("Line: %s\n",crkLine.c_str());
    size_t nPos    = crkLine.find("*")+1;
    size_t nLength = crkLine.substr(nPos).find_first_of(",");
//    if (nLength!=string::npos) nLength-=nPos;

//    mexPrintf("nPos=%d, nLength=%d, string::npos=%d\n",nPos,nLength,string::npos);
    
    rkSectionName = crkLine.substr(nPos,nLength);
    trim(rkSectionName);

//    mexPrintf("Property found: %s\n",rkSectionName.c_str());
  }
  /****************************************************************************/
  void
  inp_assembly_section::read_property(const string        &crkLine  ,
                                      std::vector<double> &rafValues)
  const
  {
//    mexPrintf("Entering inp_assembly_section::read_property\n");
//
//    mexPrintf("rafValues: ");
//    for (size_t i=0; i<rafValues.size(); ++i)
//      mexPrintf("%f ",rafValues[i]);
//    mexPrintf("\n");

    read_numeric_values(crkLine,rafValues);

//    mexPrintf("rafValues (after reading): ");
//    for (size_t i=0; i<rafValues.size(); ++i)
//      mexPrintf("%f ",rafValues[i]);
//    mexPrintf("\n");
  }
  /****************************************************************************/
  bool
  inp_assembly_section::read_structure(ifstream &rkInpFileStream)
  {
    eScanState eState=SS_in_assembly_section;
    std::streampos nPrevStreamPos=rkInpFileStream.tellg();

    size_t nCurrSetQuantity=0;
    size_t nCurrSet=0;

    string kLine;

    /* search for material data, everything else is skipped for now */
    while (getline_checked(rkInpFileStream,kLine))
    {
//      mexPrintf("read structure line  : %s.\n",kLine.c_str());

      if (!kLine.find("*MATERIAL"))
      {
        rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
        eState = SS_read_material_header;
        break;
      }

      /* next section found */
      else if (!kLine.find("***"))
      {
        rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
        eState = SS_out_assembly_section;
        break;
      }

      nPrevStreamPos=rkInpFileStream.tellg();
    }

    if (rkInpFileStream.eof())
      eState = SS_out_assembly_section;

    
    while ((eState!=SS_out_assembly_section) && getline_checked(rkInpFileStream,kLine))
    {
//      printState(eState,kLine);
      
      /* check for scanner state */
      switch (eState)
      {
        case SS_in_assembly_section:
        {
          if (!kLine.find("*MATERIAL"))
          {
            rkInpFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_material_header;
          }
          else if (!kLine.find("**"))
          {
            eState = SS_out_assembly_section;
          }

          break;
        }

        case SS_read_material_header:
        {
          m_akMaterials.resize(m_akMaterials.size()+1);
          read_material_header(kLine,m_akMaterials.back().m_kName);
          
          eState=SS_read_property_header;

          break;
        }

        case SS_read_property_header:
        {
          if (!kLine.compare("**"))
          {
            if (m_akMaterials.back().m_akProperties.size()==0)
              /* empty material found */
              throw std::string("inp_assembly_section::read_structure: empty material found");

            eState = SS_in_assembly_section;
          }
          else if (!kLine.find("*"))
          {
            m_akMaterials.back().m_akProperties.resize(m_akMaterials.back().m_akProperties.size()+1);
            read_property_header(kLine,m_akMaterials.back().m_akProperties.back().m_kName);

            eState = SS_read_property;
          }
          else
          {
            throw std::string("inp_assembly_section::read_structure: file syntax corrupted. I expected either a material property or end of section, i.e. \"**\"");
          }

          break;
        }

        case SS_read_property:
        {
          if (!kLine.find("*"))
            throw std::string("inp_assembly_section::read_structure: property of material must contain only numeric values");

          read_property(kLine,m_akMaterials.back().m_akProperties.back().m_akValues);

//          mexPrintf("Material %s: Property=%s, Values=",m_akMaterials.back().m_kName.c_str(),m_akMaterials.back().m_akProperties.back().m_kName.c_str());
//
//          for (size_t i=0; i<m_akMaterials.back().m_akProperties.back().m_akValues.size(); ++i)
//            mexPrintf("%f ",m_akMaterials.back().m_akProperties.back().m_akValues[i]);
//          mexPrintf("\n");

          eState = SS_read_property_header;

          break;
        }
          
        default:
          throw std::string("inp_assembly_section::read_structure: eState contains an undefined state");
          
      } /* switch (eState) */

      nPrevStreamPos=rkInpFileStream.tellg();
    } /* while ((eState!=SS_out_nodes_section) && getline_checked(rkInpFileStream,kLine)) */
    
    if (eState != SS_out_assembly_section)
      throw std::string("inp_assembly_section::read_structure: the end of the ***...* A S S E M B L Y *...*** section is missing");

//    vector<material>::iterator pMatIter;
//    vector<material::property>::iterator pPropIter;
//
//    for (pMatIter=m_akMaterials.begin(); pMatIter!=m_akMaterials.end(); ++pMatIter)
//    {
//      /* name of material */
//      mexPrintf("Material: Name=%s\n",pMatIter->m_kName.c_str());
//
//      for (pPropIter=pMatIter->m_akProperties.begin(); pPropIter!=pMatIter->m_akProperties.end(); ++pPropIter)
//      {
//        mexPrintf("\t%s: Values=",pPropIter->m_kName.c_str());
//        for (size_t i=0; i<pPropIter->m_akValues.size(); ++i)
//          mexPrintf("%f ",pPropIter->m_akValues[i]);
//        mexPrintf("\n");
//      }
//    }
    
    return true;
  }
  /****************************************************************************/
  mxArray*
  inp_assembly_section::scan(ifstream &rkInpFileStream)
  {
//    mexPrintf("Assembly section\n");

//    mexPrintf("stream position: %d.\n",(int) rkInpFileStream.tellg());
//    streampos nSectionStart = rkInpFileStream.tellg();

//    if (!read_structure(rkInpFileStream)) return mxCreateCellArray(0, NULL);

    read_structure(rkInpFileStream);
    
//    rkInpFileStream.clear();
//    rkInpFileStream.seekg(nSectionStart);
//
//    read_contents(rkInpFileStream);

//    return mxCreateCellArray(0, NULL);

    return create_and_fill_matlab_structure();
  }
  /****************************************************************************/
  
}; /* namespace inp_n */
