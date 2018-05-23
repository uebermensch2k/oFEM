/**
 * \file
 * \author Michael Hagel
 * \date 15.05.2018
 */

#include "msh_physical_section.hpp"
#include "msh_utils.hpp"
#include <algorithm>					// That's why we dont use namespace std!
#include <memory.h>

using namespace std;

namespace msh_n {

  const string msh_physical_section::ms_SectionHeader = "$PhysicalNames";

  /****************************************************************************/
  msh_physical_section::msh_physical_section()
  :
  m_akSetsPhysical (0)
  {
    /**/
  }
  /****************************************************************************/

  void msh_physical_section::read_physical (std::string &crkLine,
	   										std::string &physicalName,
											int &Id){
		stringstream kStream(crkLine);
		int tmp=0;

		kStream>>tmp;
		kStream>>tmp;
		Id = tmp;
		kStream>>physicalName;
		physicalName.erase(std::remove(physicalName.begin(), physicalName.end(), '"'), physicalName.end());

  }

  /****************************************************************************/

  mxArray*
  msh_physical_section::create_and_fill_matlab_structure()
  {
    const mwSize nSetsQuantity = m_akSetsPhysical.size();

    mwSize anDims[2] = { 2, nSetsQuantity };
    mxArray *pkCell  = mxCreateCellArray(2,anDims);
	int tmp;

    size_t nCurrSetPhysical=0;

    for (size_t i=0; i<nSetsQuantity; ++i)
    {

		if (m_akSetsPhysical.size())
		{
		  mxArray *pkNameMatrix = mxCreateString(m_akSetsPhysical[i].m_kPhysicalName.c_str());
		  mxSetCell(pkCell,i*2,pkNameMatrix);

		  /* content of nodeset */
		  mxArray *pkIdMatrix = mxCreateNumericMatrix(1, 1, mxINT8_CLASS, mxREAL);

		  memcpy(mxGetPr(pkIdMatrix),&m_akSetsPhysical[i].m_nPhysicalIdx, sizeof(int));

		  mxSetCell(pkCell,i*2+1,pkIdMatrix);
		  ++nCurrSetPhysical;
		}
    }

    return pkCell;
  }
  /****************************************************************************/
  bool
  msh_physical_section::read_structure(ifstream &rkMshFileStream)
  {
//    mexPrintf("Reading structure\n");

    eScanState eState=SS_read_set_header;
    std::streampos nPrevStreamPos=rkMshFileStream.tellg();

    size_t nCurrSetQuantity=0;
    size_t nCurrSet=0;

    string kLine;

    while ((eState!=SS_out_physical_section) && getline_checked(rkMshFileStream,kLine))
    {
//      mexPrintf("line  : %s.\n",kLine.c_str());

      /* check for scanner state */
      switch (eState)
      {
        case SS_in_physical_section:
        {
          if (!kLine.find("$PhysicalNames"))
          {
            rkMshFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
//            mexPrintf("Switching from SS_in_elements_section to SS_read_set_header\n");
            break;
          }

          if (!kLine.compare("$EndPhysicalNames"))
          {
            return false; /* signal empty section */
          }

          break;
        }

        case SS_read_set_header:
        {
//          mexPrintf("reading header: %s\n",kLine.c_str());
          nCurrSet         = m_akSetsPhysical.size();
		  nCurrSetQuantity = 0;
		  std::vector<double> afValues;

          std::string kSectionName, kSectionType;
          read_numeric_values(kLine,afValues);

          /* dependent on element switch state to read tri or tetra */
          if (afValues.size()==1){
			  eState = SS_read_physical;
          }
          else
            throw std::string("msh_elements_section::read_structure: element line must contain exactly 4 or exactly 5 values");

          break;
        }

        case SS_read_physical:
        {
//          mexPrintf("reading element: %s\n",kLine.c_str());

          if (!kLine.find("$PhysicalNames"))
          {
            rkMshFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
//            mexPrintf("Switching from SS_read_element to SS_read_set_header\n");
            break;
          }

          if (!kLine.compare("$EndPhysicalNames"))
          {
            //m_akSetsTri.back().m_nLength=nCurrSetQuantity;

            eState = SS_out_physical_section;
//            mexPrintf("Switching from SS_read_element to SS_out_elements_section\n");
            break;
          }
		  	m_akSetsPhysical.resize(m_akSetsPhysical.size()+1);

          break;
        }

        default:
          throw std::string("msh_elements_section::read_structure: eState contains an undefined state");

      } /* switch (eState) */

      nPrevStreamPos=rkMshFileStream.tellg();
    } /* while ((eState!=SS_out_elements_section) && getline_checked(rkMshFileStream,kLine)) */

    if (eState != SS_out_physical_section)
      throw std::string("msh_elements_section::read_structure: the end of the $Elements section is missing");

    return true;
  }
  /****************************************************************************/
  void
  msh_physical_section::read_contents(ifstream &rkMshFileStream)
  {
//    mexPrintf("Reading content\n");

    eScanState eState=SS_read_set_header;
    std::streampos nPrevStreamPos=rkMshFileStream.tellg();

    size_t nCurrSetQuantity=0;

    int Id;
	string physName;

    string kLine;

//    getline_checked(rkMshFileStream,kLine);
//    rkMshFileStream.seekg(nPrevStreamPos); /* line to be reread */
//    mexPrintf("line: %s\n",kLine.c_str());

    while ((eState!=SS_out_physical_section) && getline_checked(rkMshFileStream,kLine))
    {
      /* check for scanner state */
      switch (eState)
      {
        case SS_in_physical_section:
        {
          if (!kLine.find("$PhysicalNames"))
          {
            rkMshFileStream.seekg(nPrevStreamPos); /* line to be reread */
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
		  std::vector<double> afValues;

          std::string kSectionName, kSectionType;
          read_numeric_values(kLine,afValues);
		  //mexPrintf("%s/n",kLine.c_str());

          if (afValues.size()==1){
			  eState = SS_read_physical;
          }
          else
            throw std::string("msh_elements_section::read_structure: element line must contain exactly 4 or exactly 5 values");

          break;
        }

        case SS_read_physical:
        {
          if (!kLine.find("$PhysicalNames"))
          {
            rkMshFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
//            mexPrintf("Switching from SS_read_element to SS_read_set_header\n");
            break;
          }

          if (!kLine.compare("$EndPhysicalNames"))
          {
            eState = SS_out_physical_section;
//            mexPrintf("Switching from SS_read_element to SS_out_elements_section\n");
            break;
          }
		  else
		  {
		  	read_physical(kLine,physName,Id);
			m_akSetsPhysical[nCurrSetQuantity].m_nPhysicalIdx = Id;
			m_akSetsPhysical[nCurrSetQuantity].m_kPhysicalName = physName;
			++nCurrSetQuantity;


          break;
          }
	    }

        default:
          throw std::string("msh_elements_section::read_structure: eState contains an undefined state");

      } /* switch (eState) */

      nPrevStreamPos=rkMshFileStream.tellg();
    } /* while ((eState!=SS_out_elements_section) && getline_checked(rkMshFileStream,kLine)) */
  }
  /****************************************************************************/
  mxArray*
  msh_physical_section::scan(ifstream &rkMshFileStream)
  {
    streampos nSectionStart = rkMshFileStream.tellg();

    if (!read_structure(rkMshFileStream)) return mxCreateCellArray(0, NULL);

    rkMshFileStream.clear();
    rkMshFileStream.seekg(nSectionStart);

//    mexPrintf("Reading content\n");

    read_contents(rkMshFileStream);

    return create_and_fill_matlab_structure();
//    return mxCreateCellArray(0, NULL);
  }
  /****************************************************************************/

}; /* namespace msh_n */
