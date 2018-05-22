/**
 * \file
 * \author Michael Dudzinski
 * \date 10.09.2014
 */

#include "msh_elements_section.hpp"
#include "msh_utils.hpp"

using namespace std;

namespace msh_n {

  const string msh_elements_section::ms_SectionHeader = "$Elements";

  /****************************************************************************/
  msh_elements_section::msh_elements_section()
  :
  m_akSetsLine (0),
  m_akSetsTri  (0),
  m_akSetsTetra(0)
  {
    /**/
  }
  /****************************************************************************/
  mxArray*
  msh_elements_section::create_and_fill_matlab_structure()
  {
    const mwSize nSetsQuantity = 3;

    mwSize anDims[2] = { 2, nSetsQuantity };
    mxArray *pkCell  = mxCreateCellArray(2,anDims);

    size_t nCurrSetLine=0,nCurrSetTri=0,nCurrSetTetra=0;
	string line = "Line";
	string tri = "Tri";
	string tetra = "Tetra";

    for (size_t i=0; i<nSetsQuantity; ++i)
    {

		if (m_akSetsLine.size() && m_akSetsLine[0].m_nLength)
		{

		  mxArray *pkNameMatrix = mxCreateString(line.c_str());
	      mxSetCell(pkCell,0*2,pkNameMatrix);

		  /* content of nodeset */
		  mxArray *pkElemsMatrix = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

		  mxSetData(pkElemsMatrix, m_akSetsLine[0].m_pfBuffer);
		  mxSetM   (pkElemsMatrix, m_akSetsLine[0].m_nLength );
		  mxSetN   (pkElemsMatrix, 3                                    );

		  mxSetCell(pkCell,0*2+1,pkElemsMatrix);
		  ++i;

		  ++nCurrSetLine;
		}

      if (m_akSetsTri.size() && m_akSetsTri[0].m_nLength)
      {
		  mxArray *pkNameMatrix = mxCreateString(tri.c_str());
	      mxSetCell(pkCell,1*2,pkNameMatrix);

        /* content of nodeset */
        mxArray *pkElemsMatrix = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

        mxSetData(pkElemsMatrix, m_akSetsTri[0].m_pfBuffer);
        mxSetM   (pkElemsMatrix, m_akSetsTri[0].m_nLength );
        mxSetN   (pkElemsMatrix, 4                                  );

        mxSetCell(pkCell,1*2+1,pkElemsMatrix);
		++i;

        ++nCurrSetTri;
      }

      if (m_akSetsTetra.size() && m_akSetsTetra[0].m_nLength)
      {
		  mxArray *pkNameMatrix = mxCreateString(tetra.c_str());
	      mxSetCell(pkCell,2*2,pkNameMatrix);


        /* content of nodeset */
        mxArray *pkElemsMatrix = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

        mxSetData(pkElemsMatrix, m_akSetsTetra[0].m_pfBuffer);
        mxSetM   (pkElemsMatrix, m_akSetsTetra[0].m_nLength );
        mxSetN   (pkElemsMatrix, 5                                      );

        mxSetCell(pkCell,2*2+1,pkElemsMatrix);
		++i;

        ++nCurrSetTetra;
      }
    }

    return pkCell;
  }
  /****************************************************************************/
  void
  msh_elements_section::read_set_header(const string &crkLine      ,
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
  msh_elements_section::read_element_line(const string &crkLine,
										 double       &rfN1   ,
										 double       &rfN2   ,
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

  }

  void
  msh_elements_section::read_element_tri(const string &crkLine,
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
  msh_elements_section::read_element_tetra(const string &crkLine,
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
  msh_elements_section::read_structure(ifstream &rkMshFileStream)
  {
//    mexPrintf("Reading structure\n");

    eScanState eState=SS_read_set_header;
    std::streampos nPrevStreamPos=rkMshFileStream.tellg();

    size_t nCurrSetQuantity=0;
    size_t nCurrSet=0;
	size_t nLine=0;
	size_t nTri=0;
	size_t nTet=0;

    string kLine;

    while ((eState!=SS_out_elements_section) && getline_checked(rkMshFileStream,kLine))
    {
//      mexPrintf("line  : %s.\n",kLine.c_str());

      /* check for scanner state */
      switch (eState)
      {
        case SS_in_elements_section:
        {
          if (!kLine.find("$Elements"))
          {
            rkMshFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
//            mexPrintf("Switching from SS_in_elements_section to SS_read_set_header\n");
            break;
          }

          if (!kLine.compare("$EndElements"))
          {
            return false; /* signal empty section */
          }

          break;
        }

        case SS_read_set_header:
        {
//          mexPrintf("reading header: %s\n",kLine.c_str());
          nCurrSet         = m_akSetsLine.size()+m_akSetsTri.size()+m_akSetsTetra.size();
		  nCurrSetQuantity = 0;
		  std::vector<double> afValues;

          std::string kSectionName, kSectionType;
          read_numeric_values(kLine,afValues);

          /* dependent on element switch state to read tri or tetra */
          if (afValues.size()==1){
			  eState = SS_read_element;
          }
          else
            throw std::string("msh_elements_section::read_structure: element line must contain exactly 4 or exactly 5 values");

          break;
        }

        case SS_read_element:
        {
//          mexPrintf("reading element: %s\n",kLine.c_str());

          if (!kLine.find("$Elements"))
          {
            m_akSetsTri.back().m_nLength=nCurrSetQuantity;

            rkMshFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
//            mexPrintf("Switching from SS_read_element to SS_read_set_header\n");
            break;
          }

          if (!kLine.compare("$EndElements"))
          {
            //m_akSetsTri.back().m_nLength=nCurrSetQuantity;
			m_akSetsLine.resize(m_akSetsLine.size()+1);
			m_akSetsTri.resize(m_akSetsTri.size()+1);
			m_akSetsTetra.resize(m_akSetsTetra.size()+1);
			m_akSetsLine.back().m_nLength=nLine;
			m_akSetsTri.back().m_nLength=nTri;
			m_akSetsTetra.back().m_nLength=nTet;

            eState = SS_out_elements_section;
//            mexPrintf("Switching from SS_read_element to SS_out_elements_section\n");
            break;
          }

		  std::vector<int> afValues;
			read_numeric_values(kLine,afValues);
		  switch(afValues[1]){
			  case Type_line:
				++nLine;
			  break;
			  case Type_triangle:
				++nTri;
			  break;
			  case Type_tetra:
				++nTet;
			  break;
		  };

          break;
        }

        default:
          throw std::string("msh_elements_section::read_structure: eState contains an undefined state");

      } /* switch (eState) */

      nPrevStreamPos=rkMshFileStream.tellg();
    } /* while ((eState!=SS_out_elements_section) && getline_checked(rkMshFileStream,kLine)) */

    if (eState != SS_out_elements_section)
      throw std::string("msh_elements_section::read_structure: the end of the $Elements section is missing");

    return true;
  }
  /****************************************************************************/
  void
  msh_elements_section::read_contents(ifstream &rkMshFileStream)
  {
//    mexPrintf("Reading content\n");

    eScanState eState=SS_read_set_header;
    std::streampos nPrevStreamPos=rkMshFileStream.tellg();

    size_t nCurrSetQuantity=0;
	int    nCurrSetLine=-1;
    int    nCurrSetTri=-1;
    int    nCurrSetTetra=-1;
	int    nCurrLine=0;
	int    nCurrTri=0;
	int    nCurrTet=0;

    double n1,n2,n3,n4, ID;

    string kLine;

//    getline_checked(rkMshFileStream,kLine);
//    rkMshFileStream.seekg(nPrevStreamPos); /* line to be reread */
//    mexPrintf("line: %s\n",kLine.c_str());

    while ((eState!=SS_out_elements_section) && getline_checked(rkMshFileStream,kLine))
    {

      /* check for scanner state */
      switch (eState)
      {
        case SS_in_elements_section:
        {
          if (!kLine.find("$Elements"))
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
//		  mexPrintf("%s/n",kLine.c_str());

          if (afValues.size()==1){
			  eState = SS_read_element;
			  if(m_akSetsLine[0].m_nLength>0)
			  	create_vector_array_3d(m_akSetsLine[0].m_nLength ,
                                     m_akSetsLine[0].m_pfBuffer,
                                     m_akSetsLine[0].m_pfId    ,
                                     m_akSetsLine[0].m_pfN1    ,
                                     m_akSetsLine[0].m_pfN2    );
			if(m_akSetsTri[0].m_nLength>0)
				create_vector_array_4d(m_akSetsTri[0].m_nLength  ,
								   m_akSetsTri[0].m_pfBuffer ,
								   m_akSetsTri[0].m_pfId     ,
								   m_akSetsTri[0].m_pfN1     ,
								   m_akSetsTri[0].m_pfN2     ,
							       m_akSetsTri[0].m_pfN3     );
			if(m_akSetsTetra[0].m_nLength>0)
		   		create_vector_array_5d(m_akSetsTetra[0].m_nLength  ,
								   m_akSetsTetra[0].m_pfBuffer,
								   m_akSetsTetra[0].m_pfId    ,
								   m_akSetsTetra[0].m_pfN1    ,
								   m_akSetsTetra[0].m_pfN2    ,
								   m_akSetsTetra[0].m_pfN3    ,
								   m_akSetsTetra[0].m_pfN4    );
          }
          else
            throw std::string("msh_elements_section::read_structure: element line must contain exactly 4 or exactly 5 values");

          break;
        }

        case SS_read_element:
        {
          if (!kLine.find("$Elements"))
          {
            rkMshFileStream.seekg(nPrevStreamPos); /* line to be reread */
            eState = SS_read_set_header;
//            mexPrintf("Switching from SS_read_element to SS_read_set_header\n");
            break;
          }

          if (!kLine.compare("$EndElements"))
          {
            eState = SS_out_elements_section;
//            mexPrintf("Switching from SS_read_element to SS_out_elements_section\n");
            break;
          }
		  std::vector<int> afValues;
          read_numeric_values(kLine,afValues);
		  int nInts = int(afValues[2]);
		  switch(afValues[1]){
			  case Type_line:
				m_akSetsLine[0].m_pfN1[nCurrLine]=afValues[nInts+3];
				m_akSetsLine[0].m_pfN2[nCurrLine]=afValues[nInts+4];
				m_akSetsLine[0].m_pfId[nCurrLine]=afValues[3];
				++nCurrLine;
			  break;
			  case Type_triangle:
				  m_akSetsTri[0].m_pfN1[nCurrTri]=afValues[nInts+3];
				  m_akSetsTri[0].m_pfN2[nCurrTri]=afValues[nInts+4];
				  m_akSetsTri[0].m_pfN3[nCurrTri]=afValues[nInts+5];
				  m_akSetsTri[0].m_pfId[nCurrTri]=afValues[3];
				  ++nCurrTri;
			  break;
			  case Type_tetra:
				  m_akSetsTetra[0].m_pfN1[nCurrTet]=afValues[nInts+3];
				  m_akSetsTetra[0].m_pfN2[nCurrTet]=afValues[nInts+4];
				  m_akSetsTetra[0].m_pfN3[nCurrTet]=afValues[nInts+5];
				  m_akSetsTetra[0].m_pfN4[nCurrTet]=afValues[nInts+6];
				  m_akSetsTetra[0].m_pfId[nCurrTet]=afValues[3];
				  ++nCurrTet;
			  break;
		  }

          break;
        }

        default:
          throw std::string("msh_elements_section::read_structure: eState contains an undefined state");

      } /* switch (eState) */

      nPrevStreamPos=rkMshFileStream.tellg();
    } /* while ((eState!=SS_out_elements_section) && getline_checked(rkMshFileStream,kLine)) */
  }
  /****************************************************************************/
  mxArray*
  msh_elements_section::scan(ifstream &rkMshFileStream)
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
