/**
 * @file inp_read_nodesets.cpp
 *
 * inp_read_coords takes the name of an inp-file and returns the
 * coordinates of the nodes of the contained mesh.
 *
 * Copyright: Michael Dudzinski,
 *  Department of the Theory of Electrical Engineering,
 *  University of the Federal Armed Forces Hamburg,
 *  Hamburg,
 *  Germany
 */



//#include "matrix.h"

#include "inp_nodes_section.hpp"
#include "inp_elements_section.hpp"
#include "inp_nodesets_section.hpp"
#include "inp_sidesets_section.hpp"
#include "inp_properties_section.hpp"
#include "inp_assembly_section.hpp"
#include "inp_utils.hpp"

using namespace std;
using namespace inp_n;

const char * const g_cpccFunctionIdStr = "MATLAB:mexfunction:inp_read_file";

/*==============================================================================
 *
 *============================================================================*/
inline
void
PrintErrAndExit(const char * const cpccErrMsg)
{
  mexErrMsgIdAndTxt(g_cpccFunctionIdStr, cpccErrMsg);
}



/*==============================================================================
 *
 *============================================================================*/
mxArray* create_cell_and_fill_partially()
{
  mwSize anDims[2] = { 2, 6 };
  mxArray *pkCell  = mxCreateCellArray(2,anDims);

  // nodes section
  mxArray *pkNameMatrix = mxCreateString("coords");
  mxSetCell(pkCell,0*2,pkNameMatrix);

  // elements section
  pkNameMatrix = mxCreateString("elems");
  mxSetCell(pkCell,1*2,pkNameMatrix);

  // nodesets section
  pkNameMatrix = mxCreateString("nodesets");
  mxSetCell(pkCell,2*2,pkNameMatrix);

  // sidesets section
  pkNameMatrix = mxCreateString("sidesets");
  mxSetCell(pkCell,3*2,pkNameMatrix);

  // properties section
  pkNameMatrix = mxCreateString("properties");
  mxSetCell(pkCell,4*2,pkNameMatrix);

  // materials section
  pkNameMatrix = mxCreateString("materials");
  mxSetCell(pkCell,5*2,pkNameMatrix);

  return pkCell;
}


/*==============================================================================
 *
 *============================================================================*/
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
//  mexPrintf("Entering mexFunction\n");

  /****************************************************************************/
  /* check validity of arguments                                              */
  /****************************************************************************/
  /* check number ... */
  if (nrhs>1)
    PrintErrAndExit("Exactly one input argument expected.\n");
  
  /* ... and type of the input */
  if ( !mxIsChar(prhs[0]) || (mxGetM(prhs[0])!=1) )
    PrintErrAndExit("Input argument must be a string.\n");
  /****************************************************************************/
  /****************************************************************************/





  /****************************************************************************/
  /* open inp-file                                                            */
  /****************************************************************************/
  ifstream kInpFileStream;
  try {
    open_inp_file(prhs[0], kInpFileStream);
  }
  catch (string err) {
    PrintErrAndExit((err+"\n").c_str());
  }
  /****************************************************************************/
  /****************************************************************************/

  
  
  

  /****************************************************************************/
  /* scan section                                                             */
  /****************************************************************************/
  inp_nodes_section      nodes_section     ;
  inp_elements_section   elements_section  ;
  inp_nodesets_section   nodesets_section  ;
  inp_sidesets_section   sidesets_section  ;
  inp_properties_section properties_section;
  inp_assembly_section   assembly_section;
  string kLine;

  
  mxArray *pkCell = create_cell_and_fill_partially();

  try {
//	mexPrintf("Entering scan loop\n");
    while (getline_checked(kInpFileStream, kLine))
    {
//      mexPrintf("line  : %s\n",kLine.c_str());
//      continue;
//      mexPrintf("header: %s.\n",inp_nodes_section::ms_SectionHeader.c_str());
      
      if (!kLine.compare(inp_nodes_section::ms_SectionHeader))
      {
        mxArray *pkArray = nodes_section.scan(kInpFileStream);
//        mxArray *pkArray = mxCreateCellArray(0, NULL);
        mxSetCell(pkCell,0*2+1,pkArray);
      }
      else if (!kLine.compare(inp_elements_section::ms_SectionHeader))
      {
        mxArray *pkArray = elements_section.scan(kInpFileStream);
//        mxArray *pkArray = mxCreateCellArray(0, NULL);
        mxSetCell(pkCell,1*2+1,pkArray);
      }
      else if (!kLine.compare(inp_nodesets_section::ms_SectionHeader))
      {
        mxArray *pkArray = nodesets_section.scan(kInpFileStream);
//        mxArray *pkArray = mxCreateCellArray(0, NULL);
        mxSetCell(pkCell,2*2+1,pkArray);
      }
      else if (!kLine.compare(inp_sidesets_section::ms_SectionHeader))
      {
        mxArray *pkArray = sidesets_section.scan(kInpFileStream);
//        mxArray *pkArray = mxCreateCellArray(0, NULL);
        mxSetCell(pkCell,3*2+1,pkArray);
      }
      else if (!kLine.compare(inp_properties_section::ms_SectionHeader))
      {
        mxArray *pkArray = properties_section.scan(kInpFileStream);
//        mxArray *pkArray = mxCreateCellArray(0, NULL);
        mxSetCell(pkCell,4*2+1,pkArray);
      }
      else if (!kLine.compare(inp_assembly_section::ms_SectionHeader))
      {
        mxArray *pkArray = assembly_section.scan(kInpFileStream);
//        mxArray *pkArray = mxCreateCellArray(0, NULL);
        mxSetCell(pkCell,5*2+1,pkArray);
      }

    } /* while (getline_checked(kInpFileStream, kLine)) */
  }
  catch (string err) {
    PrintErrAndExit((err+"\n").c_str());
  }

  
  
  
  
  /****************************************************************************/
  /* nothing found                                                            */
  /****************************************************************************/
  plhs[0] = pkCell;
  /****************************************************************************/
  /****************************************************************************/
}
