/**
 * @file msh_read_file.cpp
 *
 * msh_read_file takes the name of a msh-file and returns the
 * coordinates, elements and sidesets of the contained mesh.
 *
 * based on inp_read_file.cpp
 *
 * Copyright: Michael Dudzinski,
 *  Department of the Theory of Electrical Engineering,
 *  University of the Federal Armed Forces Hamburg,
 *  Hamburg,
 *  Germany
 *
 * Modified (and renamed) by: Harald Scharf,
 *  French-German Research Institute of Saint-Louis,
 *  Saint-Louis, France
 *  13.12.2017
 */

#include "msh_nodes_section.h"
#include "msh_elements_section.h"

using namespace std;
using namespace msh_n;

const char * const g_cpccFunctionIdStr = "MATLAB:mexfunction:msh_read_file";

/*==============================================================================
 *
 *============================================================================*/
inline void PrintErrAndExit(const char * const cpccErrMsg)
{
  mexErrMsgIdAndTxt(g_cpccFunctionIdStr, cpccErrMsg);
}

/*==============================================================================
 *
 *============================================================================*/
mxArray* create_cell_and_fill_partially()
{
  mwSize anDims[2] = { 2, 3 };
  mxArray *pkCell = mxCreateCellArray(2, anDims);

  // nodes section
  mxArray *pkNameMatrix = mxCreateString("coords");
  mxSetCell(pkCell, 0*2, pkNameMatrix);

  // elements section
  pkNameMatrix = mxCreateString("elems");
  mxSetCell(pkCell, 1*2, pkNameMatrix);

  // sidesets section
  pkNameMatrix = mxCreateString("sidesets");
  mxSetCell(pkCell, 2*2, pkNameMatrix);

  return pkCell;
}

/*==============================================================================
 *
 *============================================================================*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // check validity of arguments
  // check number ...
  if (nrhs>1)
    PrintErrAndExit("Exactly one input argument expected.\n");
  
  // ... and type of the input
  if ( !mxIsChar(prhs[0]) || (mxGetM(prhs[0])!=1) )
    PrintErrAndExit("Input argument must be a string.\n");

  // open msh-file
  ifstream kMshFileStream;
  try {
    open_msh_file(prhs[0], kMshFileStream);
  }
  catch (string err) {
    PrintErrAndExit((err+"\n").c_str());
  }

  // scan section
  msh_nodes_section nodes_section;
  msh_elements_section elements_section;
  string kLine;

  mxArray *pkCell = create_cell_and_fill_partially();

  try {
    mxArray *pkArray = nodes_section.scan(kMshFileStream);
    mxSetCell(pkCell, 0*2+1, pkArray);
    
    kMshFileStream.clear();
    kMshFileStream.seekg(0);
    
    pkArray = elements_section.scan(kMshFileStream);
    mxSetCell(pkCell, 1*2+1, pkArray);
    
    pkArray = elements_section.scan_bd();
    mxSetCell(pkCell, 2*2+1, pkArray);
  }
  catch (string err) {
    PrintErrAndExit((err+"\n").c_str());
  }

  plhs[0] = pkCell;
}
