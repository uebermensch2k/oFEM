/**
 * @file inp_read_asm.cpp
 *
 * inp_read_asm takes the name of an inp-file and returns the
 * materials defined in the assembly section.
 *
 * Copyright: Michael Dudzinski,
 *  Department of the Theory of Electrical Engineering,
 *  University of the Federal Armed Forces Hamburg,
 *  Hamburg,
 *  Germany
 */



//#include "matrix.h"

#include "inp_assembly_section.hpp"
#include "inp_utils.hpp"

using namespace std;
using namespace inp_n;

const char * const g_cpccFunctionIdStr = "MATLAB:mexfunction:inp_read_asm";

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
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
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
  inp_assembly_section assembly_section;
  string kLine;

  while (getline_checked(kInpFileStream, kLine))
  {
//    mexPrintf("line  : %s.\n",kLine.c_str());
//    mexPrintf("header: %s.\n",inp_nodes_section::ms_SectionHeader.c_str());
    
    if (!kLine.compare(inp_assembly_section::ms_SectionHeader))
    {
      try {
//        mexPrintf("stream position: %d.\n",(int) kInpFileStream.tellg());
        plhs[0] = assembly_section.scan(kInpFileStream);
      }
      catch (string err) {
        PrintErrAndExit((err+"\n").c_str());
      }

      return;
    }
  } // while (bDoScanning)

  
  
  
  
  /****************************************************************************/
  /* nothing found                                                            */
  /****************************************************************************/
  plhs[0] = mxCreateCellArray(0, NULL);
  /****************************************************************************/
  /****************************************************************************/
}
