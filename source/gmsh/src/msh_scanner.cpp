/**
 * @file inp_read_coords.cpp
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "mex.h"

using namespace std;

namespace inp_n {
  
  namespace coordinates {
    
    /**
     * @enum eScanState
     *
     * States of the scanner. It is assumed that we are doing the following job:
     *
     * The scanner initializes in the state SS_searching. If this nasty
     * "***...NODES...*" thing is found the scanner switches to the
     * SS_in_nodes_section state. In the SS_in_nodes_section state the scanner
     * assumes every line to be composed of four comma seperated number. The first
     * one being an interger counting the coordinate number, the following three
     * numbers being floating points representing the 3D coordinates of the node.
     * If being in the SS_in_nodes_section state and encountering the line "**" the
     * scanner switches to the SS_out_nodes_section state. Since we assume only one
     * nodes section, this indicates the end of the scan procedure.
     *
     * TODO: what about 2D meshes? Any difference in reading coordinates?
     */
    enum eScanState {
      SS_in_nodes_section,
      SS_read_nodeset_header,
      SS_read_coordinate,
      SS_out_nodes_section
    };
    
    static const string g_ckNodesHeader   = "********************************** N O D E S **********************************";
    static const string g_ckFunctionIdStr = "MATLAB:mexfunction:inp_read_coords";
    
    /*==============================================================================
     *
     *============================================================================*/
    inline
    void
    PrintErrAndExit(const char * const cpccErrMsg)
    {
      mexErrMsgIdAndTxt(g_ckFunctionIdStr.c_str(), cpccErrMsg);
    }
    
    
    /*==============================================================================
     *
     *============================================================================*/
    void
    open_inp_file(const mxArray *prhs,
                  ifstream &kInpFileStream)
    {
      /* read file name */
      char *pcFileName = mxArrayToString(prhs);
      if (!pcFileName)
        PrintErrAndExit("I can't create the file name string. I'm out of memory.\n");
      
      string kFileName(pcFileName);
      mxFree(pcFileName);
      
      /* check for *.inp */
      if (kFileName.find(".inp")==string::npos)
        kFileName.append(".inp");
      
      /* open file */
      kInpFileStream.open(kFileName.c_str());
      
      if (!kInpFileStream.is_open())
        PrintErrAndExit("I can't open the specified file. Are you sure it exists?\n");
    }
    
    
    /*==============================================================================
     *
     *============================================================================*/
    size_t
    get_coordinate_dimension(const string &crkLine)
    {
      std::istringstream kLineStream(crkLine);
      
      // parse line
      size_t nDim=0;
      double fValue=0.;
      bool bEOF=false;
      while (!bEOF)
      {
        // read until numerical value is encountered
        while (true)
        {
          kLineStream>>fValue;
          if (kLineStream.bad())
            PrintErrAndExit("Bad stream state.\n");
          
          if (kLineStream.fail())
          {
            kLineStream.clear();
            kLineStream.ignore();
            if (kLineStream.eof())
            {
              bEOF=true;
              break;
            }
            continue;
          }
          break;
        }
        
        if (!bEOF) ++nDim;
        
        bEOF=kLineStream.eof();
      }
      
      return nDim;
    }
    
    
    /*==============================================================================
     *
     *============================================================================*/
    void
    append_coordinate(const string &crkLine   ,
                      const size_t &crnCurrRow,
                      double       *pfXCoords ,
                      double       *pfYCoords ,
                      double       *pfZCoords )
    {
      int nNodesNumber;
      double fCoord;
      string kCoord;
      istringstream kSStream;
      
      size_t nPosB = 0;
      size_t nPosE = crkLine.find(',');
      
      /* read x-coordinate */
      nPosB = nPosE+1;
      nPosE = crkLine.find(',',nPosB);
      kSStream.str(crkLine.substr(nPosB,nPosE-nPosB));
      kSStream>>fCoord;
      pfXCoords[crnCurrRow] = fCoord;
      
      /* read y-coordinate */
      nPosB = nPosE+1;
      nPosE = crkLine.find(',',nPosB);
      kSStream.clear();
      kSStream.str(crkLine.substr(nPosB,nPosE-nPosB));
      kSStream>>fCoord;
      pfYCoords[crnCurrRow] = fCoord;
      
      /* read z-coordinate */
      kSStream.clear();
      kSStream.str(crkLine.substr(nPosE+1));
      kSStream>>fCoord;
      pfZCoords[crnCurrRow] = fCoord;
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
      open_inp_file(prhs[0], kInpFileStream);
      /****************************************************************************/
      /****************************************************************************/
      
      
      
      
      
      /* number of nodes to be read */
      size_t nNodesQuantity=0;
      
      
      
      
      
      /****************************************************************************/
      /* first scan, here the number of nodes is read                             */
      /****************************************************************************/
      {
        eScanState eState=SS_searching;
        bool bDoScanning=true;
        streampos nNodeSectionPos=0;
        string kLine;
        
        while (bDoScanning && (eState!=SS_out_nodes_section))
        {
          getline(kInpFileStream, kLine);
          
          /* check stream state */
          if (!kInpFileStream.good())
          {
            if (kInpFileStream.eof())
            {
              bDoScanning=false;
              continue;
            }
            else
            {
              /* everything else is an hard error => return */
              PrintErrAndExit("Bad stream state.\n");
            }
          }
          
          /* remove possible carriage return */
          if (*(kLine.end()-1)=='\r') kLine.erase(kLine.end()-1,kLine.end());
          
          /* check for scanner state */
          switch (eState)
          {
            case SS_searching:
            {
              //          mexPrintf("state: searching.\n");
              //          mexPrintf("state: SS_searching.\n");
              if (!kLine.compare(ckNodesHeader))
              {
                /* now pointing to the first line after the ***...NODES...*** line */
                nNodeSectionPos=kInpFileStream.tellg();
                eState = SS_in_nodes_section;
              }
              
              continue;
            }
              
            case SS_in_nodes_section:
            {
              //          mexPrintf("state: SS_in_node_section.\n");
              //          mexPrintf("state: in nodes section.\n");
              /* section ends? */
              if (!kLine.compare("**"))
              {
                //            mexPrintf("state: SS_out_node_section.\n");
                /* *NODE line counted => correct */
                --nNodesQuantity;
                
                eState = SS_out_nodes_section;
                continue;
              }
              
              /* read next coordinate */
              ++nNodesQuantity;
              continue;
            }
              
              /*case SS_out_nodes_section:
               {
               bDoScanning = false;
               continue;
               }*/
              
            default:
              PrintErrAndExit("Something went badly wrong during first scan.\n");
              
          } // switch (eState)
        } // while (bDoScanning)
        
        if (eState == SS_searching)
          PrintErrAndExit("I haven't found a NODES section in the file.\n");
        
        if (eState != SS_out_nodes_section)
          PrintErrAndExit("Format of inp-file is corrupted. I haven't found the end of the NODES section in the file.\n");
        
        
        /* reset stream to point to begin of the NODES section, i.e. the first line
         * after ***...NODES...*** */
        kInpFileStream.clear();
        kInpFileStream.seekg(nNodeSectionPos);
      } // first scan
      /****************************************************************************/
      /****************************************************************************/
      
      
      
      
      /****************************************************************************/
      /* create buffers to store data                                             */
      /****************************************************************************/
      double *pfMatrixData = (double*) mxCalloc(3*nNodesQuantity, sizeof(double));
      
      if (!pfMatrixData)
        PrintErrAndExit("I was unable to allocate memory for the coordinate buffer.\n");
      
      double *pfXCoords=pfMatrixData                 ;
      double *pfYCoords=pfMatrixData+  nNodesQuantity;
      double *pfZCoords=pfMatrixData+2*nNodesQuantity;
      /****************************************************************************/
      /****************************************************************************/
      
      
      
      
      /****************************************************************************/
      /* now scan the contents                                                    */
      /****************************************************************************/
      {
        /* stream points into NODES section */
        eScanState eState=SS_in_nodes_section;
        bool bDoScanning=true;
        size_t nCurrRow=0;
        string kLine;
        
        /* fisr line is *NODE line */
        getline(kInpFileStream, kLine);
        
        while (bDoScanning && (eState!=SS_out_nodes_section))
        {
          getline(kInpFileStream, kLine);
          
          /* check for scanner state */
          if (!kInpFileStream.good())
          {
            if (kInpFileStream.eof())
            {
              bDoScanning=false;
              continue;
            }
            else
            {
              /* everything else is an hard error => return */
              mxFree(pfMatrixData);
              PrintErrAndExit("Bad stream state.\n");
            }
          }
          
          /* remove possible carriage return */
          if (*(kLine.end()-1)=='\r') kLine.erase(kLine.end()-1,kLine.end());
          
          /* check contents of line and decide scanners state */
          switch (eState)
          {
            case SS_in_nodes_section:
            {
              //          mexPrintf("state: in nodes section.\n");
              
              /* section end ? */
              if (!kLine.compare("**"))
              {
                eState = SS_out_nodes_section;
                continue;
              }
              
              /* read next coordinate */
              append_coordinate(kLine,nCurrRow,pfXCoords,pfYCoords,pfZCoords);
              ++nCurrRow;
              
              continue;
            }
              
              /*case SS_out_nodes_section:
               {
               bDoScanning = false;
               continue;
               }*/
              
            default:
              mxFree(pfMatrixData);
              PrintErrAndExit("Something went badly wrong in main scan loop.\n");
              
          } // switch (eState)
        } // while (bDoScanning)
      } // scan contents
      /****************************************************************************/
      /****************************************************************************/
      
      
      
      
      /****************************************************************************/
      /* nodes read successfully => create output */
      /****************************************************************************/
      plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
      
      if (!plhs[0])
      {
        mxFree(pfMatrixData);
        PrintErrAndExit("I was unable to create the coordinate matrix.\n");
      }
      
      mxSetData(plhs[0], pfMatrixData  );
      mxSetM   (plhs[0], nNodesQuantity);
      mxSetN   (plhs[0], 3             );
      /****************************************************************************/
      /****************************************************************************/
    }
    
  }
  
} /* namespace inp_n */
