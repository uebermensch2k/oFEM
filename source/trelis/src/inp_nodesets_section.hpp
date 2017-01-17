/**
 * inp_nodesets_section.hpp
 *
 * Copyright: Michael Dudzinski,
 *  Department of the Theory of Electrical Engineering,
 *  University of the Federal Armed Forces Hamburg,
 *  Hamburg,
 *  Germany
 */

/**
 * \file
 * \author Michael Dudzinski
 * \date 10.09.2014
 */

#ifndef __INP_NODESETS_SECTION_HEADER_FILE_INCLUDED__
#define __INP_NODESETS_SECTION_HEADER_FILE_INCLUDED__

#include <fstream>
#include <string>
#include <vector>

#include "mex.h"

namespace inp_n {

  class inp_nodesets_section {
//    struct nodeset {
//      std::string m_kSectionName;
//
//      std::vector<size_t> m_kNodeIDs;
//
//      size_t      m_nLength     ;
//      double     *m_pfBuffer    ;
//
//      nodeset()
//      :
//      m_kSectionName(""),
//      m_kNodeIDs(0),
//      m_nLength(0),
//      m_pfBuffer(NULL)
//      {
//        /**/
//      }
//    }; /* struct nodeset */
    struct nodeset {
      std::string           m_kSectionName;
      std::vector<uint32_t> m_kNodeIDs;
      
      nodeset()
      :
      m_kSectionName(""),
      m_kNodeIDs(0)
      {
        /**/
      }
    }; /* struct nodeset */

    std::vector<nodeset> m_akNodesets;


    mxArray *create_and_fill_matlab_structure();

    /* scanner */
  public:
    static const std::string ms_SectionHeader;

  private:
    enum eScanState {
      SS_in_nodesets_section,
      SS_read_set_header,
      SS_read_node_ids,
      SS_out_nodesets_section
    }; /* enum eScanState */

    void printState      (const eScanState &eState, const std::string &crkLine) const;
    void read_set_header (const std::string &crkLine, std::string &rkSectionName) const;
    void read_node_ids   (const std::string &crkLine, std::vector<uint32_t> &rkNodeIDs);

    void read_contents (std::ifstream &rkInpFileStream);

  public:
    inp_nodesets_section();

    mxArray* scan(std::ifstream &rkInpFileStream);
  }; /* class inp_nodesets_section */

}; /* namespace inp_n */

#endif /* __INP_NODESETS_SECTION_HEADER_FILE_INCLUDED__ */
