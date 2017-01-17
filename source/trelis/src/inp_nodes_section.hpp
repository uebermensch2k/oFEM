/**
 * inp_nodes_section.hpp
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

#ifndef __INP_NODES_SECTION_HEADER_FILE_INCLUDED__
#define __INP_NODES_SECTION_HEADER_FILE_INCLUDED__

#include <fstream>
#include <string>
#include <vector>

#include "mex.h"

namespace inp_n {

  class inp_nodes_section {
    struct nodes_2d {
      size_t      m_nSectionIdx;
      std::string m_kSectionName;
      size_t      m_nLength     ;

      double     *m_pfBuffer    ;
      double     *m_pfX         ;
      double     *m_pfY         ;

      nodes_2d()
      :
      m_nSectionIdx(0),
      m_kSectionName(""),
      m_nLength(0),
      m_pfBuffer(NULL),
      m_pfX(NULL),
      m_pfY(NULL)
      {
        /**/
      }
    }; /* struct coordinates_2d */

    struct nodes_3d {
      size_t      m_nSectionIdx;
      std::string m_kSectionName;
      size_t      m_nLength     ;

      double     *m_pfBuffer    ;
      double     *m_pfX         ;
      double     *m_pfY         ;
      double     *m_pfZ         ;

      nodes_3d()
      :
      m_nSectionIdx(0),
      m_kSectionName(""),
      m_nLength(0),
      m_pfBuffer(NULL),
      m_pfX(NULL),
      m_pfY(NULL),
      m_pfZ(NULL)
      {
        /**/
      }
    }; /* struct coordinates_3d */

    std::vector<nodes_2d> m_akSets2D;
    std::vector<nodes_3d> m_akSets3D;

//    std::streampos m_nCoordinateSectionStart;


    mxArray *create_and_fill_matlab_structure();

    /* scanner */
  public:
    static const std::string ms_SectionHeader;

  private:
    enum eScanState {
      SS_in_nodes_section,
      SS_read_set_header,
      SS_read_node_2d,
      SS_read_node_3d,
      SS_out_nodes_section
    }; /* enum eScanState */

    void printState     (const eScanState &eState, const std::string &crkLine) const;
    void read_set_header(const std::string &crkLine, std::string &rkSectionName) const;
    void read_node_2d   (const std::string &crkLine, double &rfX, double &rfY) const;
    void read_node_3d   (const std::string &crkLine, double &rfX, double &rfY, double &rfZ) const;

    bool read_structure(std::ifstream &rkInpFileStream);
    void read_contents (std::ifstream &rkInpFileStream);

  public:
    inp_nodes_section();

    mxArray* scan(std::ifstream &rkInpFileStream);
  }; /* class inp_nodes_section */

}; /* namespace inp_n */

#endif /* __INP_NODES_SECTION_HEADER_FILE_INCLUDED__ */
