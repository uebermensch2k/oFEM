/**
 * msh_nodes_section.h
 *
 * (based on inp_nodes_section.hpp by Michael Dudzinski)
 *
 * Copyright:
 *  Michael Dudzinski,
 *  Department of the Theory of Electrical Engineering,
 *  University of the Federal Armed Forces Hamburg,
 *  Hamburg,
 *  Germany
 *
 *  Harald Scharf,
 *  French-German Research Institute of Saint-Louis
 *  Saint-Louis, France
 */

/**
 * \file
 * \author Harald Scharf
 * \date 13.12.2017
 */

#ifndef __MSH_NODES_SECTION_H__
#define __MSH_NODES_SECTION_H__

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include "mex.h"
#include "msh_util.h"

namespace msh_n {
  
  class msh_nodes_section {
    struct nodes_2d {
      size_t      m_nSectionIdx;
      std::string m_kSectionName;
      size_t      m_nLength     ;

      double     *m_pfBuffer    ;
      double     *m_pfX         ;
      double     *m_pfY         ;

      nodes_2d() :
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

      nodes_3d() :
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
    }; // struct coordinates_3d

    std::vector<nodes_2d> m_akSets2D;
    std::vector<nodes_3d> m_akSets3D;

    mxArray *create_and_fill_matlab_structure();

  // scanner
  public:
    static const std::string ms_DimHeader;
    static const std::string ms_NodesHeader;
    static const std::string ms_NodesEnd;

  private:
    enum eScanState {
      SS_read_header_dim,
      SS_read_nb_dim,
      SS_read_header_nodes,
      SS_read_nb_nodes,
      SS_read_node,
      SS_out_nodes_section
    }; // enum eScanState
    
    size_t nGeometryDim;

    void read_node_2d   (const std::string &crkLine, double &rfX, double &rfY) const;
    void read_node_3d   (const std::string &crkLine, double &rfX, double &rfY, double &rfZ) const;
    void read_contents (std::ifstream &rkMshFileStream);

  public:
    msh_nodes_section();

    mxArray* scan(std::ifstream &rkMshFileStream);
  }; // class msh_nodes_section

};  // namespace msh_n

#endif // __MSH_NODES_SECTION_H__
