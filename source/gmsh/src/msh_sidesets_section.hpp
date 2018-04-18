/**
 * inp_sidesets_section.hpp
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

#ifndef __INP_SIDESETS_SECTION_HEADER_FILE_INCLUDED__
#define __INP_SIDESETS_SECTION_HEADER_FILE_INCLUDED__

#include <fstream>
#include <string>
#include <vector>

#include "mex.h"

namespace inp_n {

  class inp_sidesets_section {

    struct sideset {
      enum type_e {
        curve   = 0,
        surface = 1,
        edge    = 2,
        face    = 3,
        tri     = 4,
        invalid
      }; /* enum type_e */
      
      type_e      m_eType;

      std::string m_kSectionName;

      size_t      m_nID;

      size_t      m_nLength     ;
      double     *m_pfBuffer    ;

      sideset()
      :
      m_eType(invalid),
      m_kSectionName(""),
      m_nID(0),
      m_nLength(0),
      m_pfBuffer(NULL)
      {
        /**/
      }
    }; /* struct sideset */

    struct surface {

      std::string m_kSectionName;

      std::vector<sideset> m_akSidesets;

      surface()
      :
      m_kSectionName(""),
      m_akSidesets(0)
      {
        /**/
      }
    }; /* struct surface */

    std::vector<surface> m_akSurfaces;


    mxArray *create_and_fill_matlab_structure();

    /* scanner */
  public:
    static const std::string ms_SectionHeader;

  private:
    enum eScanState {
      SS_in_sidesets_section,
      SS_read_elementset_header,
      SS_read_element_ids,
      SS_read_surface_header,
      SS_read_elementset_name_and_type,
      SS_out_sidesets_section
    }; /* enum eScanState */

    void read_elementset_header(const std::string &crkLine, std::string &rkSectionName) const;
    void read_element_ids      (const std::string &crkLine, double *pfElementIDs, size_t &rnLastIdx);
    size_t count_element_ids   (const std::string &crkLine);
    void read_surface_header   (const std::string &crkLine, std::string &rkName) const;
    void read_elementset_name_type_and_id(const std::string &crkLine, std::string &rkName, sideset::type_e &reType, size_t &rkID) const;

    bool read_structure(std::ifstream &rkInpFileStream);
    void read_contents (std::ifstream &rkInpFileStream);

  public:
    inp_sidesets_section();

    mxArray* scan(std::ifstream &rkInpFileStream);
  }; /* class inp_sidesets_section */

}; /* namespace inp_n */

#endif /* __INP_SIDESETS_SECTION_HEADER_FILE_INCLUDED__ */
