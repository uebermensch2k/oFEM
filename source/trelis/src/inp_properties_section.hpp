/**
 * inp_properties_section.hpp
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

#ifndef __INP_PROPERTIES_SECTION_HEADER_FILE_INCLUDED__
#define __INP_PROPERTIES_SECTION_HEADER_FILE_INCLUDED__

#include <fstream>
#include <string>
#include <vector>

#include "mex.h"

namespace inp_n {

  class inp_properties_section {

    struct property {
      std::string m_kPropertyName;
      std::string m_kSetName;
      std::string m_kMaterialName;

      property()
      :
      m_kPropertyName(""),
      m_kSetName(""),
      m_kMaterialName("")
      {
        /**/
      }
    }; /* struct property */

    std::vector<property> m_akProperties;


    mxArray *create_and_fill_matlab_structure();

    /* scanner */
  public:
    static const std::string ms_SectionHeader;

  private:
    enum eScanState {
      SS_in_properties_section,
      SS_read_solid_section_header,
      SS_read_shell_section_header,
      SS_out_properties_section
    }; /* enum eScanState */

    void printState               (const eScanState &eState, const std::string &crkLine) const;
    void read_solid_section_header(const std::string &crkLine, std::string &rkSetName, std::string &rkMaterialName) const;
    void read_shell_section_header(const std::string &crkLine, std::string &rkSetName, std::string &rkMaterialName, std::string &rkIntegrationName) const;

    bool read_contents (std::ifstream &rkInpFileStream);

  public:
    inp_properties_section();

    mxArray* scan(std::ifstream &rkInpFileStream);
  }; /* class inp_properties_section */

}; /* namespace inp_n */

#endif /* __INP_PROPERTIES_SECTION_HEADER_FILE_INCLUDED__ */
