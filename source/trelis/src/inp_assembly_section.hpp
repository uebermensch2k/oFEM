/**
 * inp_assembly_section.hpp
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

#ifndef __INP_ASSEMBLY_SECTION_HEADER_FILE_INCLUDED__
#define __INP_ASSEMBLY_SECTION_HEADER_FILE_INCLUDED__

#include <fstream>
#include <string>
#include <vector>

#include "mex.h"

namespace inp_n {

  class inp_assembly_section {
    struct material {
      struct property {
        std::string         m_kName;
        std::vector<double> m_akValues;

        property()
        :
        m_kName(""),
        m_akValues(0)
        {
          /**/
        }
      }; /* property */

      std::string           m_kName;
      std::vector<property> m_akProperties;

      material()
      :
      m_kName(""),
      m_akProperties(0)
      {
        /**/
      }
    }; /* struct material */

    std::vector<material> m_akMaterials;

    mxArray *create_and_fill_matlab_structure();

    /* scanner */
  public:
    static const std::string ms_SectionHeader;

  private:
    enum eScanState {
      SS_in_assembly_section,
      SS_read_material_header,
      SS_read_property_header,
      SS_read_property,
      SS_out_assembly_section
    }; /* enum eScanState */

    void printState(const eScanState &eState,const std::string &crkLine) const;

    void read_material_header(const std::string &crkLine, std::string &rkSectionName) const;
    void read_property_header(const std::string &crkLine, std::string &rkSectionName) const;
    void read_property       (const std::string &crkLine, std::vector<double> &rafValues) const;

    bool read_structure(std::ifstream &rkInpFileStream);
//    void read_contents (std::ifstream &rkInpFileStream);

  public:
    inp_assembly_section();

    mxArray* scan(std::ifstream &rkInpFileStream);
  }; /* class inp_assembly_section */

}; /* namespace inp_n */

#endif /* __INP_ASSEMBLY_SECTION_HEADER_FILE_INCLUDED__ */
