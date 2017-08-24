/**
 * inp_elements_section.hpp
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

#ifndef __INP_ELEMENTS_SECTION_HEADER_FILE_INCLUDED__
#define __INP_ELEMENTS_SECTION_HEADER_FILE_INCLUDED__

#include <fstream>
#include <string>
#include <vector>

#include "mex.h"

namespace inp_n {

  class inp_elements_section {
    struct elements_tri {
      size_t      m_nSectionIdx;
      std::string m_kSectionName;
      size_t      m_nLength     ;

      double     *m_pfBuffer    ;
      double     *m_pfN1        ;
      double     *m_pfN2        ;
      double     *m_pfN3        ;
      double     *m_pfId        ;

      elements_tri()
      :
      m_nSectionIdx(0),
      m_kSectionName(""),
      m_nLength(0),
      m_pfBuffer(NULL),
      m_pfN1(NULL),
      m_pfN2(NULL),
      m_pfN3(NULL),
      m_pfId(NULL)
      {
        /**/
      }
    }; /* struct elements_tri */

    struct elements_tetra {
      size_t      m_nSectionIdx;
      std::string m_kSectionName;
      size_t      m_nLength     ;

      double     *m_pfBuffer    ;
      double     *m_pfN1        ;
      double     *m_pfN2        ;
      double     *m_pfN3        ;
      double     *m_pfN4        ;
      double     *m_pfId        ;

      elements_tetra()
      :
      m_nSectionIdx(0),
      m_kSectionName(""),
      m_nLength(0),
      m_pfBuffer(NULL),
      m_pfN1(NULL),
      m_pfN2(NULL),
      m_pfN3(NULL),
      m_pfN4(NULL),
      m_pfId(NULL)
      {
        /**/
      }
    }; /* struct elements_tetra */

    std::vector<elements_tri  > m_akSetsTri  ;
    std::vector<elements_tetra> m_akSetsTetra;


    mxArray *create_and_fill_matlab_structure();

    /* scanner */
  public:
    static const std::string ms_SectionHeader;

  private:
    enum eScanState {
      SS_in_elements_section,
      SS_read_set_header,
      SS_read_element_tri,
      SS_read_element_tetra,
      SS_out_elements_section
    }; /* enum eScanState */

    void read_set_header   (const std::string &crkLine, std::string &rkSectionName, std::string &rkSectionType) const;
    void read_element_tri  (const std::string &crkLine, double &rfN1, double &rfN2, double &rfN3, double &rfID) const;
    void read_element_tetra(const std::string &crkLine, double &rfN1, double &rfN2, double &rfN3, double &rfN4, double &rfID) const;

    bool read_structure(std::ifstream &rkInpFileStream);
    void read_contents (std::ifstream &rkInpFileStream);

  public:
    inp_elements_section();

    mxArray* scan(std::ifstream &rkInpFileStream);
  }; /* class inp_elements_section */

}; /* namespace inp_n */

#endif /* __INP_ELEMENTS_SECTION_HEADER_FILE_INCLUDED__ */
