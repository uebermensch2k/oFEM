/**
 * msh_elements_section.hpp
 *
 * Copyright: Michael Hagel,
 *  Department of the Theory of Electrical Engineering,
 *  University of the Federal Armed Forces Hamburg,
 *  Hamburg,
 *  Germany
 */

/**
 * \file
 * \author Michael Hagel
 * \date 15.05.2108
 */

#ifndef __MSH_PHYSICAL_SECTION_HEADER_FILE_INCLUDED__
#define __MSH_PHYSICAL_SECTION_HEADER_FILE_INCLUDED__

#include <fstream>
#include <string>
#include <vector>

#include "mex.h"

namespace msh_n {

  class msh_physical_section {

	struct physical{
		int 		m_nPhysicalIdx;
		std::string m_kPhysicalName;

		physical()
		:
		m_nPhysicalIdx(0),
		m_kPhysicalName("")
		{
		  /**/
		}
	};

	std::vector<physical > m_akSetsPhysical ;

    mxArray *create_and_fill_matlab_structure();

    /* scanner */
  public:
    static const std::string ms_SectionHeader;

  private:
    enum eScanState {
      SS_in_physical_section,
      SS_read_set_header,
	  SS_read_physical,
      SS_out_physical_section
    }; /* enum eScanState */

	void read_physical (std::string &crkLine, std::string &physicalName, int &Id);
    bool read_structure(std::ifstream &rkMshFileStream);
    void read_contents (std::ifstream &rkMshFileStream);

  public:
    msh_physical_section();

    mxArray* scan(std::ifstream &rkMshFileStream);
  }; /* class msh_elements_section */

}; /* namespace msh_n */

#endif /* __MSH_PHYSICAL_SECTION_HEADER_FILE_INCLUDED__ */
