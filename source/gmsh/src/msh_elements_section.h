/**
 * msh_elements_section.h
 *
 * (based on inp_elements_section.hpp by Michael Dudzinski)
 *
 * Copyright: Michael Dudzinski,
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
 * \date 14.12.2017
 */

#ifndef __MSH_ELEMENTS_SECTION_H__
#define __MSH_ELEMENTS_SECTION_H__

#include <stdio.h>
#include <fstream>
#include <string.h>
#include <vector>
#include <unordered_set>

#include "mex.h"
#include "msh_util.h"

namespace msh_n {

  class msh_elements_section {
    
    struct elements_line {
      size_t      m_nSectionIdx;
      std::string m_kSectionName;
      size_t      m_nLength     ;

      double     *m_pfBuffer    ;
      double     *m_pfN1        ;
      double     *m_pfN2        ;
      double     *m_pfId        ;

      elements_line()
      :
      m_nSectionIdx(0),
      m_kSectionName(""),
      m_nLength(0),
      m_pfBuffer(NULL),
      m_pfN1(NULL),
      m_pfN2(NULL),
      m_pfId(NULL)
      {
        /**/
      }
    }; // struct elements_line

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
    }; // struct elements_tri

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
    }; // struct elements_tetra

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
    
    std::vector<elements_line> m_akSetsLine;
    std::vector<elements_tri> m_akSetsTri;
    std::vector<elements_tetra> m_akSetsTetra;
    std::vector<surface> m_akSurfaces;

    mxArray *create_and_fill_matlab_structure();
    mxArray *create_and_fill_matlab_structure_bd();

  public:
    static const std::string ms_DimHeader;
    static const std::string ms_PhysicalEntitiesHeader;
    static const std::string ms_PhysicalEntitiesEnd;
    static const std::string ms_ElementsHeader;
    static const std::string ms_ElementsEnd;

  private:
    enum eScanState {
      SS_read_header_dim,
      SS_read_nb_dim,
      SS_read_header_entities,
      SS_read_nb_entities,
      SS_read_entity,
      SS_read_header_elements,
      SS_read_nb_elements,
      SS_read_element,
      SS_out_elements_section
    };

    size_t nGeometryDim;
    std::vector<size_t> m_anPartMapping;
    std::vector<size_t> m_anBdPartMapping;
    
    inline int findIndex(const std::vector<size_t> &v, const size_t nId);
    void findElementByLine(const size_t nLineN1, const size_t nLineN2, std::vector<double> &rvnE1, std::vector<double> &rvnE2, std::vector<double> &rvnE3);
    
    void read_entity(const std::string &crkLine, size_t &rnId, size_t &rnDim, std::string &rkName) const;
    void read_element(const std::string &crkLine, size_t &rType, size_t &rN1, size_t &rN2, size_t &rN3, size_t &rN4, size_t &rEntity, size_t &rID) const;
    
    bool read_structure(std::ifstream &rkMshFileStream);
    void read_contents (std::ifstream &rkMshFileStream);
    void find_sidesets();

  public:
    msh_elements_section();

    mxArray* scan(std::ifstream &rkMshFileStream);
    mxArray* scan_bd();
  }; // class msh_elements_section

}; // namespace msh_n

#endif // __MSH_ELEMENTS_SECTION_H__
