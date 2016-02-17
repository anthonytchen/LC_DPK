/* The header file for stratum corneum  */
#ifndef _H_STRACORN_
#define _H_STRACORN_

#include "Comp.h"

class StraCorn : public Comp
{
 public:
  double m_rou_lipid, m_rou_keratin, m_rou_water, // the density of lipid, keratin and water
    m_T, m_eta, // temperature, viscosity of water
    m_V_mortar, m_V_brick, m_V_all; // the volume of mortar, brick and sum (all) 
                                    //	in each element of stratum corneum dimensions
  int m_nx_grids_lipid, m_nx_grids_cc, // number of grids for each lipid/corneocyte layer in the x direction
    m_ny_grids_lipid, m_ny_grids_cc_dn, // number of grids for each lipid/corneocyte layer in the y direction
    m_n_layer_x; // number of layers in the x direction
  double m_geom_g, m_geom_d, m_geom_s, m_geom_t, m_geom_dm, m_geom_dn, m_w, // geometry parameters
    m_offset_y; // offset at the left simulation boundary, relative to corneocyte
  //  double m_mass_in, m_mass_out; // the mass transferred in and out of stratum corneum

public:
  StraCorn(void) {};	
  ~StraCorn(void) {};
  void Init(double, double, double, double, double, int, int, double, CoordSys, BdyCond, BdyCond, BdyCond, BdyCond);
  void Release();
  
  void createGrids(Chemical, double, double, double);
 
  // I/O functions
  void getAmount(double*, double*, double*);
  void comp1DConc();
  void saveCoord(const char [], const char []);

};

#endif
