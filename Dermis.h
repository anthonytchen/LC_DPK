/* The header file for dermis  */
#ifndef _H_DERMIS_
#define _H_DERMIS_

#include "Comp.h"

class Dermis : public Comp
{
 public:
  // double m_mass_in, m_mass_out; // mass transferred in and out of VE
  bool m_bToBlood;

  // blood involved: 
  double m_dermis_totalV, // total volume of dermis
    m_bld_skin_flow, // total blood flow rate in skin
    m_bld_concChem, // current concentration in blood
    m_bld_fu, // fraction of unbounded solute in blood
    m_par_de2blood, // partition coefficient from dermis to blood
    m_mass_into_dermis, m_mass_outof_dermis; // totoal mass (or molar) flow rates

public:
  Dermis(void) {};	
  ~Dermis(void) {};
  void Init(double, double, double, int, int, bool, CoordSys, BdyCond, BdyCond, BdyCond, BdyCond);
  void InitDermisBlood(double, double, double, double bld_concChem=0, double skin_area=1.8);
  
  void createGrids(Chemical, double);
  void updateBlood(double);
	
  // Functions needed for computing ODE's right hand side (i.e. dy/dt)
  void compODE_dydt (double, const double[], double []);
  static void* static_compODE_dydt_block_threads(void *);
  void compODE_dydt_block (double, const double[], double [], int, int, int, int);
  
  // I/O functions
  void saveCoord(const char [], const char []);
};

#endif
