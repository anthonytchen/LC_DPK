/* The header file for sebum  */
#ifndef _H_SEBUM_
#define _H_SEBUM_

#include "Comp.h"

class Sebum : public Comp
{
 public:
  double m_T, m_eta, // temperature, sebum viscosity at the temperature
    m_K_sw,  // partition between sebum and water
    m_D, // diffusion coefficient in sebum
    m_init_concChem; // initial sebum concentration

public:
  Sebum(void) {};	
  ~Sebum(void) {};
  void Init(double, double, double, int, int, CoordSys, BdyCond, BdyCond, BdyCond, BdyCond, double K_sw=-1, double D_sebum=-1);
  
  void createGrids(Chemical, double, double);
 
  // I/O functions
  void saveCoord(const char [], const char []);
};

#endif
