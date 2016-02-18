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
  void Init(double, double, double, int, int, double, double, double, CoordSys, BdyCond, BdyCond, BdyCond, BdyCond);
  
  void createGrids(Chemical, double, double);
 
  // I/O functions
  void saveCoord(const char [], const char []);
};

#endif
