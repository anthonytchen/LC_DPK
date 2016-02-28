#include "stdafx.h"
#include "SurSebum.h"

void SurSebum::Init(double x_length, double y_length, double dz_dtheta, int n_grids_x, int n_grids_y,
		    CoordSys coord_sys, BdyCond bdy_cond_up, BdyCond bdy_cond_left, BdyCond bdy_cond_right, BdyCond bdy_cond_down,
		    double init_mass_solid, double k_disv_per_area, double k_rect, double Csat)
{
  Sebum::Init(x_length, y_length, dz_dtheta, n_grids_x, n_grids_y, 
	      coord_sys, bdy_cond_up, bdy_cond_left, bdy_cond_right, bdy_cond_down);

  m_mass_solid = init_mass_solid;
  if (m_mass_solid>0) m_b_has_solid=true;
  m_k_disv_per_area = k_disv_per_area;
  updateKdisv(Cube);

  m_k_rect = k_rect;
  if (m_k_rect>0) m_b_has_react=true;

  m_Csat = Csat;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++
functions for computing the right-hand size of the odes */

/* ++++++++++
   Special treatment is needed because of possible mass transfer
    between liquid sebum and solids
 */
void SurSebum::compODE_dydt (double t, const double y[], double f[])
{

  /* first call the function in base class to calculate diffusion */
  Comp::compODE_dydt(t, y, f); 

  /* then calculate the mass transfer between sebum and solid */
  int i, j, idx;

  if (m_b_has_solid){
    idx = 0; // assume the solid only disolves into the left-most grid of sebum
    f[idx] += m_k_disv * (m_Csat - y[idx]);
  }

  if (m_b_has_react){
    for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
      for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right
	
	idx = i*m_ny + j;
	f[idx] -= m_k_rect * y[idx];
 
      } // for j
    } // for i
  }
  // TODO: handle solid mass balance somewhere
}
/* ----------------- */

void SurSebum::updateKdisv(CryShape shape)
{
  double volume, area, len;

  volume = m_mass_solid / m_rho_solid;

  switch (shape) {
  case Sphere :
    SayBye("Crystal shape not implemented yet");
    break;
  case Cube :
    len = pow(volume, 1.0/3);
    area = len*len*3; // only half of the surface, assuming axis-symmetric at the centre of the crystal
    m_k_disv = m_k_disv_per_area * area;
    break;
  default :
    SayBye("Crystal shape not implemented");
    break;
  }

}
