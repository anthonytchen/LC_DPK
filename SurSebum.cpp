#include "stdafx.h"
#include "SurSebum.h"

/*!
  Initialise surface sebum (SurSebum).
  Args:
      init_mass_solid: initial mass of solid particle; if negative it doesn't exist
 */
void SurSebum::Init(double x_length, double y_length, double dz_dtheta, int n_grids_x, int n_grids_y,
		    CoordSys coord_sys, BdyCond bdy_cond_up, BdyCond bdy_cond_left, BdyCond bdy_cond_right, BdyCond bdy_cond_down,
		    Crystal crystal, double init_mass_solid, double k_disv_per_area, double k_rect, double Csat)
{
  Sebum::Init(x_length, y_length, dz_dtheta, n_grids_x, n_grids_y, 
	      coord_sys, bdy_cond_up, bdy_cond_left, bdy_cond_right, bdy_cond_down);

  m_mass_solid = init_mass_solid;
  if (m_mass_solid>0) m_b_has_solid=true;

  m_k_rect = k_rect;
  if (m_k_rect>0) m_b_has_react=true;

  m_Csat = Csat;

  if (m_b_has_solid) {
    // add one dimension for the solid dissolution
    //   currently only implemented to dissolve into the first grid in the surface sebum
    m_crystal.shape = crystal.shape;
    m_crystal.density = crystal.density;
    m_crystal.dim = crystal.dim;
    m_crystal.area = crystal.area;
    for (int i=0; i<crystal.dim; i++)
      m_crystal.len[i] = crystal.len[i];
    
    m_k_disv_per_area = k_disv_per_area;
    updateKdisv(m_crystal.shape);

    m_dim += 1; 
  }
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
  double V, Cdiff;

  if (m_b_has_solid){
    idx = m_dim-1; // this is the solid index
    updateKdisv(m_crystal.shape, y[idx]);

    idx = 0; // assume the solid only disolves into the left-most grid of sebum
    // idx = m_dim-2; // assume the solid only disolves into the right-most grid of sebum
    Cdiff =  m_Csat - y[idx];
    V = compVolume(m_grids[idx]);
    f[idx] += m_k_disv * Cdiff / V;
    
    idx = m_dim-1; // this is the solid index
    f[idx] = -m_k_disv * Cdiff; // mass balance for the solid
  }

  if (m_b_has_react){
    for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
      for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right
	
	idx = i*m_ny + j;
	f[idx] -= m_k_rect * y[idx];
 
      } // for j
    } // for i
  }

}
/* ----------------- */

void SurSebum::updateKdisv(CryShape shape, double mass_solid)
{
  double volume, area, len1, len, tmp, factor;

  if ( mass_solid < 0 ) {
    m_k_disv = .0;
    return;
  }
  
  volume = mass_solid / m_crystal.density;

  switch (shape) {
  case Sphere :
    SayBye("Crystal shape not implemented yet");
    break;
    
  case Cube :
    SayBye ("Crystal shape out of date; please use HyperRect");
    len = sqrt(volume/m_dz_dtheta);
    area = len*len*2 + len*m_dz_dtheta*4;
    break;
    
  case HyperRect : // hyper-rectangular

    if (m_coord_sys == Cartesian){
      factor = sqrt( volume/m_dz_dtheta / m_crystal.area);
      len = m_crystal.len[0]*factor;
      len1 = m_crystal.len[1]*factor;
      area = len*len1*2 + (len+len1)*m_dz_dtheta*2;
    }
    else if (m_coord_sys == Cylindrical ) {
      SayBye("Cylindrical coordinate for hyper-reactangular solids is not implmented");
    }
    else
      SayBye("Coordinate system not implemented");

    break;

  case BottomOnly : // dissolution only depends on the bottom area which is fixed

    if (m_coord_sys == Cartesian)
      area = m_crystal.len[0] * m_dz_dtheta;
    else if (m_coord_sys == Cylindrical )
      area = M_PI * (m_crystal.len[0]*m_crystal.len[0]) * m_dz_dtheta / 360; // top view area, same as bottom area
    else
      SayBye("Coordinate system not implemented");
    
    break;

  default :
    SayBye("Crystal shape not implemented");
    break;
  }
  // printf("area = %.5e\n", area);

  m_k_disv = m_k_disv_per_area * area;

}


void SurSebum::updateKdisv(CryShape shape)
{
  updateKdisv(shape, m_mass_solid);
}
/*
void SurSebum::setGridsConc(const double fGridsConc[], int dim)
{

  if (m_b_has_solid) {
    Comp::setGridsConc(fGridsConc, dim-1);
    m_mass_solid = fGridsConc[dim-1];
  }
  else 
    Comp::setGridsConc(fGridsConc, dim);
}

void SurSebum::getGridsConc(double *fGridsConc, int dim)
{

  if (m_b_has_solid) {
    Comp::getGridsConc(fGridsConc, dim-1);
    fGridsConc[dim-1] = m_mass_solid;
  }
  else 
    Comp::setGridsConc(fGridsConc, dim);
}
*/

void SurSebum::saveGrids(bool b_1st_time, const char fn[])
{
  FILE *file = NULL;

  if (m_b_has_solid) {
    if ( b_1st_time ) {
      file = fopen(fn, "w");
      b_1st_time = false;
    }
    else 
      file = fopen(fn, "a");

    fprintf(file, "%.5e\n", m_mass_solid);
    fclose(file);
#ifdef _DEBUG_
      printf( "\tSolid\t\t %.3e\n", m_mass_solid );
#endif
  }

  Comp::saveGrids(b_1st_time, fn);

}
