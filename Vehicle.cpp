#include "stdafx.h"
#include "Vehicle.h"

/*
 */
void Vehicle::Init(double x_length, double y_length, double dz_dtheta, int n_grids_x, int n_grids_y,
		   double init_concChem,
		  CoordSys coord_sys, BdyCond bdy_cond_up, BdyCond bdy_cond_left, BdyCond bdy_cond_right, BdyCond bdy_cond_down)
{	
  // call Init of the base class Comp
  Comp::Init (coord_sys, dz_dtheta, bdy_cond_up, bdy_cond_left, bdy_cond_right, bdy_cond_down);

  m_x_length = x_length;
  m_y_length = y_length;
  m_dz_dtheta = dz_dtheta;
	
  m_nx = n_grids_x;
  m_ny = n_grids_y;
  m_init_concChem = init_concChem;

  m_T = 309; // temperature (Kelvin)
  m_eta = 7.1E-4; // water viscosity at above temperature (Pa s),
  m_K_vw = 1; // partition coefficient between vehicle and water (1 means vehicle is water)
}


/*
 */
void Vehicle::createGrids(Chemical chem, double coord_x_now)
{
  int i, j, idx, idx_x, idx_y;
  double dx, dy, coord_x, coord_y;

  dx = m_x_length / m_nx;
  dy = m_y_length / m_ny;
	
  m_grids = new Grid[m_nx*m_ny]; // organised in row dominant

  coord_x = coord_x_now;   coord_y = 0;
  struct Point current_point;
  setPoint(current_point, coord_x, coord_y, dx, dy, "VH", "VH");
    
  idx_x = idx_y = 0;

  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right
			
      idx = i*m_ny + j;

      m_grids[idx].InitVH("VH", chem, m_init_concChem, current_point.x_coord, current_point.y_coord, current_point.dx, current_point.dy, m_dz_dtheta, m_T, m_eta, m_K_vw);

      // update current_point
      if (j==m_ny-1) { // last element in the lateral direction, move down
	coord_x += dx; 	coord_y = 0;
      } else { // not the last element in the lateral direction, thus move to the right
	coord_y += dy;
      }

      setPoint(current_point, coord_x, coord_y, dx, dy, "VH", "VH");
      
    } // for j
  } // for i
}

	
/*  +++  I/O functions +++++++++ */

void Vehicle::saveCoord(const char fn_x[], const char fn_y[])
{
  Comp::saveCoord(fn_x, fn_y, ".vh");
}
/*  ------------ END <I/O functions> -------------------- */
