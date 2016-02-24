#include "stdafx.h"
#include "Dermis.h"

/*
 */
void Dermis::Init(double x_length, double y_length, double dz_dtheta, int n_grids_x, int n_grids_y, bool b_has_blood,
		  CoordSys coord_sys, BdyCond bdy_cond_up, BdyCond bdy_cond_left, BdyCond bdy_cond_right, BdyCond bdy_cond_down)
{	
  // call Init of the base class Comp
  Comp::Init (coord_sys, dz_dtheta, bdy_cond_up, bdy_cond_left, bdy_cond_right, bdy_cond_down);

  m_bToBlood = b_has_blood; // whether dermis is connected to blood

  m_x_length = x_length;
  m_y_length = y_length;
  m_dz_dtheta = dz_dtheta;
	
  m_nx = n_grids_x;
  m_ny = n_grids_y;
  m_dim = m_nx * m_ny;
}

void Dermis::InitDermisBlood(double bld_skin_flow, double bld_fu, double par_de2blood, double bld_concChem, double skin_area)
{
  assert(m_bToBlood);

  m_bld_skin_flow = bld_skin_flow; // total blood flow rate in skin
  m_bld_fu = bld_fu; // fraction of unbounded solute in blood
  m_par_de2blood = par_de2blood; // partition coefficient from dermis to blood
  m_bld_concChem = bld_concChem; // initial concentration in blood
  m_dermis_totalV = m_x_length*skin_area; // total volume of dermis
}


/*
  Chemical parameters:
    MW: molecular weight
    Kow: partition coefficient between octanol and water
    pKa: the ionisation of the chemical
    frac_non_ion: fraction of non-ionisation
    frac_unbound: fraction of unbound to 2.7% albumin
 */
void Dermis::createGrids(Chemical chem, double coord_x_start, double coord_y_start)
{
  int i, j, idx, idx_x, idx_y, gsl_errno;
  double dx, dy, coord_x, coord_y;

  dx = m_x_length / m_nx;
  dy = m_y_length / m_ny;
	
  m_grids = new Grid[m_nx*m_ny]; // organised in row dominant

  coord_x = coord_x_start;   coord_y = coord_y_start;
  struct Point current_point;
  setPoint(current_point, coord_x, coord_y, dx, dy, "DE", "DE");
    
  idx_x = idx_y = 0;

  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right
			
      idx = i*m_ny + j;

      // For now, diffusion and partition coefficients in dermis are the same as those in viable epidermis
      m_grids[idx].InitVE_DE("DE", chem, 0, current_point.x_coord, current_point.y_coord, current_point.dx, current_point.dy, m_dz_dtheta);

      // update current_point
      if (j==m_ny-1) { // last element in the lateral direction, move down
	coord_x += dx; 	coord_y = coord_y_start;
      } else { // not the last element in the lateral direction, thus move to the right
	coord_y += dy;
      }

      setPoint(current_point, coord_x, coord_y, dx, dy, "DE", "DE");
      
    } // for j
  } // for i

}


void Dermis::updateBlood(double concChem)
{
  m_bld_concChem = concChem;
}

/* functions for computing the right-hand size of the odes */

void Dermis::compODE_dydt (double t, const double y[], double f[])
{

  /* first call the function in base class to calculate diffusion */
  Comp::compODE_dydt(t, y, f); 

  /* then calculate the convective terms into blood */
  int i, rc;

  if (m_bToBlood){ 
    
    m_mass_into_dermis = m_mass_outof_dermis = .0; // reset both for later calculation for blood

    if (NTHREADS==1 || NTHREADS > m_ny) {
      compODE_dydt_block (t, y, f, 0, m_nx, 0, m_ny);
    } 
    else {		
      struct pthread_struct p[NTHREADS];
      pthread_t threads[NTHREADS];
			
      for ( i=0; i < NTHREADS; i++ ) {
	p[i].comp_obj = this;
	p[i].t=t; p[i].y=y; p[i].f=f;
	p[i].x_start=0; p[i].x_end=m_nx;
			
	p[i].y_start=i*m_ny/NTHREADS; p[i].y_end=(i+1)*m_ny/NTHREADS;
	pthread_create(&threads[i], NULL, static_compODE_dydt_block_threads, (void *) &p[i]);
      }
		
      for ( i=0; i<NTHREADS; i++ ) {
	rc = pthread_join(threads[i], NULL);
	assert(rc==0);
      }		
    }
  } // only do above if convection into blood
}

// the static container function needed for using multipe threads
void* Dermis::static_compODE_dydt_block_threads(void *paras)
{
  struct pthread_struct p = *((struct pthread_struct *) paras);
  p.comp_obj->compODE_dydt_block (p.t, p.y, p.f, p.x_start, p.x_end, p.y_start, p.y_end);
}

// the actual funtion to calculate dy/dy
void Dermis::compODE_dydt_block (double t, const double y[], double f[], 
				 int idx_x_start, int idx_x_end, int idx_y_start, int idx_y_end)
{
  assert(m_bToBlood);

  int i, j, idx_this;
  double conc_this, volume_this, flow_this_grid, fin, fout; 

  assert(idx_x_start>=0 && idx_x_start<=m_nx);
  assert(idx_x_end>idx_x_start && idx_x_end<=m_nx);
  assert(idx_y_start>=0 && idx_y_start<=m_ny);
  assert(idx_y_end>idx_y_start && idx_y_end<=m_ny);
	
  Grid *gridThiis = NULL;
	
  // Calculate mass flow from blood to this grid

  for ( i=idx_x_start; i<idx_x_end; i++ ) { // x direction up to down
    for ( j=idx_y_start; j<idx_y_end; j++ ) { // y direction left to right
				
      idx_this = i*m_ny+j;
			
      gridThiis = &m_grids[idx_this];
      conc_this = y[idx_this];
      volume_this = compVolume(*gridThiis);

      flow_this_grid = m_bld_skin_flow * volume_this/m_dermis_totalV;
      fin = flow_this_grid * m_bld_concChem;
      fout = flow_this_grid * conc_this * (gridThiis->m_chemical.m_frac_unbound/m_bld_fu) / m_par_de2blood;

      f[idx_this] += (fin-fout)/volume_this;

      m_mass_into_dermis += fin;
      m_mass_outof_dermis += fout;
 
    } // for j
  } // for i
	
}
/* ----------------- */

	
/*  +++  I/O functions +++++++++ */

void Dermis::saveCoord(const char fn_x[], const char fn_y[])
{
  Comp::saveCoord(fn_x, fn_y, ".de");  
}
/*  ------------ END <I/O functions> -------------------- */
