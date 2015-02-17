#include "stdafx.h"
#include "ViaEpd.h"

/* Structure and definition for parallel computing */
#define NTHREADS 8 // number of threads for parallel computing
struct pthread_struct {
	ViaEpd *ve_obj;
	double t;
	const double *y;
	double *f;
	int x_start, x_end, y_start, y_end;
};
/* ------------------- */



/*
 */
void ViaEpd::Init(int x_length, int y_length, double dz, int n_layer_x)
{	

  /* set up some constant values */

  // density of lipid, keratin and water
  m_rou_lipid = 1e3; // kg m^{-3}
  m_rou_keratin = 1.2e3; // kg m^{-3}
  m_rou_water = 1e3; // kg m^{-3}

  m_boundary_cond = 1; // boundary condition for left/right; 
                       //  0: zero flux; 1: periodic

  /* ---- */

  m_grids = NULL;
  m_ode_Jacobian = NULL;

  m_x_length = x_length;
  m_y_length = y_length;
  m_dz = dz;
	
  m_nx = n_layer_x;
  m_ny = 1; // lateral direction only 1 grid
}

void ViaEpd::Release()
{
  int i, j, idx;
  if (!m_grids) {
    for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
      for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right
	idx = i*m_ny + j;
	m_grids[idx].Release();
      }
    }
    delete [] m_grids;
  }
  m_gridBdyUp.Release();
  m_gridBdyDown.Release();
  m_gridBdyLeft.Release();
  m_gridBdyRight.Release();
}

/*
  Chemical parameters:
    MW: molecular weight
    Kow: partition coefficient between octanol and water
    pKa: the ionisation of the chemical
 */
void ViaEpd::createGrids(double MW, double Kow, double pKa, double coord_x_now)
{
  int i, j, idx, idx_x, idx_y, gsl_errno;
  double dx, dy, coord_x, coord_y;

  // initialise boundary grids
  m_gridBdyUp.Init("SC",    0, Kow, 0, 0, m_dz); // infinite source, can be updated to diminishing vehicle by using updateBoundary()
  m_gridBdyDown.Init("SK",  0, Kow, 0, 0, m_dz); // infinite sink, can be updated to underlying viable epidermis by using updateBoundary()
  m_gridBdyLeft.Init("SK",  0, Kow, 0, 0, m_dz); // infinite sink
  m_gridBdyRight.Init("SK", 0, Kow, 0, 0, m_dz); // infinite sink

  dx = m_x_length / m_nx;
  dy = m_y_length / m_ny;
	
  m_grids = new Grid[m_nx*m_ny]; // organised in row dominant

  coord_x = coord_x_now;   coord_y = 0;
  struct Point current_point;
  setPoint(current_point, coord_x, coord_y, dx, dy, "VE", "VE");
    
  idx_x = idx_y = 0;

  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right
			
      idx = i*m_ny + j;

      m_grids[idx].InitVE(MW, Kow, pKa, current_point.x_coord, current_point.y_coord, current_point.dx, current_point.dy, m_dz);

      // update current_point
      if (j==m_ny-1) { // last element in the lateral direction, move down
	coord_x += dx; 	coord_y = 0;
      } else { // not the last element in the lateral direction, thus move to the right
	coord_y += dy;
      }

      setPoint(current_point, coord_x, coord_y, dx, dy, "VE", "VE");
      
    } // for j
  } // for i

}

void ViaEpd::updateBoundary(Grid* up, Grid* down, Grid* left, Grid* right)
{
  if (up!=NULL) m_gridBdyUp.set(up);
  if (down!=NULL) m_gridBdyDown.set(down);
  if (left!=NULL) m_gridBdyLeft.set(left);
  if (right!=NULL) m_gridBdyRight.set(right);
}

/* functions for computing the right-hand size of the odes */

// todo: remove repeated calculation of mass transfer
void ViaEpd::compODE_dydt (double t, const double y[], double f[])
{
  int i, rc;
	
  if (NTHREADS==1) {
    compODE_dydt_block (t, y, f, 0, m_nx, 0, m_ny);
  } else {		
    struct pthread_struct p[NTHREADS];
    pthread_t threads[NTHREADS];
			
    for ( i=0; i < NTHREADS; i++ ) {
      p[i].ve_obj = this;
      p[i].t=t; p[i].y=y; p[i].f=f;			
      p[i].y_start=0; p[i].y_end=m_ny;
			
      p[i].x_start=i*m_nx/NTHREADS; p[i].x_end=(i+1)*m_nx/NTHREADS;
      pthread_create(&threads[i], NULL, static_compODE_dydt_block_threads, (void *) &p[i]);
    }
		
    for ( i=0; i<NTHREADS; i++ ) {
      rc = pthread_join(threads[i], NULL);
      assert(rc==0);
    }		
  }
}

// the static container function needed for using multipe threads
void* ViaEpd::static_compODE_dydt_block_threads(void *paras)
{
  struct pthread_struct p = *((struct pthread_struct *) paras);
  p.ve_obj->compODE_dydt_block (p.t, p.y, p.f, p.x_start, p.x_end, p.y_start, p.y_end);
}

// the actual funtion to calculate dy/dy
void ViaEpd::compODE_dydt_block (double t, const double y[], double f[], 
				 int idx_x_start, int idx_x_end, int idx_y_start, int idx_y_end)
{
  int i, j, idx_this, idx_other, dim;
  double flux, mass_transfer_rate, conc_this, conc_other, deriv_this, deriv_other,
    deriv_this_sum, volume_this;
	
  dim = m_nx*m_ny;
	
  assert(idx_x_start>=0 && idx_x_start<=m_nx);
  assert(idx_x_end>idx_x_start && idx_x_end<=m_nx);
  assert(idx_y_start>=0 && idx_y_start<=m_ny);
  assert(idx_y_end>idx_y_start && idx_y_end<=m_ny);
	
  /* Assumptions:
     -- The vehicle is a infinite source.
     -- Left/right/down boundaries are sink (concentration always zero)
  */
  Grid *gridThiis, *gridUp, *gridLeft, *gridRight, *gridDown;
	
  gridThiis = gridUp = gridLeft = gridRight = gridDown = NULL;
	
  // Calculate diffused mass
  for ( i=idx_x_start; i<idx_x_end; i++ ) { // x direction up to down
    for ( j=idx_y_start; j<idx_y_end; j++ ) { // y direction left to right
				
      mass_transfer_rate = 0;
      idx_this = i*m_ny+j;
			
      gridThiis = &m_grids[idx_this];
      conc_this = y[idx_this];
      volume_this = gridThiis->m_dx * gridThiis->m_dy * gridThiis->m_dz;
			
      // Setup the neighbouring grids and calculate the mass transfer rate.
      // If Jacobian required, calculate the Jacobian in the following order:
      //	up, left, self, right, down.
			
      // diffusion from up

      if ( i==0 ) { // topmost layer, its top is up boundary
	gridUp = &m_gridBdyUp;
	conc_other = m_gridBdyUp.m_concChem;
      } else {
	idx_other = (i-1)*m_ny+j;
	gridUp = &m_grids[idx_other];
	conc_other = y[idx_other];
      }
      flux = gridThiis->compFlux( gridUp, conc_this, conc_other,
				  gridThiis->m_dx/2, gridUp->m_dx/2, &deriv_this, &deriv_other);
      mass_transfer_rate += gridThiis->m_dy*gridThiis->m_dz * flux;			
      if (m_ode_Jacobian!=NULL) {
	deriv_this_sum = deriv_this / gridThiis->m_dx;
	if ( i!=0 ) 
	  m_ode_Jacobian[ idx_this*dim + idx_other ] = deriv_other / gridThiis->m_dx;
      }

      // diffusion from left

      if ( j==0 ) { // leftmost layer, its left is sink
	
	if ( m_boundary_cond == 0 ) {
	  gridLeft = &m_gridBdyLeft;
	  conc_other = m_gridBdyLeft.m_concChem;
	  flux = 0; // left impermeable
	} else if ( m_boundary_cond == 1 ) {
	  idx_other = i*m_ny+m_ny-1;
	  gridLeft = &m_grids[idx_other];
	  conc_other = y[idx_other];
	  flux = gridThiis->compFlux( gridLeft, conc_this, conc_other, 
				      gridThiis->m_dy/2, gridLeft->m_dy/2, &deriv_this, &deriv_other);	  
	}
	
      } else {
	idx_other = i*m_ny+j-1;
	gridLeft = &m_grids[idx_other];
	conc_other = y[idx_other];
	flux = gridThiis->compFlux( gridLeft, conc_this, conc_other, 
				    gridThiis->m_dy/2, gridLeft->m_dy/2, &deriv_this, &deriv_other);
      }	
      
	
      mass_transfer_rate += gridThiis->m_dx*gridThiis->m_dz * flux;
      if (m_ode_Jacobian!=NULL) {
	deriv_this_sum += deriv_this / gridThiis->m_dy;
	if ( j!=0 ) 
	  m_ode_Jacobian[ idx_this*dim + idx_other ] = deriv_other / gridThiis->m_dy;
      }

      // diffusion from right

      if ( j==m_ny-1 ) { // right layer, its right is sink

	if ( m_boundary_cond == 0 ) {
	  gridRight = &m_gridBdyRight;
	  conc_other = m_gridBdyRight.m_concChem;
	  flux = 0; // left impermeable
	} else if ( m_boundary_cond == 1 ) {
	  idx_other = i*m_ny;
	  gridRight = &m_grids[idx_other];
	  conc_other = y[idx_other];
	  flux = gridThiis->compFlux( gridRight, conc_this, conc_other, 
				      gridThiis->m_dy/2, gridRight->m_dy/2, &deriv_this, &deriv_other);
	}
	
      } else {
	idx_other = i*m_ny+j+1;
	gridRight = &m_grids[idx_other];
	conc_other = y[idx_other];
	flux = gridThiis->compFlux( gridRight, conc_this, conc_other, 
				    gridThiis->m_dy/2, gridRight->m_dy/2, &deriv_this, &deriv_other);
      }	

      mass_transfer_rate += gridThiis->m_dx*gridThiis->m_dz * flux;
      if (m_ode_Jacobian!=NULL) {
	deriv_this_sum += deriv_this / gridThiis->m_dy;
	if ( j!=m_ny-1 ) 
	  m_ode_Jacobian[ idx_this*dim + idx_other ] = deriv_other / gridThiis->m_dy;
      }


      // diffusion from down
			
      if ( i==m_nx-1 ) { // bottom layer, its down is down boundary
	gridDown = &m_gridBdyDown;
	conc_other = m_gridBdyDown.m_concChem;
      } else {
	idx_other = (i+1)*m_ny+j;
	gridDown = &m_grids[idx_other];
	conc_other = y[idx_other];
      }
      flux = gridThiis->compFlux( gridDown, conc_this, conc_other, 
				  gridThiis->m_dx/2, gridDown->m_dx/2, &deriv_this, &deriv_other);

      mass_transfer_rate += gridThiis->m_dy*gridThiis->m_dz * flux;
      if (m_ode_Jacobian!=NULL) {
	deriv_this_sum += deriv_this / gridThiis->m_dx;
	if ( i!=m_nx-1 ) 
	  m_ode_Jacobian[ idx_this*dim + idx_other ] = deriv_other / gridThiis->m_dx;
      }

			
      f[idx_this] = mass_transfer_rate / volume_this;
      if (m_ode_Jacobian!=NULL) 
	m_ode_Jacobian[ idx_this*dim + idx_this ] = deriv_this_sum;

    } // for j
    //printf("\n");
  } // for i
	
}
/* ----------------- */

	
/*  +++  I/O functions +++++++++ */

void ViaEpd::displayGrids()
{
  assert( m_grids );

  int i, j, idx;
  printf("# of grids: [x] %d, [y] %d in viable epidermis\n", m_nx, m_ny);

  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right	
      printf("V ");
    }
    printf("\n");
  }
  fflush(stdout);
}

void ViaEpd::saveGrids(bool b_1st_time, const char fn[])
{
  assert( m_grids );

  FILE *file = NULL;
  int i, j, idx;

  // save grids
  if ( b_1st_time )
    file = fopen(fn, "w");
  else 
    file = fopen(fn, "a");
	
  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right		
      idx = i*m_ny + j;
      fprintf(file, "%.5e\t", m_grids[idx].getConcChem());
    } // for j
    fprintf(file, "\n");
  } // for i

  fclose(file);
}

void ViaEpd::saveCoord(const char fn_x[], const char fn_y[])
{
  assert( m_grids );

  FILE *file_x, *file_y;
  int i, j, idx;

  // save grids
  file_x = fopen(fn_x, "w");
  file_y = fopen(fn_y, "w");

  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right		
      idx = i*m_ny + j;
      fprintf(file_x, "%.5e\t", m_grids[idx].m_x_coord);
      fprintf(file_y, "%.5e\t", m_grids[idx].m_y_coord);
    } // for j
    fprintf(file_x, "\n");
    fprintf(file_y, "\n");
  } // for i

  fclose(file_x);
  fclose(file_y);
}
/*  ------------ END <I/O functions> -------------------- */
