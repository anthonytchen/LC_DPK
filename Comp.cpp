#include "stdafx.h"
#include "Comp.h"

/* Structure and definition for parallel computing */
#define NTHREADS 1 // number of threads for parallel computing
struct pthread_struct {
	Comp *comp_obj;
	double t;
	const double *y;
	double *f;
	int x_start, x_end, y_start, y_end;
};
/* ------------------- */



/* This function will be called at the beginning of the compartment-specific Init function
   thus it only sets up member variables generic to all types of compartment
 */
void Comp::Init( CoordSys coord_sys, double dz_dtheta,
		 BdyCond bdy_cond_up, BdyCond bdy_cond_left, BdyCond bdy_cond_right, BdyCond bdy_cond_down)
{	
  m_coord_sys = coord_sys;
  m_dz_dtheta = dz_dtheta;

  m_BdyCond_up = bdy_cond_up;
  m_BdyCond_left = bdy_cond_left;
  m_BdyCond_right = bdy_cond_right;
  m_BdyCond_down = bdy_cond_down;

  m_gridsBdyUp = NULL;
  m_gridsBdyLeft = NULL;
  m_gridsBdyRight = NULL;
  m_gridsBdyDown = NULL;
  m_MassIn_up = m_MassIn_left = m_MassOut_right = m_MassOut_down = NULL;
  m_grids = NULL;
  m_conc1D = m_coord1D = NULL;
  m_ode_Jacobian = NULL;
}

/*
void StraCorn::updateBoundary(Grid* up, Grid* down, Grid* left, Grid* right)
{
  if (up!=NULL) m_gridBdyUp.set(up);
  if (down!=NULL) m_gridBdyDown.set(down);
  if (left!=NULL) m_gridBdyLeft.set(left);
  if (right!=NULL) m_gridBdyRight.set(right);
}
*/


/* Calculate the interfacial area between gridThiis and a neighbouring grid
   direction: [0] = up; [1] = left; [2] = right; [3] = down
*/
double Comp::compInterArea(Grid gridThiis, int direction)
{
  double area, r1, r2, pi_alpha_360;

  switch (m_coord_sys) {

  case Cartesian :

    if (direction==0 || direction==3)
      area = gridThiis.m_dy * gridThiis.m_dz;
    else if (direction==1 || direction==2)
      area = gridThiis.m_dx * gridThiis.m_dz;
    else
      SayBye("Option not valid");
    break;

  case Cylindrical :

    r1 = gridThiis.m_y_coord;
    r2 = gridThiis.m_y_coord + gridThiis.m_dy;
    pi_alpha_360 = M_PI * gridThiis.m_dz / 360;

    if (direction==0 || direction==3)
      area = pi_alpha_360 * (r2*r2 - r1*r1);
    else if (direction==1)
      area = gridThiis.m_dx * pi_alpha_360 * 2 * r1;
    else if (direction==2)
      area = gridThiis.m_dx * pi_alpha_360 * 2 * r2;
    else
      SayBye("Option not valid");
    break;

  default :
    SayBye("Coordinate system not implemented");
    break;
  }

}

/* compute the volume of this grid */
double Comp::compVolume(Grid gridThiis)
{
  double volume, r1, r2;

  switch (m_coord_sys) {

  case Cartesian :
    volume = gridThiis.m_dx * gridThiis.m_dy * gridThiis.m_dz;
    break;
  case  Cylindrical :
    r1 = gridThiis.m_y_coord;
    r2 = gridThiis.m_y_coord + gridThiis.m_dy;
    volume = gridThiis.m_dx * (M_PI*gridThiis.m_dz/360) * (r2*r2 - r1*r1);
    break;
  default :
    SayBye("Coordinate system not implemented");
    break;

  }
}

/* compute mass transfer between gridThhis and neighbouring grids to the right
   whose meshing does not match exactly to gridTThis, e.g. 
   gridTThis may interface with multiple grids to the right */
double Comp::compMassIrregGridsRight(Grid gridThhis, double conc_this)
{
  int i;
  double x1_this, x2_this, currentX, nextX, x_length, z_length, thd, deriv_this, deriv_other, area;
  double massIntoThis, mass;
  Grid *gridOther = NULL;
  
  thd = gridThhis.m_dx * 1e-3; // to compare whether two real numbers are the same

  x1_this = gridTThis.m_x_coord;
  x2_this = x1_this + gridTThis.m_dx;

  massTrans = currentX = 0;

  z_length = compInterArea(gridThiis, 2) / gridThhis.m_dx; // interfacial z_length

  for (i=0; i<m_n_gridsBdyRight; i++) {
    currentX = m_gridsBdyRight[i].m_x_coord;
    nextX = currentX + m_gridsBdyRight[i].m_dx;

    if ( x1_this > currentX-thd ) {

      if ( x2_this < nextX+thd ) { // gridThhis is contained between currentX & nextX

	x_length = x2_this - x1_this;

	gridOther = &m_gridsBdyRight[i];
	conc_other = gridOther->m_concChem;
     
	flux = gridThiis.compFlux( gridOther, conc_this, conc_other, 
				    gridThiis.m_dx/2, gridOther->m_dx/2, &deriv_this, &deriv_other);
	area = z_length * x_length;
	massTrans += area * flux;

      }
    }

  }

  return massTrans;
}


/* compute mass transfer between gridThhis and neighbouring grids downward
   whose meshing does not match exactly to gridTThis, e.g. 
   gridTThis may interface with multiple grids downward */
double Comp::compMassIrregGridsDown(Grid gridThhis)
{
  double massTransfer;
  return massTransfer;
}


/* functions for computing the right-hand size of the odes */

void Comp::compODE_dydt(double t, const double y[], double f[])
{
  int i, rc;
	
  if (NTHREADS==1) {
    compODE_dydt_block (t, y, f, 0, m_nx, 0, m_ny);
  } else {		
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
}

// the static container function needed for using multipe threads
void* Comp::static_compODE_dydt_block_threads(void *paras)
{
  struct pthread_struct p = *((struct pthread_struct *) paras);
  p.comp_obj->compODE_dydt_block (p.t, p.y, p.f, p.x_start, p.x_end, p.y_start, p.y_end);
}

// the actual funtion to calculate dy/dy
void Comp::compODE_dydt_block (double t, const double y[], double f[], 
			       int idx_x_start, int idx_x_end, int idx_y_start, int idx_y_end)
{
  int i, j, idx_this, idx_other, dim;
  double flux, mass, mass_transfer_rate, conc_this, conc_other, deriv_this, deriv_other,
    deriv_this_sum, volume_this, area;
	
  dim = m_nx*m_ny;
	
  assert(idx_x_start>=0 && idx_x_start<=m_nx);
  assert(idx_x_end>idx_x_start && idx_x_end<=m_nx);
  assert(idx_y_start>=0 && idx_y_start<=m_ny);
  assert(idx_y_end>idx_y_start && idx_y_end<=m_ny);
	
  Grid *gridThiis, *gridUp, *gridLeft, *gridRight, *gridDown;
	
  gridThiis = gridUp = gridLeft = gridRight = gridDown = NULL;

  /* todo: initialise / reset the boundary grids !!! */
  //  if ( idx_x_start == 0 ) m_mass_in = 0; // re-set mass transferred into SC, ready for calculation
  //if ( idx_x_end == m_nx ) m_mass_out = 0; // re-set mass transferred out of SC, ready for calculation

  // Calculate diffused mass
  for ( i=idx_x_start; i<idx_x_end; i++ ) { // x direction up to down
    for ( j=idx_y_start; j<idx_y_end; j++ ) { // y direction left to right
				
      mass_transfer_rate = 0;
      idx_this = i*m_ny+j;
			
      gridThiis = &m_grids[idx_this];
      conc_this = y[idx_this];
      volume_this = compVolume(*gridThiis);
			
      // Setup the neighbouring grids and calculate the mass transfer rate.
      // in the following order: up, left, (self), right, down.
			
      /* diffusion from up */

      area = compInterArea(*gridThiis, 0); // interfacial area

      if ( i==0 ) { // topmost layer, its top is up boundary
	switch (m_BdyCond_up) {
	case ZeroFlux :
	  break;
	case FromOther :
	  mass_transfer_rate += m_MassIn_up[j];
	  break;
	}
      }
      else { // not topmost layer
	idx_other = (i-1)*m_ny+j;
	gridUp = &m_grids[idx_other];
	conc_other = y[idx_other];

	flux = gridThiis->compFlux( gridUp, conc_this, conc_other,
				    gridThiis->m_dx/2, gridUp->m_dx/2, &deriv_this, &deriv_other);
	mass_transfer_rate += area * flux;
      }


      /* diffusion from left */

      area = compInterArea(*gridThiis, 1); // interfacial area

      if ( j==0 ) { // leftmost layer, its left is left boundary
	switch (m_BdyCond_left) {
	case ZeroFlux : 
	  break;
	case ZeroConc :
	  SayBye("to be implemented");
	  break;
	case Periodic :
	  idx_other = i*m_ny+m_ny-1;
	  gridLeft = &m_grids[idx_other];
	  conc_other = y[idx_other];
	  flux = gridThiis->compFlux( gridLeft, conc_this, conc_other, 
				      gridThiis->m_dy/2, gridLeft->m_dy/2, &deriv_this, &deriv_other);	  
	  mass_transfer_rate += area * flux;	
	  break;
	case FromOther :
	  mass_transfer_rate += m_MassIn_left[i];
	  break;
	}
      } 
      else { // not leftmost layer
	idx_other = i*m_ny+j-1;
	gridLeft = &m_grids[idx_other];
	conc_other = y[idx_other];
	flux = gridThiis->compFlux( gridLeft, conc_this, conc_other, 
				    gridThiis->m_dy/2, gridLeft->m_dy/2, &deriv_this, &deriv_other);
	
	mass_transfer_rate += area * flux;
      }

      /* diffusion from right */

      area = compInterArea(*gridThiis, 2); // interfacial area

      if ( j==m_ny-1 ) { // rightmost layer

	switch (m_BdyCond_right) {
	case ZeroFlux : 
	  break;
	case ZeroConc :
	  SayBye("to be implemented");
	  break;
	case Periodic :
	  idx_other = i*m_ny;
	  gridRight = &m_grids[idx_other];
	  conc_other = y[idx_other];
	  flux = gridThiis->compFlux( gridRight, conc_this, conc_other, 
				      gridThiis->m_dy/2, gridRight->m_dy/2, &deriv_this, &deriv_other);
	  mass_transfer_rate += area * flux;
	  break;
	case FromOther :
	  mass = compMassIrregGridsRight(*gridThiis);
	  mass_transfer_rate += mass;
	  break;
	}

      } 
      else { // not rightmost layer
	idx_other = i*m_ny+j+1;
	gridRight = &m_grids[idx_other];
	conc_other = y[idx_other];
	flux = gridThiis->compFlux( gridRight, conc_this, conc_other, 
				    gridThiis->m_dy/2, gridRight->m_dy/2, &deriv_this, &deriv_other);
	mass_transfer_rate += area * flux;
      }	


      /* diffusion from down */

      area = compInterArea(*gridThiis, 2); // interfacial area	

      if ( i==m_nx-1 ) { // bottom layer

	switch (m_BdyCond_down) {
	case ZeroFlux : 
	  break;
	case ZeroConc :
	  SayBye("to be implemented");
	  break;
	case FromOther :
	  mass = compMassIrregGridsDown(*gridThiis);
	  mass_transfer_rate += mass;
	  break;
	}

      } 
      else { // not downest layer
	idx_other = (i+1)*m_ny+j;
	gridDown = &m_grids[idx_other];
	conc_other = y[idx_other];      
	flux = gridThiis->compFlux( gridDown, conc_this, conc_other, 
				    gridThiis->m_dx/2, gridDown->m_dx/2, &deriv_this, &deriv_other);
	mass_transfer_rate += area * flux;
      }
		
      f[idx_this] = mass_transfer_rate / volume_this;

    } // for j
  } // for i
	
}
/* ----------------- */

	
/*  ++++++++++++++++++++++++++++++++++
	I/O functions
	++++++++++++++++++++++++++++++++++ */

void Comp::displayGrids()
{
  assert( m_grids );

  int i, j, idx, gsl_errno;;
  printf("# of grids: [x] %d, [y] %d in stratum corneum\n", m_nx, m_ny);

  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right	

      idx = i*m_ny + j;			
      if ( !strcmp(m_grids[idx].m_name, "LP") )
	printf("L ");
      else if ( !strcmp(m_grids[idx].m_name, "CC") )
	printf("C ");
      else
	gsl_error ("subtype name unknown", __FILE__, __LINE__, gsl_errno); 
				
    } // for j
    printf("\n");
  } // for i
  fflush(stdout);
}

void Comp::getGridsConc(double *fGridsConc, int dim)
{
  // Return concentration at the grids in fGridsConc
  assert( m_grids && fGridsConc && dim==m_nx*m_ny);

  int i, j, idx;
	
  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right
       idx = i*m_ny + j;
       fGridsConc[idx] = m_grids[idx].getConcChem();
     }
  }
}


void Comp::saveGrids(bool b_1st_time, const char fn[])
{
  assert( m_grids );

  FILE *file = NULL;
  int i, j, idx;

  if ( b_1st_time )
    file = fopen(fn, "w");
  else 
    file = fopen(fn, "a");
	
  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right		
      idx = i*m_ny + j;
      fprintf(file, "%.5e\t", m_grids[idx].getConcChem());
    }
    fprintf(file, "\n");
  }

  fclose(file);
}

void Comp::getXCoord(double *coord_x, int dim)
{
  assert( m_grids && coord_x && dim==m_nx*m_ny );

  int i, j, idx;

  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right		
      idx = i*m_ny + j;
      coord_x[idx] = m_grids[idx].m_x_coord;
    }
  }

}

void Comp::getYCoord(double *coord_y, int dim)
{
  assert( m_grids && coord_y && dim==m_nx*m_ny );

  int i, j, idx;

  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right		
      idx = i*m_ny + j;
      coord_y[idx] = m_grids[idx].m_y_coord;
    }
  }
}

void Comp::saveCoord(const char fn_x[], const char fn_y[])
{
  assert( m_grids );

  FILE *file_x, *file_y;
  char fn1[1024], fn2[1024];
  int i, j, idx;

  // save grids
  strcpy(fn1, fn_x); strcat(fn1, ".sc");
  strcpy(fn2, fn_y); strcat(fn2, ".sc");

  file_x = fopen(fn1, "w");
  file_y = fopen(fn2, "w");

  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right		
      idx = i*m_ny + j;
      fprintf(file_x, "%.5e\t", m_grids[idx].m_x_coord);
      fprintf(file_y, "%.5e\t", m_grids[idx].m_y_coord);			
    }
    fprintf(file_x, "\n");
    fprintf(file_y, "\n");
  }

  fclose(file_x);
  fclose(file_y);
}
/*  ---- END <I/O functions> ---- */
