#include "stdafx.h"
#include "Comp.h"

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
  m_n_gridsBdyRight = 0;
  m_n_gridsBdyDown = 0;

  // m_gridsBdyUp = NULL;
  // m_gridsBdyLeft = NULL;
  m_gridsBdyRight = NULL;
  m_gridsBdyDown = NULL;
  m_MassIn_up = m_MassIn_left = m_MassOut_right = m_MassOut_down = NULL;
  m_grids = NULL;
  m_conc1D = m_coord1D = NULL;

  m_gridSink.InitSK(); // setup a sink grid
}

void Comp::Release()
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

  if (m_BdyCond_up == FromOther)
    delete [] m_MassIn_up;

  if (m_BdyCond_left == FromOther)
    delete [] m_MassIn_left;

  if (m_BdyCond_right == FromOther){
    delete [] m_MassOut_right;
    delete [] m_gridsBdyRight;
  }
  if (m_BdyCond_down == FromOther){
    delete [] m_MassOut_down;
    delete [] m_gridsBdyDown;
  }
}

void Comp::createBoundary(int n_gridsBdyRight, int n_gridsBdyDown)
{
  if (m_BdyCond_up == FromOther)
    m_MassIn_up = new double[m_ny];

  if (m_BdyCond_left == FromOther)
    m_MassIn_left = new double[m_nx];

  if (m_BdyCond_right == FromOther) {
    assert (n_gridsBdyRight > 0);
    m_n_gridsBdyRight = n_gridsBdyRight;
    m_MassOut_right = new double [n_gridsBdyRight];
    m_gridsBdyRight = new Grid[n_gridsBdyRight];
  }

  if (m_BdyCond_down == FromOther) {
    assert (n_gridsBdyDown > 0);
    m_n_gridsBdyDown = n_gridsBdyDown;
    m_MassOut_down = new double [n_gridsBdyDown];
    m_gridsBdyDown = new Grid[n_gridsBdyDown];
  }
}

void Comp::setBoundaryGrids(Grid *gridsBdyRight, Grid *gridsBdyDown)
{
  int i;

  if (m_n_gridsBdyRight) {
    for (i=0; i< m_n_gridsBdyRight; i++)
      m_gridsBdyRight[i] = gridsBdyRight[i];
  }

  if (m_n_gridsBdyDown) {
    for (i=0; i< m_n_gridsBdyDown; i++)
      m_gridsBdyDown[i] = gridsBdyDown[i];
  }
}

void Comp::setBoundaryConc(double *concBdyRight, double *concBdyDown)
{
  int i;

  if (m_n_gridsBdyRight) {
    for (i=0; i< m_n_gridsBdyRight; i++)
      m_gridsBdyRight[i].m_concChem = concBdyRight[i];
  }

  if (m_n_gridsBdyDown) {
    for (i=0; i< m_n_gridsBdyDown; i++)
      m_gridsBdyDown[i].m_concChem = concBdyDown[i];
  }
}

void Comp::setBdyMassInOutZero()
{
  int i;
  if (m_MassIn_up)
    memset(m_MassIn_up, 0, sizeof(double)*m_ny);
  if (m_MassIn_left)
    memset(m_MassIn_left, 0, sizeof(double)*m_nx);
  if (m_MassOut_right)
    memset(m_MassOut_right, 0, sizeof(double)*m_n_gridsBdyRight);
  if (m_MassOut_down)
    memset(m_MassOut_down, 0, sizeof(double)*m_n_gridsBdyDown);
}

void Comp::passBdyMassOut(Comp *bdyRight, Comp *bdyDown)
{
  assert(!bdyRight || m_n_gridsBdyRight==bdyRight->m_nx);
  assert(!bdyDown || m_n_gridsBdyDown==bdyDown->m_ny);

  if (bdyRight)
    memcpy(bdyRight->m_MassIn_left, m_MassOut_right, sizeof(double)*m_n_gridsBdyRight);
  if (bdyDown)
    memcpy(bdyDown->m_MassIn_up, m_MassOut_down, sizeof(double)*m_n_gridsBdyDown);
}

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

  return area;

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

  return volume;
}

/* compute this compartment's total area to a certain direction defined below
   direction: [0] = up; [1] = left; [2] = right; [3] = down */
double Comp::compTotalArea(int direction)
{
  double area, r1, r2, pi_alpha_360;

  switch (m_coord_sys) {

  case Cartesian :

    if (direction==0 || direction==3)
      area = m_y_length * m_dz_dtheta;
    else if (direction==1 || direction==2)
      area = m_x_length * m_dz_dtheta;
    else
      SayBye("Option not valid");
    break;

  case Cylindrical :

    pi_alpha_360 = M_PI * m_dz_dtheta / 360;

    if (direction==0 || direction==3)
      area = pi_alpha_360 * m_y_length * m_y_length;
    else if (direction==1)
      area = 0;
    else if (direction==2)
      area = m_x_length * pi_alpha_360 * 2 * m_y_length;
    else
      SayBye("Option not valid");
    break;

  default :
    SayBye("Coordinate system not implemented");
    break;
  }

  return area;
}

/* compute this compartment's total volume */
double Comp::compTotalVolume()
{
  double volume;
  SayBye("not implemented yet");
  return volume;
}


/* compute mass transfer between gridThhis and neighbouring grids to the right
   whose meshing does not match exactly to gridThiis, e.g. 
   gridThiis may interface with multiple grids to the right */
double Comp::compMassIrregGridsRight(Grid gridThiis, double conc_this)
{
  int i;
  bool bDone = false;
  double x1_this, x2_this, currentX, nextX, x_length, z_length, thd, deriv_this, deriv_other, area;
  double massIntoThis, mass, conc_other, flux;
  Grid *gridOther = NULL;
  
  //  thd = gridThiis.m_dx * 1e-3; // to compare whether two real numbers are the same

  massIntoThis = 0;

  z_length = compInterArea(gridThiis, 2) / gridThiis.m_dx; // interfacial z_length

  for (i=0; i<m_n_gridsBdyRight; i++) {

    x1_this = gridThiis.m_x_coord;     // x coordiantes of this
    x2_this = x1_this + gridThiis.m_dx;//

    currentX = m_gridsBdyRight[i].m_x_coord;    // x coordinates of the neighbouring grid
    nextX = currentX + m_gridsBdyRight[i].m_dx; //

    // to compare whether two real numbers are the same
    thd = gridThiis.m_dx<m_gridsBdyDown[i].m_dx ? gridThiis.m_dx : m_gridsBdyDown[i].m_dx;
    thd *= 1e-3;

    if ( x1_this > currentX-thd ) {

      if ( x2_this < nextX+thd ) { // gridThiis is contained between currentX & nextX
	x_length = x2_this - x1_this;
	bDone = true;
      }
      else { // gridThiis extends beyond nextX
	x_length = nextX - x1_this;
      }

      gridOther = &m_gridsBdyRight[i];
      conc_other = gridOther->m_concChem;
     
      flux = gridThiis.compFlux( gridOther, conc_this, conc_other, 
				 gridThiis.m_dy/2, gridOther->m_dy/2, &deriv_this, &deriv_other);
      area = z_length * x_length;
      mass = area * flux;

      massIntoThis += mass;
      m_MassOut_right[i] += -mass; // flux into gridThiis is negative flux into the neighbour

      if (bDone)
	break;
      else {
	gridThiis.m_x_coord = nextX;
	gridThiis.m_dx -= x_length;
      }
    } // if

  } // for i

  return massIntoThis;
}


/* compute mass transfer between gridThhis and neighbouring grids downward
   whose meshing does not match exactly to gridTThis, e.g. 
   gridTThis may interface with multiple grids downward */
double Comp::compMassIrregGridsDown(Grid gridThiis, double conc_this)
{
  int i;
  bool bDone = false;
  double y1_this, y2_this, currentY, nextY, y_length, z_length, thd, deriv_this, deriv_other, area;
  double massIntoThis, mass, conc_other, flux;
  Grid *gridOther = NULL;
  
  //  thd = gridThiis.m_dy * 1e-3; // to compare whether two real numbers are the same
  massIntoThis = 0;

  for (i=0; i<m_n_gridsBdyDown; i++) {

    y1_this = gridThiis.m_y_coord;      // y coordinates of this grid
    y2_this = y1_this + gridThiis.m_dy; //

    currentY = m_gridsBdyDown[i].m_y_coord;    // y coordinates of the neighbouring grid
    nextY = currentY + m_gridsBdyDown[i].m_dy; //

    // to compare whether two real numbers are the same
    thd = gridThiis.m_dy<m_gridsBdyDown[i].m_dy ? gridThiis.m_dy : m_gridsBdyDown[i].m_dy;
    thd *= 1e-3;

    if ( y1_this > currentY-thd ) {

      if ( y2_this < nextY+thd ) { // gridThiis is contained between currentY & nextY
	y_length = y2_this - y1_this;
	bDone = true;
      }
      else { // gridThiis extends beyond nextY
	y_length = nextY - y1_this;
      }

      gridOther = &m_gridsBdyDown[i];
      conc_other = gridOther->m_concChem;
     
      flux = gridThiis.compFlux( gridOther, conc_this, conc_other, 
				 gridThiis.m_dx/2, gridOther->m_dx/2, &deriv_this, &deriv_other);

      gridThiis.m_dy = y_length; // this is needed to calculate the correct interfacial area, especially for cylindrical coordinate
      area = compInterArea(gridThiis, 3); // interfacial area
      mass = area * flux;

      massIntoThis += mass;
      m_MassOut_down[i] += -mass; // flux into gridThiis is negative flux into the neighbour

      if (bDone)
	break;
      else {
	gridThiis.m_y_coord = nextY;
	gridThiis.m_dy = y2_this - y1_this - y_length;
      }
    } // if

  } // for i

  return massIntoThis;
}


/* functions for computing the right-hand size of the odes */

void Comp::compODE_dydt(double t, const double y[], double f[])
{
  int i, rc;
	
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
	  flux = gridThiis->compFlux( &m_gridSink, conc_this, 0, gridThiis->m_dy/2, 0, &deriv_this, &deriv_other);	  
	  mass_transfer_rate += area * flux;
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
	  flux = gridThiis->compFlux( &m_gridSink, conc_this, 0, gridThiis->m_dy/2, 0, &deriv_this, &deriv_other);	  
	  mass_transfer_rate += area * flux;
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
	  mass = compMassIrregGridsRight(*gridThiis, conc_this);
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

      area = compInterArea(*gridThiis, 3); // interfacial area

      if ( i==m_nx-1 ) { // bottom layer

	switch (m_BdyCond_down) {
	case ZeroFlux : 
	  break;
	case ZeroConc :
	  flux = gridThiis->compFlux( &m_gridSink, conc_this, 0, gridThiis->m_dx/2, 0, &deriv_this, &deriv_other);	  
	  mass_transfer_rate += area * flux;
	  break;
	case FromOther :
	  mass = compMassIrregGridsDown(*gridThiis, conc_this);
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
  printf("# of grids: [x] %d, [y] %d\n", m_nx, m_ny);

  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right	

      idx = i*m_ny + j;			
      if ( !strcmp(m_grids[idx].m_name, "LP") )
	printf("L ");
      else if ( !strcmp(m_grids[idx].m_name, "CC") )
	printf("C ");
      else if ( !strcmp(m_grids[idx].m_name, "VE") )
	printf("V ");
      else
	SayBye ("subtype name unknown" ); 
				
    } // for j
    printf("\n");
  } // for i
  fflush(stdout);
}


double Comp::getAmount()
{
  assert( m_grids );

  int i, j, idx;
  double volume, amount;
	
  amount = .0;
  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right
       idx = i*m_ny + j;

       volume = compVolume(m_grids[idx]);
       amount += m_grids[idx].getConcChem() * volume;
     }
  }

  return amount;
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

void Comp::saveCoord(const char fn_x[], const char fn_y[], const char fn_suffix[])
{
  assert( m_grids );

  FILE *file_x, *file_y;
  char fn1[1024], fn2[1024];
  int i, j, idx;

  // save grids
  strcpy(fn1, fn_x); strcat(fn1, fn_suffix);
  strcpy(fn2, fn_y); strcat(fn2, fn_suffix);

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
