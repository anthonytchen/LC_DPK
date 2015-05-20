#include "stdafx.h"
#include "StraCorn.h"

/* Structure and definition for parallel computing */
#define NTHREADS 1 // number of threads for parallel computing
struct pthread_struct {
	StraCorn *sc_obj;
	double t;
	const double *y;
	double *f;
	int x_start, x_end, y_start, y_end;
};
/* ------------------- */



/*
  Geometric dimensions:
    g: vertical size of lipid layer (typical: 0.075 micron)
    d: lateral size of corneocyte (typical: 40 micron)
    s: lateral space between two corneocytes (typical: 0.075 micron)
    t: vertical size of corneocyte (typical: 0.8 micron)
    dz: depth (perpendicular to x-y domain) of the simulation domain
 */
void StraCorn::Init(double g, double d, double s, double t, double dz,
		    int n_layer_x, int n_layer_y, double offset_y, int bdy_cond)
{	

  /* set up some constant values */

  // density of lipid, keratin and water
  m_rou_lipid = 1e3; // kg m^{-3}
  m_rou_keratin = 1.2e3; // kg m^{-3}
  m_rou_water = 1e3; // kg m^{-3}

  // dimension related; c.f. Readme.docx for more details
  m_w = 8.0; // offset ratio, 8.0
  m_nx_grids_lipid = 2; // # of x-grids for lipid layer, 2
  m_nx_grids_cc = 4; // # of x-grids for corneocyte layer, 4
  m_ny_grids_lipid = 2; // # of y-grids for lipid layer, 2
  m_ny_grids_cc_dn = 2; // # of y-grids for dn-part of the offset corneocyte layer, 2

  m_T = 309; // temperature (Kelvin)
  m_eta = 7.1E-4; // water viscosity at above temperature (Pa s),

  m_boundary_cond = bdy_cond; // boundary condition for left/right; 

  /* ---- */

  m_grids = NULL;
  m_ode_Jacobian = NULL;

  assert( s < d ); // inter-corneocyte gap must be less than the corneocyte width
  m_geom_g = g; 
  m_geom_d = d; 
  m_geom_s = s; 
  m_geom_t = t; 

  m_geom_dm = m_w*(m_geom_d-m_geom_s)/(1+m_w);
  m_geom_dn = m_geom_d-m_geom_s-m_geom_dm;
  m_dz = dz;
	
  // Vertical direction, lipid layer is at both top and bottom of the stratum corneum
  m_nx = (m_nx_grids_lipid+m_nx_grids_cc)*n_layer_x + m_nx_grids_lipid; 	
  // Lateral direction, [dh] [s] [dm] [s], here d=dh+dm+s, w=dm/dh
  m_ny = (int) ( m_ny_grids_lipid*2 + m_ny_grids_cc_dn + m_ny_grids_cc_dn*m_w ) * n_layer_y;

  m_x_length = n_layer_x*(g+t)+g;
  m_y_length = n_layer_y*(d+s);

  m_V_mortar = ( g*(d+s)+t*s ) * dz;	
  m_V_brick = d*t * dz;
  m_V_all = m_V_mortar + m_V_brick;

  m_offset_y =  offset_y;

  // setup the array for 1D concentration and coordinates
  m_conc1D = new double [n_layer_x];
  m_coord1D = new double [n_layer_x];
  for ( int i = 0; i < n_layer_x-1; i ++ )
    m_coord1D[i] = 0.5*( i*(g+t) + (i+1)*(g+t) );
  m_coord1D[n_layer_x-1] = 0.5*( (n_layer_x-1)*(g+t) + n_layer_x*(g+t)+g );
  m_n_layer_x = n_layer_x;
}

void StraCorn::Release()
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
  delete [] m_conc1D;
  delete [] m_coord1D;

  m_gridBdyUp.Release();
  m_gridBdyDown.Release();
  m_gridBdyLeft.Release();
  m_gridBdyRight.Release();
}


void StraCorn::createGrids(double MW, double Kow, double water_frac, double conc_vehicle, double diffu_vehicle)
{
  bool bOffset = false;
  int i, j, idx, idx_x, idx_y, idx_y_offset, cc_subtype_offset, gsl_errno;
  int cc_subtype = 0; // 0 = d_n; 1 = s; 2 = d_m, 3 = s;
  double dx_lipid, dx_cc, dy_lipid, dy_cc, dy_offset, dy_dm, coord_x, coord_y;

  // initialise boundary grids
  m_gridBdyUp.Init("SC", conc_vehicle, Kow, 0, 0, m_dz, diffu_vehicle); // infinite source, can be updated to diminishing vehicle by using updateBoundary()
  m_gridBdyDown.Init("SK",  0,         Kow, 0, 0, m_dz); // infinite sink, can be updated to underlying viable epidermis by using updateBoundary()
  m_gridBdyLeft.Init("SK",  0,         Kow, 0, 0, m_dz); // infinite sink
  m_gridBdyRight.Init("SK", 0,         Kow, 0, 0, m_dz); // infinite sink

	
  dx_lipid = m_geom_g/m_nx_grids_lipid;
  dx_cc = m_geom_t/m_nx_grids_cc;
  dy_lipid = m_geom_s/m_ny_grids_lipid;
  dy_cc = m_geom_dn/m_ny_grids_cc_dn;
	
  m_grids = new Grid[m_nx*m_ny]; // organised in row dominant

  // work out the starting point from given offset
  
  double len = m_geom_d+m_geom_s, len_vec[4];
  len_vec[0]=m_geom_dn; len_vec[1]=m_geom_s; len_vec[2]=m_geom_dm; len_vec[3]=m_geom_s;

  while (m_offset_y > len) m_offset_y -= len;
  len = 0;
  for ( i = 0; i < 4; i ++ ){
    len += len_vec[i];
    if (m_offset_y < len) break;
  }
  cc_subtype_offset = i;
  len = m_offset_y - len + len_vec[i];

    
  switch (cc_subtype_offset) {
  case 0 :
  case 2 : // corneocyte width 
    idx_y_offset = (int) floor( len / dy_cc );
    dy_offset = dy_cc;
    break;
  case 1 :
  case 3 : // lipid width
    idx_y_offset = (int) floor( len / dy_lipid );
    dy_offset = dy_lipid;
    break;
  default :
    gsl_error ("subtype name unknown", __FILE__, __LINE__, gsl_errno);    
  }

  struct Point current_point;
  setPoint(current_point, 0, 0, dx_lipid, dy_offset, "LP", "LP"); // starting from lipid on the top layer
    
  idx_x = 0; idx_y = idx_y_offset;
  coord_x = coord_y = 0;
  cc_subtype = cc_subtype_offset;

  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right
			
      idx = i*m_ny + j;

      // assign type
      if ( !strcmp(current_point.x_type, "LP") || !strcmp(current_point.y_type, "LP") ) { 
	// entire lipid layer (1st strcmp) or lateral lipid between two coreneocytes
	m_grids[idx].Init("LP", MW, water_frac, water_frac, m_V_mortar, m_V_brick, m_V_all,
			  m_rou_lipid, m_rou_keratin, m_rou_water, m_T, m_eta, Kow, 
			  current_point.x_coord, current_point.y_coord, current_point.dx, current_point.dy, m_dz);
      } else {
	m_grids[idx].Init("CC", MW, water_frac, water_frac, m_V_mortar, m_V_brick, m_V_all,
			  m_rou_lipid, m_rou_keratin, m_rou_water, m_T, m_eta, Kow, 
			  current_point.x_coord, current_point.y_coord, current_point.dx, current_point.dy, m_dz);
      }

      // update current_point
      if (j==m_ny-1) { // last element in the lateral direction, move down
				
	idx_x ++;

	if ( !strcmp(current_point.x_type, "LP") ) {
	  coord_x += dx_lipid;
	  if ( idx_x == m_nx_grids_lipid ) {
	    switch (cc_subtype_offset) {
	    case 0 :
	    case 2 :
	      setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "CC");
	      break;
	    case 1 :
	      if (!bOffset)   setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "CC");
	      else            setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "LP");
	      break;
	    case 3 :
	      if (!bOffset)   setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "LP");
	      else            setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "CC");
	      break;
	    default :
	      break;
	    }	      
	    idx_x = 0;
	  } else {
	    setPoint(current_point, coord_x, 0, dx_lipid, dy_offset, "LP", "LP");
	  }
	} else if ( !strcmp(current_point.x_type, "CC") ) {
	  coord_x += dx_cc;
	  if ( idx_x == m_nx_grids_cc ) {
	    setPoint(current_point, coord_x, 0, dx_lipid, dy_offset, "LP", "LP");
	    idx_x = 0;
	    bOffset = !bOffset;
	  } else {
	    switch (cc_subtype_offset) {
	    case 0 :
	    case 2 :
	      setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "CC");
	      break;
	    case 1 :
	      if (!bOffset)   setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "CC");
	      else            setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "LP");
	      break;
	    case 3 :
	      if (!bOffset)   setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "LP");
	      else            setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "CC");
	      break;
	    default :
	      break;
	    }
	  }
	}
	coord_y = 0;
	idx_y = idx_y_offset;
	cc_subtype = cc_subtype_offset;

      } else { // not the last element in the lateral direction, thus move to the right
			
	idx_y ++;
				
	if ( !strcmp(current_point.x_type, "LP") ) { // current row is lipid					
	  switch (cc_subtype) {
	  case 0 : // now within dh
	    coord_y += dy_cc;
	    if (idx_y==m_ny_grids_cc_dn){
	      setPoint(current_point, coord_x, coord_y, dx_lipid, dy_lipid, "LP", "LP");
	      idx_y = 0;
	      cc_subtype ++;
	    } else {
	      setPoint(current_point, coord_x, coord_y, dx_lipid, dy_cc, "LP", "LP");
	    }
	    break;
	  case 1 : // now within s
	    coord_y += dy_lipid;
	    if (idx_y==m_ny_grids_lipid){
	      setPoint(current_point, coord_x, coord_y, dx_lipid, dy_cc, "LP", "LP");
	      idx_y = 0;
	      cc_subtype ++;
	    } else {
	      setPoint(current_point, coord_x, coord_y, dx_lipid, dy_lipid, "LP", "LP");
	    }
	    break;
	  case 2 : // now wtihin dm
	    coord_y += dy_cc;
	    if (idx_y==m_ny_grids_cc_dn*m_w){
	      setPoint(current_point, coord_x, coord_y, dx_lipid, dy_lipid, "LP", "LP");
	      idx_y = 0;
	      cc_subtype ++;
	    } else {
	      setPoint(current_point, coord_x, coord_y, dx_lipid, dy_cc, "LP", "LP");
	    }
	    break;
	  case 3 : // now within the 2nd s
	    coord_y += dy_lipid;
	    if (idx_y==m_ny_grids_lipid){
	      setPoint(current_point, coord_x, coord_y, dx_lipid, dy_cc, "LP", "LP");
	      idx_y = 0;
	      cc_subtype = 0;
	    } else {
	      setPoint(current_point, coord_x, coord_y, dx_lipid, dy_lipid, "LP", "LP");
	    }
	    break;
	  default :
	    gsl_error("cc_subtype not implemented", __FILE__, __LINE__, gsl_errno);
	    exit(-1);
	  } // switch-case
	} else if ( !strcmp(current_point.x_type, "CC") ) { // current row is corneocyte
	  switch (cc_subtype) {
	  case 0 : // now within dh
	    coord_y += dy_cc;
	    if (idx_y==m_ny_grids_cc_dn){
	      if (bOffset)
		setPoint(current_point, coord_x, coord_y, dx_cc, dy_lipid, "CC", "LP");
	      else
		setPoint(current_point, coord_x, coord_y, dx_cc, dy_lipid, "CC", "CC");
	      idx_y = 0;
	      cc_subtype ++;
	    } else {
	      setPoint(current_point, coord_x, coord_y, dx_cc, dy_cc, "CC", "CC");
	    }
	    break;
	  case 1 : // now within s
	    coord_y += dy_lipid;
	    if (idx_y==m_ny_grids_lipid){
	      setPoint(current_point, coord_x, coord_y, dx_cc, dy_cc, "CC", "CC");
	      idx_y = 0;
	      cc_subtype ++;
	    } else {
	      if (bOffset)
		setPoint(current_point, coord_x, coord_y, dx_cc, dy_lipid, "CC", "LP");
	      else
		setPoint(current_point, coord_x, coord_y, dx_cc, dy_lipid, "CC", "CC");
	    }
	    break;
	  case 2 : // now wtihin dm
	    coord_y += dy_cc;
	    if (idx_y==m_ny_grids_cc_dn*m_w){
	      if (bOffset)
		setPoint(current_point, coord_x, coord_y, dx_cc, dy_lipid, "CC", "CC");
	      else
		setPoint(current_point, coord_x, coord_y, dx_cc, dy_lipid, "CC", "LP");
	      idx_y = 0;
	      cc_subtype ++;
	    } else {
	      setPoint(current_point, coord_x, coord_y, dx_cc, dy_cc, "CC", "CC");
	    }
	    break;
	  case 3 : // now within the 2nd s
	    coord_y += dy_lipid;
	    if (idx_y==m_ny_grids_lipid){
	      setPoint(current_point, coord_x, coord_y, dx_cc, dy_cc, "CC", "CC");
	      idx_y = 0;
	      cc_subtype = 0;
	    } else {
	      if (bOffset)
		setPoint(current_point, coord_x, coord_y, dx_cc, dy_lipid, "CC", "CC");
	      else
		setPoint(current_point, coord_x, coord_y, dx_cc, dy_lipid, "CC", "LP");
	    }
	    break;
	  default :
	    printf("error: cc_subtype not implemented\n");
	    exit(0);
	  } // switch-case
	} // if-else-if
      } // if (j==m_ny-1)
      
    } // for j
  } // for i

}

void StraCorn::updateBoundary(Grid* up, Grid* down, Grid* left, Grid* right)
{
  if (up!=NULL) m_gridBdyUp.set(up);
  if (down!=NULL) m_gridBdyDown.set(down);
  if (left!=NULL) m_gridBdyLeft.set(left);
  if (right!=NULL) m_gridBdyRight.set(right);
}

/* functions for computing the right-hand size of the odes */

// todo: remove repeated calculation of mass transfer
void StraCorn::compODE_dydt (double t, const double y[], double f[])
{
  int i, rc;
	
  if (NTHREADS==1) {
    compODE_dydt_block (t, y, f, 0, m_nx, 0, m_ny);
  } else {		
    struct pthread_struct p[NTHREADS];
    pthread_t threads[NTHREADS];
			
    for ( i=0; i < NTHREADS; i++ ) {
      p[i].sc_obj = this;
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
void* StraCorn::static_compODE_dydt_block_threads(void *paras)
{
  struct pthread_struct p = *((struct pthread_struct *) paras);
  p.sc_obj->compODE_dydt_block (p.t, p.y, p.f, p.x_start, p.x_end, p.y_start, p.y_end);
}

// the actual funtion to calculate dy/dy
void StraCorn::compODE_dydt_block (double t, const double y[], double f[], 
				   int idx_x_start, int idx_x_end, int idx_y_start, int idx_y_end)
{
  int i, j, idx_this, idx_other, dim;
  double flux, mass, mass_transfer_rate, conc_this, conc_other, deriv_this, deriv_other,
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
  if ( idx_x_start == 0 ) m_mass_in = 0; // re-set mass transferred into SC, ready for calculation
  if ( idx_x_end == m_nx ) m_mass_out = 0; // re-set mass transferred out of SC, ready for calculation

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
      mass = gridThiis->m_dy*gridThiis->m_dz * flux;

      if ( i==0 ) m_mass_in += mass;
      mass_transfer_rate +=  mass;
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
      mass = gridThiis->m_dy*gridThiis->m_dz * flux;

      if ( i==m_nx-1 ) m_mass_out += -mass;
      mass_transfer_rate += mass;
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

	
/*  ++++++++++++++++++++++++++++++++++
	I/O functions
	++++++++++++++++++++++++++++++++++ */

void StraCorn::displayGrids()
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

void StraCorn::getGridsConc(double *fGridsConc, int dim)
{
  // Return concentration at the grids in fGridsConc
  assert( m_grids && fGridsConc && dim==m_nx*m_ny);

  int i, j, idx;
	
  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right
       idx = i*m_ny + j;
       fGridsConc[idx] = m_grids[idx].getConcChem();
     }
  } // for i
}

void StraCorn::comp1DConc()
{
  int i, j, idx, nx_grids_lipid, nx_grids_cc, idx_conc1D;
  double conc_sum, area, area_sum;

  assert( m_n_layer_x > 1 );

  idx_conc1D = 0;

  conc_sum = 0; area_sum = 0;
  nx_grids_lipid = m_nx_grids_lipid;
  nx_grids_cc = m_nx_grids_cc;

  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right

       idx = i*m_ny + j;

       area = m_grids[idx].m_dx * m_grids[idx].m_dy;
       conc_sum += m_grids[idx].getConcChem() * area;
       area_sum += area;

    } // for j

    if ( nx_grids_lipid > 0 )
      nx_grids_lipid --;
    else
      nx_grids_cc --;

    if ( !nx_grids_lipid && !nx_grids_cc ) { // finished calculation for 1 lipid and 1 corneocyte layer
      m_conc1D[idx_conc1D] = conc_sum / area_sum;
      idx_conc1D ++;
      conc_sum = 0; area_sum = 0;
      if ( idx_conc1D == m_n_layer_x-1 )
	nx_grids_lipid = m_nx_grids_lipid * 2; // also count the botton lipid layer
      else 
	nx_grids_lipid = m_nx_grids_lipid;
      nx_grids_cc = m_nx_grids_cc;
    }

  } // for i

}

void StraCorn::saveGrids(bool b_1st_time, const char fn[])
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

void StraCorn::getXCoord(double *coord_x, int dim)
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

void StraCorn::getYCoord(double *coord_y, int dim)
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

void StraCorn::saveCoord(const char fn_x[], const char fn_y[])
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
