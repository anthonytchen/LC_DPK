#include "stdafx.h"
#include "Skin.h"

#define NTHREADS 8

struct pthread_struct {
	Skin *skin_obj;
	double t;
	const double *y;
	double *f;
	int x_start, x_end, y_start, y_end;
};

void Skin::Init(Chemical chemSolute, double conc_vehicle, double diffu_vehicle,	double partition_vehicle, 
		double dx_vehicle, double area_vehicle, 
		int n_layer_x_sc, int n_layer_y_sc, int n_grids_x_ve, int n_grids_x_de, double offset_y_sc,
		bool bInfSrc)
{

  m_dz = 0.01; // fixing dz, the dimension perpendicular to x-y domain

  m_boundary_cond = 1; // Boundary condition, [0] - zero flux at left/right sides
	               //    [1] - mirror flux at left/right sides (i.e. periodic boundary condition)

  /* set up stratum corneum using fixed geometry */

  double g, d, s, t, water_frac;
  g=.075e-6; d=40e-6; s=0.075e-6; t=0.8e-6;
  water_frac = 0.55; // mass fraction of water in stratum corneum

  m_StraCorn.Init(g, d, s, t, m_dz, n_layer_x_sc, n_layer_y_sc, offset_y_sc, m_boundary_cond);
  m_StraCorn.createGrids(chemSolute.m_mw, chemSolute.m_K_ow, water_frac, conc_vehicle, diffu_vehicle);


  /* set up viable epidermis using fixed geometry */

  double x_len_sc, x_len_ve, y_len_ve;
  x_len_ve = 20e-6;
  x_len_sc = n_layer_x_sc*(g+t) + g;
  y_len_ve = n_layer_y_sc*(d+s); // depends on the setup for stratum corneum

  m_ViaEpd.Init(x_len_ve, y_len_ve, m_dz, n_grids_x_ve);
  m_ViaEpd.createGrids(chemSolute.m_mw, chemSolute.m_K_ow, chemSolute.m_pKa, chemSolute.m_acid_base, x_len_sc);

  /* set up dermis using fixed geometry */

  double x_len_de, y_len_de;
  x_len_de = 900e-6;
  y_len_de = y_len_ve; // depends on the setup for ve (thus also on stratum corneum)

  m_Dermis.Init(x_len_de, y_len_de, m_dz, n_grids_x_de);
  m_Dermis.createGrids(chemSolute.m_mw, chemSolute.m_K_ow, chemSolute.m_pKa, chemSolute.m_acid_base, x_len_sc+x_len_ve);

  /* set up blood compartment, then set up the blood properties in dermis */
  m_Blood.Init(m_Dermis.m_grids->m_ve_fu, 70, 'M');
  double par_de2blood = 1/pow(10, 0.04); // log_P blood:skin is 0.04 for nicotine
  m_Dermis.InitDermisBlood(m_Blood.m_flow_capil, m_Blood.m_f_unbound, par_de2blood);

  /* set up vehicle using fixed geometry */
  m_Vehicle_area = area_vehicle;
  m_gridVehicle.Init("SC", conc_vehicle, chemSolute.m_K_ow, dx_vehicle, y_len_ve, m_dz, diffu_vehicle, partition_vehicle);
  m_gridSink.Init("SK", 0, chemSolute.m_K_ow, 0, 0, m_dz);
  m_bInfSrc = bInfSrc; // whether the vehicle is a infinite source

  // setup the dimensions
  m_dim_vh = 1;
  m_dim_sc = m_StraCorn.m_nx*m_StraCorn.m_ny;
  m_dim_ve = m_ViaEpd.m_nx*m_ViaEpd.m_ny;
  m_dim_de = m_Dermis.m_nx*m_Dermis.m_ny;
  m_dim_bd = 1;
  m_dim_all =  m_dim_vh + m_dim_sc + m_dim_ve + m_dim_de +  m_dim_bd;

}

void Skin::Release(void)
{
  m_StraCorn.Release();
  m_ViaEpd.Release();
  m_Dermis.Release();
  m_Blood.Release();
}


int Skin::static_cvODE (double t, N_Vector y, N_Vector dydt, void *paras)
{
  double *p_y, *p_dydt;
	
  p_y = NV_DATA_S(y);
  p_dydt = NV_DATA_S(dydt);
  return ((Skin*)paras)->compODE_dydt(t, p_y, p_dydt);
}


int Skin::compODE_dydt (double t, const double y[], double f[])
{
  int i, j;
  Grid gridUp, gridDown;

  /* y and f are organised as:
     vehicle (1 grid), sc (row dominant), ve (row dominant) */

  /* compute for stratum corneum */
  gridUp.set(&m_gridVehicle);  gridUp.m_concChem = y[0];
  gridDown.set(m_ViaEpd.m_grids); gridDown.m_concChem = y[m_dim_vh+m_dim_sc];
  m_StraCorn.updateBoundary(&gridUp, &gridDown, NULL, NULL);  // update top (vehicle) and bottom (ve) boundary
  m_StraCorn.compODE_dydt(t, y+m_dim_vh, f+m_dim_vh);

  /* compute for vehicle */
  //     this is after the computation for SC because mass transfer from vehicle to SC is calculated in the SC
  if (m_bInfSrc)
    f[0] = 0; // infinite source, concentration does not change
  else
    f[0] = -m_StraCorn.m_mass_in/(m_gridVehicle.m_dx*m_gridVehicle.m_dy*m_gridVehicle.m_dz);

  /* compute for viable epidermis */
  gridDown.set(m_Dermis.m_grids); gridDown.m_concChem = y[m_dim_vh+m_dim_sc+m_dim_ve];
  m_ViaEpd.updateBoundary(NULL, &gridDown, NULL, NULL, m_StraCorn.m_mass_out); // update top (sc) and bottom (de) boundary
  m_ViaEpd.compODE_dydt(t, y+m_dim_vh+m_dim_sc, f+m_dim_vh+m_dim_sc);

  /* compute for dermis */
  gridDown.set(&m_gridSink);
  m_Dermis.updateBoundary(NULL, &gridDown, NULL, NULL, m_ViaEpd.m_mass_out); // update top (ve) and bottom (sink) boundary
  m_Dermis.updateBlood(y[m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de]);
  m_Dermis.compODE_dydt(t, y+m_dim_vh+m_dim_sc+m_dim_ve, f+m_dim_vh+m_dim_sc+m_dim_ve);

  /* compute for blood */

  // simulation is for a small skin area, but needs to multiple
  //  the mass transport due to blood flow by the actual topical application area
  double factor = m_Vehicle_area / (m_Dermis.m_y_length*m_Dermis.m_dz);

  m_Blood.updateMassInOutDermis(m_Dermis.m_mass_into_dermis, m_Dermis.m_mass_outof_dermis, factor);
  m_Blood.compODE_dydt(t, y+m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de, f+m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de);

  return GSL_SUCCESS;	
}

/* Diffusion using method of lines (MoL) */
void Skin::diffuseMoL(double t_start, double t_end)
{		
  bool bJacobianRequired = false;
  int i, j, gsl_status, dim, dim_vh, dim_sc, dim_ve, dim_de, dim_bld, flag;
  double reltol, abstol, t;
  double *y = NULL;
  N_Vector y0;
  void *cvode_mem = NULL;

  dim_vh = 1;
  dim_sc = m_StraCorn.m_nx*m_StraCorn.m_ny;
  dim_ve = m_ViaEpd.m_nx*m_ViaEpd.m_ny;
  dim_de = m_Dermis.m_nx*m_Dermis.m_ny;
  dim_bld = 1;
  dim = dim_vh + dim_sc + dim_ve + dim_de + dim_bld;

  /* get current concentration, and set as initial conditions */

  y = new double[dim];
  
  // from vehicle
  y[0] = m_gridVehicle.m_concChem; 

  // from stratum corneum
  for ( i=0; i<m_StraCorn.m_nx; i++ ) // x direction up to down
    for ( j=0; j<m_StraCorn.m_ny; j++ ) // y direction left to right	
      y[dim_vh + i*m_StraCorn.m_ny + j] = m_StraCorn.m_grids[i*m_StraCorn.m_ny + j].m_concChem;


  // from viable epidermis
  for ( i=0; i<m_ViaEpd.m_nx; i++ ) // x direction up to down
    for ( j=0; j<m_ViaEpd.m_ny; j++ ) // y direction left to right	
      y[dim_vh+dim_sc + i*m_ViaEpd.m_ny + j] = m_ViaEpd.m_grids[i*m_ViaEpd.m_ny + j].m_concChem;

  // from dermis
  for ( i=0; i<m_Dermis.m_nx; i++ ) // x direction up to down
    for ( j=0; j<m_Dermis.m_ny; j++ ) // y direction left to right	
      y[dim_vh+dim_sc+dim_ve + i*m_Dermis.m_ny + j] = m_Dermis.m_grids[i*m_Dermis.m_ny + j].m_concChem;

  // from blood
  y[dim_vh+dim_sc+dim_ve+dim_de] = m_Blood.m_concChem;

  /* ------------- */

  y0 = N_VMake_Serial(dim, y);
	
  // Call CVodeCreate to create the solver memory and specify the 
  //	 Backward Differentiation Formula and the use of a Newton iteration
  // cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
		
  // Call CVodeInit to initialize the integrator memory and specify the
  //	user's right hand side function in y'=f(t,y), the inital time t_start, and
  //	the initial condition y0.
  CVodeInit(cvode_mem, static_cvODE, t_start, y0);

  // Call CVodeSStolerances to specify the scalar relative tolerance
  //	 and scalar absolute tolerance.
  reltol=1e-6; abstol=1e-10;
  CVodeSStolerances(cvode_mem, reltol, abstol);
	
  // Set the pointer to user-defined data
  CVodeSetUserData(cvode_mem, this);
  CVodeSetMaxNumSteps(cvode_mem, 10000*5);
	
  // Call CVBand to specify the CVBAND band linear solver
  //	m_ny is the bandwidth of the banded Jacobian matrix.
  //	Then setup the Jacobian function
  //CVBand(cvode_mem, dim, m_ny, m_ny);
  //CVDlsSetBandJacFn(cvode_mem, static_cvJacobian);
  
  CVDense(cvode_mem, dim);
  // CVDlsSetDenseJacFn(cvode_mem, static_cvJacobian);
  
  // prepare for the memory
  // m_gsl_ode_Jacobian = new double [dim*dim];
  // memset(m_gsl_ode_Jacobian, 0, sizeof(double)*dim*dim);
  
  CVode(cvode_mem, t_end, y0, &t, CV_NORMAL);
	
	
  /* Extract calculated values */
		  
  // for vehicle
  m_gridVehicle.m_concChem = NV_Ith_S(y0, 0); 

  // for stratum corneum
  for ( i=0; i<m_StraCorn.m_nx; i++ ) // x direction up to down
    for ( j=0; j<m_StraCorn.m_ny; j++ ) // y direction left to right	
      m_StraCorn.m_grids[i*m_StraCorn.m_ny + j].m_concChem = NV_Ith_S(y0, dim_vh + i*m_StraCorn.m_ny + j);


  // for viable epidermis
  for ( i=0; i<m_ViaEpd.m_nx; i++ ) // x direction up to down
    for ( j=0; j<m_ViaEpd.m_ny; j++ ) // y direction left to right	
      m_ViaEpd.m_grids[i*m_ViaEpd.m_ny + j].m_concChem = NV_Ith_S(y0, dim_vh+dim_sc + i*m_ViaEpd.m_ny + j);

  // for viable epidermis
  for ( i=0; i<m_Dermis.m_nx; i++ ) // x direction up to down
    for ( j=0; j<m_Dermis.m_ny; j++ ) // y direction left to right	
      m_Dermis.m_grids[i*m_Dermis.m_ny + j].m_concChem = NV_Ith_S(y0, dim_vh+dim_sc+dim_ve + i*m_Dermis.m_ny + j);

  // for blood
  m_Blood.m_concChem = NV_Ith_S(y0, dim_vh+dim_sc+dim_ve+dim_de);

  /* ------------------------- */
  

  N_VDestroy_Serial(y0);  
  CVodeFree(&cvode_mem);
  //  delete [] m_gsl_ode_Jacobian;
  delete [] y;
}

 // reset vehicle concentration, partition coefficient, diffusivity
void Skin::resetVehicle(double concChem, double partition_coef, double diffu_coef)
{
  m_gridVehicle.Init("SC", concChem, m_gridVehicle.m_K_ow, m_gridVehicle.m_dx, m_gridVehicle.m_dy, m_gridVehicle.m_dz, 
		     diffu_coef, partition_coef);
}

void Skin::removeVehicle()
{
  resetVehicle(m_gridVehicle.m_concChem, 1e10, 1e-100);
}

// Compute flux into stratum corneum
double Skin::compFlux_2sc()
{
  int j;
  double flux, conc_this, conc_other, deriv_this, deriv_other, mass_transfer_rate, area, total_area;
  Grid *gridThiis = NULL;

  mass_transfer_rate = total_area = 0;
  
  for ( j=0; j<m_StraCorn.m_ny; j++ ) { // y direction left to right
				
    gridThiis = &m_StraCorn.m_grids[j]; 
    conc_this = gridThiis->m_concChem;
    conc_other = m_gridVehicle.m_concChem;
      
    flux = gridThiis->compFlux( &m_gridVehicle, conc_this, conc_other, 
				gridThiis->m_dx/2, m_gridVehicle.m_dx/2, &deriv_this, &deriv_other);
    area = gridThiis->m_dy*gridThiis->m_dz;
    
    total_area += area;
    mass_transfer_rate += flux*area;
    //printf("j %d, area %e, total area %e\n", j, area, total_area);
    
  }

  // printf("\t total area 2sc %e\n", total_area);

  flux = mass_transfer_rate / total_area;
  return flux;
}

// Compute flux from stratum corneum to the layer down (e.g. viable epidermis)
double Skin::compFlux_sc2down()
{
  int j, idx;
  double flux, conc_this, conc_other, deriv_this, deriv_other, mass_transfer_rate, area, total_area;
  Grid *gridThiis, *gridOther;

  idx = (m_StraCorn.m_nx-1)*m_StraCorn.m_ny;
  mass_transfer_rate = total_area = 0;

  //  printf("\n compflux_sc2ve\n");
  assert(m_ViaEpd.m_ny==1);
  gridOther = m_ViaEpd.m_grids;
  conc_other = gridOther->m_concChem;

  for ( j=0; j<m_StraCorn.m_ny; j++ ) { // y direction left to right
    gridThiis = &m_StraCorn.m_grids[idx+j]; 
    conc_this = gridThiis->m_concChem;   

    flux = gridThiis->compFlux( gridOther, conc_this, conc_other, 
				gridThiis->m_dx/2, gridOther->m_dx/2, &deriv_this, &deriv_other);
    area = gridThiis->m_dy*gridThiis->m_dz;

    total_area += area;
    mass_transfer_rate += flux*area;
  }

  //  printf("\t total area sc2down %e\n", total_area);

  flux = -mass_transfer_rate / total_area;
  return flux;
}
/*
double Skin::compFlux_sc2down()
{
  int j, idx;
  double flux, area, total_area;
  Grid *gridThiis = NULL;

  idx = (m_StraCorn.m_nx-1)*m_StraCorn.m_ny;
  total_area = 0;

  //  printf("\n compflux_sc2ve\n");

  for ( j=0; j<m_StraCorn.m_ny; j++ ) { // y direction left to right
    gridThiis = &m_StraCorn.m_grids[idx+j];
    area = gridThiis->m_dy*gridThiis->m_dz;
    total_area += area;
  }
  flux = m_StraCorn.m_mass_out / total_area;
  return flux;
  }*/

// Compute flux from viable epidermis to the layer down (e.g. dermis)
double Skin::compFlux_ve2down()
{
  int j, idx;
  double flux, conc_this, conc_other, deriv_this, deriv_other, mass_transfer_rate, area, total_area;
  Grid *gridThiis, *gridOther;

  idx = (m_ViaEpd.m_nx-1)*m_ViaEpd.m_ny;
  mass_transfer_rate = total_area = 0;

  //  printf("\n compflux_ve2sk\n");
  assert(m_Dermis.m_ny==1);
  gridOther = m_Dermis.m_grids;
  conc_other = gridOther->m_concChem;

  for ( j=0; j<m_ViaEpd.m_ny; j++ ) { // y direction left to right				
    gridThiis = &m_ViaEpd.m_grids[idx+j];
    conc_this = gridThiis->m_concChem;

    flux = gridThiis->compFlux( gridOther, conc_this, conc_other, 
				gridThiis->m_dx/2, gridOther->m_dx/2, &deriv_this, &deriv_other);
    area = gridThiis->m_dy*gridThiis->m_dz;

    total_area += area;
    mass_transfer_rate += flux*area;
  }
  //  printf("\t total area ve2down %e\n", total_area);

  flux = -mass_transfer_rate / total_area;
  return flux;
}
/*double Skin::compFlux_ve2down()
{
  int j, idx;
  double flux, area, total_area;
  Grid *gridThiis = NULL;

  idx = (m_ViaEpd.m_nx-1)*m_ViaEpd.m_ny;
  total_area = 0;

  //  printf("\n compflux_ve2sk\n");

  for ( j=0; j<m_ViaEpd.m_ny; j++ ) { // y direction left to right
				
    gridThiis = &m_ViaEpd.m_grids[idx+j];
    area = gridThiis->m_dy*gridThiis->m_dz;
    total_area += area;
  }
  flux = m_ViaEpd.m_mass_out / total_area;
  return flux;
}
*/

// Compute flux from dermis to the layer down (e.g. deep tissue)
double Skin::compFlux_de2down()
{
  int j, idx;
  double flux, conc_this, conc_other, deriv_this, deriv_other, mass_transfer_rate, area, total_area;
  Grid *gridThiis, *gridOther;

  idx = (m_Dermis.m_nx-1)*m_Dermis.m_ny;
  mass_transfer_rate = total_area = 0;

  assert(m_Dermis.m_ny==1);
  gridOther = &m_gridSink;
  conc_other = gridOther->m_concChem;

  for ( j=0; j<m_Dermis.m_ny; j++ ) { // y direction left to right				
    gridThiis = &m_Dermis.m_grids[idx+j];
    conc_this = gridThiis->m_concChem;

    flux = gridThiis->compFlux( gridOther, conc_this, conc_other, 
				gridThiis->m_dx/2, gridOther->m_dx/2, &deriv_this, &deriv_other);
    area = gridThiis->m_dy*gridThiis->m_dz;

    total_area += area;
    mass_transfer_rate += flux*area;  
  }

  // printf("\t total area de2down %e\n", total_area);

  flux = -mass_transfer_rate / total_area;
  return flux;
}
/*double Skin::compFlux_de2down()
{
  int j, idx;
  double flux, area, total_area;
  Grid *gridThiis = NULL;

  idx = (m_Dermis.m_nx-1)*m_Dermis.m_ny;
  total_area = 0;

  for ( j=0; j<m_Dermis.m_ny; j++ ) { // y direction left to right
				
    gridThiis = &m_Dermis.m_grids[idx+j];
    area = gridThiis->m_dy*gridThiis->m_dz;
    total_area += area;
  }
  flux = m_Dermis.m_mass_out / total_area;
  return flux;
}*/

	
/*  ++++++++++++++++++++++++++++++++++
	I/O functions
	++++++++++++++++++++++++++++++++++ */

void Skin::getGridsConc(double *fGridsConc, int dim)
{
  assert( fGridsConc && dim==m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de+m_dim_bd );

  // extracting concentration values from the various compartments
  fGridsConc[0] = m_gridVehicle.getConcChem();
  m_StraCorn.getGridsConc( fGridsConc+m_dim_vh, m_dim_sc );
  m_ViaEpd.getGridsConc( fGridsConc+m_dim_vh+m_dim_sc, m_dim_ve );
  m_Dermis.getGridsConc( fGridsConc+m_dim_vh+m_dim_sc+m_dim_ve, m_dim_de );
  fGridsConc[m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de] = m_Blood.getConcChem();

}

void Skin::get1DConcSC(double *ret, int dim_ret)
{
  assert( dim_ret == m_StraCorn.m_n_layer_x );
  m_StraCorn.comp1DConc();
  for ( int i = 0; i < dim_ret; i ++ )
    ret[i] = m_StraCorn.m_conc1D[i];
}


void Skin::get1DCoordSC(double *ret, int dim_ret)
{
  assert( dim_ret == m_StraCorn.m_n_layer_x );
  for ( int i = 0; i < dim_ret; i ++ )
    ret[i] = m_StraCorn.m_coord1D[i];
}


void Skin::displayGrids()
{
  m_StraCorn.displayGrids();
  m_ViaEpd.displayGrids();
  m_Dermis.displayGrids();
  m_Blood.displayGrids();

  fflush(stdout);
}


void Skin::saveGrids(bool b_1st_time, const char fn[])
{
  char fn_tmp[1024];

  strcpy(fn_tmp, fn); strcat(fn_tmp, "_vh.txt");
  saveVehicle(b_1st_time, fn_tmp);

  strcpy(fn_tmp, fn); strcat(fn_tmp, "_sc.txt");
  m_StraCorn.saveGrids(b_1st_time, fn_tmp);

  strcpy(fn_tmp, fn); strcat(fn_tmp, "_ve.txt");
  m_ViaEpd.saveGrids(b_1st_time, fn_tmp);

  strcpy(fn_tmp, fn); strcat(fn_tmp, "_de.txt");
  m_Dermis.saveGrids(b_1st_time, fn_tmp);

  strcpy(fn_tmp, fn); strcat(fn_tmp, "_bd.txt");
  m_Blood.saveConc(b_1st_time, fn_tmp);
}

void Skin::saveVehicle(bool b_1st_time, const char fn[])
{
  FILE *file = NULL;
  if ( b_1st_time )
    file = fopen(fn, "w");
  else 
    file = fopen(fn, "a");
	
  fprintf(file, "%.5e\n", m_gridVehicle.getConcChem());
  fclose(file);
}

void Skin::getXCoord(double *coord_x, int dim)
{
  assert( coord_x && dim==m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de+m_dim_bd );

  // extracting X coordinates from the various compartments

  coord_x[0] = -1e10; // vehicle
  m_StraCorn.getXCoord( coord_x+m_dim_vh, m_dim_sc );
  m_ViaEpd.getXCoord( coord_x+m_dim_vh+m_dim_sc, m_dim_ve );
  m_Dermis.getXCoord( coord_x+m_dim_vh+m_dim_sc+m_dim_ve, m_dim_de );
  coord_x[m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de] = 1e10; // blood
}

void Skin::getYCoord(double *coord_y, int dim)
{
  assert( coord_y && dim==m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de+m_dim_bd );

  // extracting Y coordinates from the various compartments

  coord_y[0] = 0; // vehicle
  m_StraCorn.getYCoord( coord_y+m_dim_vh, m_dim_sc );
  m_ViaEpd.getYCoord( coord_y+m_dim_vh+m_dim_sc, m_dim_ve );
  m_Dermis.getYCoord( coord_y+m_dim_vh+m_dim_sc+m_dim_ve, m_dim_de );
  coord_y[m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de] = 0; // blood
}

void Skin::saveCoord(const char fn_x[], const char fn_y[])
{
  m_StraCorn.saveCoord(fn_x, fn_y);
  m_ViaEpd.saveCoord(fn_x, fn_y); 
  m_Dermis.saveCoord(fn_x, fn_y); 
}
/*  END <I/O functions>
	------------------------------ */
