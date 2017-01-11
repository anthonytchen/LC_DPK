#include "stdafx.h"
#include "Skin.h"

using namespace std;

/* Initialisation for possibly multiple chemicals */
void Skin::Init()
{
  m_dz_dtheta = 0.01; // fixing dz, the dimension perpendicular to x-y domain

  m_amount = NULL;
  
  m_Vehicle = NULL;
  m_SurSebum = NULL;
  m_Sebum = NULL;
  m_StraCorn = NULL;
  m_ViaEpd = NULL;
  m_Dermis = NULL;
  m_Blood = NULL;

  m_nVehicle = m_nSurSebum = m_nSebum = m_nStraCorn = m_nViaEpd = m_nDermis = 0;
  m_b_has_blood = false;
  m_nxComp = m_nyComp = 0;
  m_CompIdx = NULL;

  m_coord_sys = Cartesian;

  // If InitReaction() is not called, set m_React.idx_substrate to -1 to indicate no reaction
  m_React.idx_substrate = -1;
}

void Skin::InitReaction(int idx_substrate, int idx_product, double Vmax, double Km)
{
  /* Now only one reaction is implemented. 
     m_React needs to be expanded into an array to consider multiple reactions */

  m_React.idx_substrate = idx_substrate;
  m_React.idx_product = idx_product;
  m_React.Vmax = Vmax;
  m_React.Km = Km;
}


void Skin::Release(void)
{
  int i;

  if (!m_amount) delete [] m_amount;

  if(m_nVehicle>0){
    delete [] m_concVehicleInit;
    for (i=0; i<m_nChem*m_nVehicle; i++)
      m_Vehicle[i].Release();
    delete [] m_Vehicle;
  }

  if(m_nStraCorn>0) {
    for (i=0; i<m_nChem*m_nStraCorn; i++) 
      m_StraCorn[i].Release();
    delete [] m_StraCorn;
  }

  if(m_nViaEpd>0) {
    for (i=0; i<m_nChem*m_nViaEpd; i++)
      m_ViaEpd[i].Release();
    delete [] m_ViaEpd;
  }

  if(m_nDermis>0) {
    for (i=0; i<m_nChem*m_nDermis; i++)
      m_Dermis[i].Release();
    delete [] m_Dermis;
  }

  if (m_nBlood>0) {
    for (i=0; i<m_nChem*m_nBlood; i++)
      m_Blood[i].Release();
    delete [] m_Blood;
  }

  releaseCompMatrix();
}

/* ++++++++++++++++++++++++++
   Functions to create individual compartments 
*/

void Skin::createCompMatrix(int nxComp, int nyComp)
{
  m_nxComp = nxComp;
  m_nyComp = nyComp;
  m_CompIdx = new CompIdx*[nxComp];
  for (int i = 0; i < nxComp; i++)
    m_CompIdx[i] = new CompIdx[nyComp];

  if ( m_nBlood > 0 ) // [total] [compartments] [blood] [blood_cleared]
    m_n_amount = m_nxComp*m_nyComp+3;
  else // [total] [compartments] [sink]
    m_n_amount = m_nxComp*m_nyComp+2;
  m_amount = new double[m_n_amount];
}

void Skin::releaseCompMatrix()
{
  for (int i = 0; i < m_nxComp; i++)
    delete [] m_CompIdx[i];
  delete [] m_CompIdx;
}

void Skin::createVH(const Chemical *chemSolute, const double *conc_vehicle, const double *partition_vehicle, const double *diffu_vehicle, 
		    double coord_x_start, double coord_y_start, double xlen, double ylen, double area_vehicle, 
		    bool bInfSrc, BdyCondStr bdys,
		    double *coord_x_end, double *coord_y_end)
{
  int i, nx, ny;

  nx = ny = 1; // vehicle as a single grid

  m_concVehicleInit = new double[m_nChem];
  m_Vehicle = new Vehicle[m_nChem];

  m_Vehicle_area = area_vehicle;
  m_bInfSrc = bInfSrc; // whether the vehicle is a infinite source

  for (i=0; i<m_nChem; i++) {
    m_concVehicleInit[i] = conc_vehicle[i];
    m_Vehicle[i].Init(xlen, ylen, m_dz_dtheta, nx, ny, conc_vehicle[i], partition_vehicle[i], diffu_vehicle[i],
		      m_coord_sys, bdys.up, bdys.left, bdys.right, bdys.down);
    m_Vehicle[i].createGrids(chemSolute[i], coord_x_start, coord_y_start);
  }
  m_dim_vh = m_Vehicle[0].m_nx * m_Vehicle[0].m_ny;
  
  *coord_x_end = coord_x_start + xlen;
  *coord_y_end = coord_y_start + ylen;

}

/* surface sebum */
void Skin::createSurSB(const Chemical *chemSolute,  
		       double coord_x_start, double coord_y_start, double xlen, double ylen, int n_grids_x, int n_grids_y,
		       BdyCondStr bdys,
		       double *coord_x_end, double *coord_y_end, int idx_sursb,
		       Crystal crystal, double init_mass_solid, double k_disv_per_area, double k_rect, double Csat)
{
  int i, idx;

  assert (m_SurSebum);

  // should be done outside this function
  //  m_Sebum_Sur = new SurSebum[m_nChem];


  for (i=0; i<m_nChem; i++) {
    idx = i + idx_sursb*m_nChem;
    m_SurSebum[idx].Init(xlen, ylen, m_dz_dtheta, n_grids_x, n_grids_y,
			 m_coord_sys, bdys.up, bdys.left, bdys.right, bdys.down, 
			 crystal, init_mass_solid, k_disv_per_area, k_rect, Csat);
    m_SurSebum[idx].createGrids(chemSolute[i], coord_x_start, coord_y_start);
  }

  m_dim_sb_sur = n_grids_x*n_grids_y;

  *coord_x_end = coord_x_start + xlen;
  *coord_y_end = coord_y_start + ylen;
}


void Skin::createSB(const Chemical *chemSolute,  
		    double coord_x_start, double coord_y_start, double xlen, double ylen, int n_grids_x, int n_grids_y,
		    BdyCondStr bdys,
		    double *coord_x_end, double *coord_y_end, int idx_harsb)
{
  int i, idx;

  assert (m_Sebum);

  // should be done outside the function
  //  m_Sebum_Hur = new SurSebum[m_nChem];

  for (i=0; i<m_nChem; i++) {
    idx = i + idx_harsb*m_nChem;
    m_Sebum[idx].Init(xlen, ylen, m_dz_dtheta, n_grids_x, n_grids_y,
		      m_coord_sys, bdys.up, bdys.left, bdys.right, bdys.down);
    m_Sebum[idx].createGrids(chemSolute[i], coord_x_start, coord_y_start);
  }

  m_dim_sb_har = n_grids_x*n_grids_y;

  *coord_x_end = coord_x_start + xlen;
  *coord_y_end = coord_y_start + ylen;

}

void Skin::createSC(const Chemical *chemSolute,
		    double coord_x_start, double coord_y_start, int n_layer_x_sc, int n_layer_y_sc, double offset_y_sc,
		    BdyCondStr bdys,
		    double *coord_x_end, double *coord_y_end)
{
  int i;
  double xlen, ylen;

  m_StraCorn = new StraCorn[m_nChem];

  double g, d, s, t, water_frac_surface;
  g=.075e-6; d=40e-6; s=0.075e-6; t=0.8e-6;
  water_frac_surface = 0.55; // mass fraction of water on the surface of stratum corneum;
                             // 55% -- fully hydrated; 20% -- dry

  for (i=0; i<m_nChem; i++) {
    m_StraCorn[i].Init(g, d, s, t, m_dz_dtheta, n_layer_x_sc, n_layer_y_sc, offset_y_sc,
		       m_coord_sys, bdys.up, bdys.left, bdys.right, bdys.down); // bdy conditions: u/l/r/d
    m_StraCorn[i].createGrids(chemSolute[i], water_frac_surface, coord_x_start, coord_y_start);     
  }

  m_dim_sc = m_StraCorn[0].m_nx * m_StraCorn[0].m_ny;
  xlen = n_layer_x_sc*(g+t) + g;
  ylen = n_layer_y_sc*(d+s);

  *coord_x_end = coord_x_start + xlen;
  *coord_y_end = coord_y_start + ylen;
}

void Skin::createVE(const Chemical *chemSolute,  
		    double coord_x_start, double coord_y_start, double xlen, double ylen, int n_grids_x, int n_grids_y,
		    BdyCondStr bdys,
		    double *coord_x_end, double *coord_y_end)
{

  int i;

  m_ViaEpd = new ViaEpd[m_nChem];

  for (i=0; i<m_nChem; i++) {
    m_ViaEpd[i].Init(xlen, ylen, m_dz_dtheta, n_grids_x, n_grids_y,
		     m_coord_sys, bdys.up, bdys.left, bdys.right, bdys.down);
    m_ViaEpd[i].createGrids(chemSolute[i], coord_x_start, coord_y_start);
  }

  m_dim_ve = m_ViaEpd[0].m_nx * m_ViaEpd[0].m_ny;

  *coord_x_end = coord_x_start + xlen;
  *coord_y_end = coord_y_start + ylen;
}

void Skin::createDE(const Chemical *chemSolute,  
		    double coord_x_start, double coord_y_start, double xlen, double ylen, int n_grids_x, int n_grids_y,
		    bool b_has_blood, BdyCondStr bdys,
		    double *coord_x_end, double *coord_y_end)
{
  int i;

  m_Dermis = new Dermis[m_nChem];

  for (i=0; i<m_nChem; i++) {
    m_Dermis[i].Init(xlen, ylen, m_dz_dtheta, n_grids_x, n_grids_y, b_has_blood,
		     m_coord_sys, bdys.up, bdys.left, bdys.right, bdys.down);
    m_Dermis[i].createGrids(chemSolute[i], coord_x_start, coord_y_start);
  }

  m_dim_de = m_Dermis[0].m_nx * m_Dermis[0].m_ny;

  *coord_x_end = coord_x_start + xlen;
  *coord_y_end = coord_y_start + ylen;
}

void Skin::createBD(const double *par_dermis2blood, const double *blood_k_clear)
{
  int i;

  assert(m_b_has_blood);
  m_Blood = new Blood[m_nChem];

  for (i=0; i<m_nChem; i++) {
    m_Blood[i].Init(m_Dermis[i].m_grids->m_chemical.m_frac_unbound, blood_k_clear[i], 70, 'M');
    m_Dermis[i].InitDermisBlood(m_Blood[i].m_flow_capil, m_Blood[i].m_f_unbound, par_dermis2blood[i]);
  }

  m_dim_bd = 2;
}

int Skin::getSizeBdyRight(int chemIdx, int cIdx_i, int cIdx_j)
{
  int i, idx, size;
  Comp *pCompThis, *pCompBdy;
  
  pCompThis = m_CompIdx[cIdx_i][cIdx_j].pComp[chemIdx];

  if (cIdx_j == m_nyComp-1){ // rightmost compartment, no right boundary
    assert( pCompThis->m_BdyCond_right!=FromOther );
    size = 0;
  }
  else {
    pCompBdy = m_CompIdx[cIdx_i][cIdx_j+1].pComp[chemIdx];
    size = pCompBdy->m_nx;
  }

  return size;
}

Comp* Skin::getConcBdyRight(const double y[], int idx_up2now, int chemIdx, int cIdx_i, int cIdx_j, double concRight[])
{
  int i, idx;
  Comp *pCompThis, *pCompBdy;

  assert(concRight!=NULL);
  
  pCompThis = m_CompIdx[cIdx_i][cIdx_j].pComp[chemIdx];
  pCompBdy = m_CompIdx[cIdx_i][cIdx_j+1].pComp[chemIdx];

  idx = idx_up2now + pCompThis->m_dim;
  for (i=0; i<pCompBdy->m_nx; i++ )
    concRight[i] = y[idx+i*pCompBdy->m_ny];

  return pCompBdy;
}

int Skin::getSizeBdyDown(int chemIdx, int cIdx_i, int cIdx_j)
{
  int i, idx, size;
  Comp *pCompThis, *pCompBdy;
  
  pCompThis = m_CompIdx[cIdx_i][cIdx_j].pComp[chemIdx];

  if (cIdx_i == m_nxComp-1){ // downmost compartment, no down boundary
    assert( pCompThis->m_BdyCond_down!=FromOther );
    size = 0;
  }
  else {
    pCompBdy = m_CompIdx[cIdx_i+1][cIdx_j].pComp[chemIdx];
    size = pCompBdy->m_ny;
  }

  return size;
}

Comp* Skin::getConcBdyDown(const double y[], int idx_up2now, int chemIdx, int cIdx_i, int cIdx_j, double* concDown)
{
  int i, j, idx;
  Comp *pCompThis, *pCompBdy, *p;
  
  pCompThis = m_CompIdx[cIdx_i][cIdx_j].pComp[chemIdx];

  if (cIdx_i == m_nxComp-1){ // downmost compartment, no down boundary
    assert( pCompThis->m_BdyCond_down!=FromOther );
    concDown = NULL;
    return NULL;
  }

  pCompBdy = m_CompIdx[cIdx_i+1][cIdx_j].pComp[chemIdx];
  // concDown = new double[pCompBdy->m_ny];

  // work out the index for the down boundary
  idx = idx_up2now;
  i = cIdx_i; 
  for (j=cIdx_j; j<m_nyComp; j++) {
      p = m_CompIdx[i][j].pComp[chemIdx];
      idx += p->m_dim;
  }
  i = cIdx_i+1;
  for (j=0; j<cIdx_j; j++) {
      p = m_CompIdx[i][j].pComp[chemIdx];
      idx += p->m_dim;
  }

  // fill in concDown from y[]
  for (j=0; j<pCompBdy->m_ny; j++ )
    concDown[j] = y[idx+j];

  return pCompBdy;

}

// ----------------------------------------------------

int Skin::static_cvODE (double t, N_Vector y, N_Vector dydt, void *paras)
{
  double *p_y, *p_dydt;
	
  p_y = NV_DATA_S(y);
  p_dydt = NV_DATA_S(dydt);
  return ((Skin*)paras)->compODE_dydt(t, p_y, p_dydt);
}


int Skin::compODE_dydt (double t, const double y[], double f[])
{
  int i, j, k, idx, dim;
  double *concBdy, *concBdyRight, *concBdyDown;
  concBdy = concBdyRight = concBdyDown = NULL;

  Comp *pComp, *pCompBdyRight, *pCompBdyDown;
  pComp = pCompBdyRight = pCompBdyDown = NULL;

  CompType cType;

  /*
  if (t> 3e-12) {
    printf("\n");
  }
  */
  for (k=0; k<m_nChem; k++) { // for each chemical

    idx = k*m_dim_all;

    // for each of the compartments in m_CompIdx

    //  1. reset mass transfer between compartments
    for (i=0; i<m_nxComp; i++) {
      for (j=0; j<m_nyComp; j++) {
	pComp = m_CompIdx[i][j].pComp[k];
	pComp->setBdyMassInOutZero();
      }
    }

    //  2. calculate right hand side of differential equations
    for (i=0; i<m_nxComp; i++) {
      for (j=0; j<m_nyComp; j++) {

	pComp = m_CompIdx[i][j].pComp[k];

	// update the concentration in boundary grids according to y[]

	pCompBdyRight = pCompBdyDown = NULL;
	concBdyRight = concBdyDown = NULL;

	dim = getSizeBdyRight(k, i, j);
	if (dim) {
	  concBdyRight = new double[dim];
	  pCompBdyRight = getConcBdyRight(y, idx, k, i, j, concBdyRight);
	}

	dim = getSizeBdyDown(k, i, j);
	if (dim) {
	  concBdyDown = new double[dim];
	  pCompBdyDown = getConcBdyDown(y, idx, k, i, j, concBdyDown);
	}

	pComp->setBoundaryConc(concBdyRight, concBdyDown);
	if (concBdyRight != NULL) delete [] concBdyRight;
	if (concBdyDown != NULL) delete [] concBdyDown;

	cType = m_CompIdx[i][j].type;

	switch (cType) {

	case emVH :
	  ((Vehicle *)pComp)->compODE_dydt(t, y+idx, f+idx); // calculate dydt
	  if (m_bInfSrc)                              // infinite source, concentration doesn't change, but above dydt calculation is still needed
	    memset(f+idx, 0, sizeof(double)*((Vehicle *)pComp)->m_dim);  //  since calling compODE_dydt will calculate the flux across boundaries properly
	  break;

	case emSB :
	  ((Sebum *)pComp)->compODE_dydt(t, y+idx, f+idx);
	  //((Sebum *)pComp)->displayGridsConc(f+idx);
	  //printf("\n");
	  break;

	case emSurSB :
	  ((SurSebum *)pComp)->compODE_dydt(t, y+idx, f+idx);
	  //printf("y, t = %e, \n", t);
	  //((SurSebum *)pComp)->displayGridsConc(y+idx);
	  //printf("f, t = %e, \n", t);
	  //((SurSebum *)pComp)->displayGridsConc(f+idx);
	  //printf("\n\n");
	  break;

	case emSC :
	  ((StraCorn *)pComp)->compODE_dydt(t, y+idx, f+idx);
	  //((StraCorn *)pComp)->displayGridsConc(f+idx);
	  //printf("\n");
	  break;

	case emVE :
	  ((ViaEpd *)pComp)->compODE_dydt(t, y+idx, f+idx);
	  break;

	case emDE :
	  if (m_b_has_blood) // y[ k*m_dim_all+m_dim_all-1 ] contains the blood concentration, i.e. the last term in the differential equations
	    ((Dermis *)pComp)->updateBlood( y[ k*m_dim_all+m_dim_all-1 ] );
	  ((Dermis *)pComp)->compODE_dydt(t, y+idx, f+idx);
	  break;

	default :
	  SayBye("Compartment type not implemented");
	  break;
	}	

	pComp->passBdyMassOut(pCompBdyRight, pCompBdyDown);
	//pComp->passBdyMassOut(NULL, pCompBdyDown);
	idx += pComp->m_dim;

      } // for j
    } // for i
    
    // compute for blood 
    //exit(0);

    if (m_b_has_blood){
      // simulation is for a small skin area, but needs to multiple
      //  the mass transport due to blood flow by the actual topical application area
      double factor = m_Vehicle_area / m_Dermis[k].compTotalArea(0);

      // todo: when more than one dermis compartments are involved, need to collate mass in-out of all dermis compartments

      m_Blood[k].updateMassInOutDermis(m_Dermis[k].m_mass_into_dermis, m_Dermis[k].m_mass_outof_dermis, factor);
      m_Blood[k].compODE_dydt(t, y+idx, f+idx);
    }
    
  } // for k, each chemical species


  /* Compute for reaction 
   todo: the following has not been fully tested and may contain bugs */

  int idx_s, idx_p;
  double C_s, C_p, mass_reacted;

  if (m_React.idx_substrate != -1) { // if idx_substrate not -1, reaction needs to be considered

    idx_s = m_React.idx_substrate * m_dim_all + m_dim_vh + m_dim_sc;
    idx_p = m_React.idx_product * m_dim_all + m_dim_vh + m_dim_sc;

    // For each grid in viable epidermis
    for ( i = 0; i < m_ViaEpd[0].m_nx; i++ ){ // verticle direction up to down
      for ( j = 0; j < m_ViaEpd[0].m_ny; j++ ){ // lateral direction left to right

	idx = i*m_ViaEpd[0].m_ny + j;
	C_s = y[idx_s+idx];
	C_p = y[idx_p+idx];

	mass_reacted = m_React.Vmax * C_s / (m_React.Km + C_s);
	f[idx_s+idx] -= mass_reacted;
	f[idx_p+idx] += mass_reacted;
      }
    }


    idx_s += m_dim_ve;
    idx_p += m_dim_ve;

    // For each grid in dermis
    for ( i = 0; i < m_Dermis[0].m_nx; i++ ){ // verticle direction up to down
      for ( j = 0; j < m_Dermis[0].m_ny; j++ ){ // lateral direction left to right

	idx = i*m_Dermis[0].m_ny + j;
	C_s = y[idx_s+idx];
	C_p = y[idx_p+idx];

	mass_reacted = m_React.Vmax * C_s / (m_React.Km + C_s);
	f[idx_s+idx] -= mass_reacted;
	f[idx_p+idx] += mass_reacted;
      }
    }

  }

  return 0;
}

/* add reaction with diffusion */
void Skin::compReaction()
{
  SayBye("not implemented yet");
}


/*!
  Diffusion using method of lines (MoL) to discretise the
  PDE into ODEs
*/
void Skin::diffuseMoL(double t_start, double t_end)
{		
  int i, j, k, idx;
  double reltol, abstol, t;
  double *y = NULL;
  double *pNVs = NULL;
  N_Vector y0;
  void *cvode_mem = NULL;

  Comp *pComp = NULL;

  /* get current concentration, and set as initial conditions */

  y = new double[m_dim_all*m_nChem];

  for (k=0; k<m_nChem; k++) {
  
    idx = k*m_dim_all;

    for (i=0; i<m_nxComp; i++) {
      for (j=0; j<m_nyComp; j++) {
	pComp = m_CompIdx[i][j].pComp[k];
	pComp->getGridsConc(y+idx, pComp->m_dim);
	idx += pComp->m_dim;
      } // for j
    } // for i

    if (m_b_has_blood)
      m_Blood[k].getGridsConc(y+idx, m_Blood[k].m_dim);

  } // for k, each chemical

  /* ------------- */

  y0 = N_VMake_Serial(m_dim_all*m_nChem, y);
	
  // Call CVodeCreate to create the solver memory and specify the 
  //	 Backward Differentiation Formula and the use of a Newton iteration
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
		
  // Call CVodeInit to initialize the integrator memory and specify the
  //	right hand side function in y'=f(t,y), the inital time t_start, and the initial condition y0.
  CVodeInit(cvode_mem, static_cvODE, t_start, y0);

  // Call CVodeSStolerances to specify the scalar relative tolerance and scalar absolute tolerance.
  reltol=1e-6; abstol=1e-15;
  CVodeSStolerances(cvode_mem, reltol, abstol);
	
  // Set the pointer to user-defined data
  CVodeSetUserData(cvode_mem, this);
  CVodeSetMaxNumSteps(cvode_mem, 10000*5);
 
  CVDense(cvode_mem, m_dim_all*m_nChem);
  // CVDlsSetDenseJacFn(cvode_mem, static_cvJacobian); 
  
  CVode(cvode_mem, t_end, y0, &t, CV_NORMAL);	
	
  /* Extract calculated values */
		  
  for (k=0; k<m_nChem; k++) {

    idx = k*m_dim_all;

    for (i=0; i<m_nxComp; i++) {
      for (j=0; j<m_nyComp; j++) {
	pNVs = NV_DATA_S(y0);
	pComp = m_CompIdx[i][j].pComp[k];
	pComp->setGridsConc(pNVs+idx, pComp->m_dim);
	idx += pComp->m_dim;
      } // for j
    } // for i

    if (m_b_has_blood)
      m_Blood[k].setGridsConc(pNVs+idx, m_Blood[k].m_dim);

  } // for k, each chemical

  /* ------------------------- */
  

  N_VDestroy_Serial(y0);  
  CVodeFree(&cvode_mem);
  delete [] y;
}

 // reset vehicle concentration, partition coefficient, diffusivity
void Skin::resetVehicle(double concChem[], double partition_coef[], double diffu_coef[])
{
  int i, j, k, idx;
  for (k=0; k<m_nChem; k++) {
    
    for (i=0; i<m_Vehicle[k].m_nx; i++) {
      for (j=0; j<m_Vehicle[k].m_ny; j++) {
	
	idx = i*m_Vehicle[k].m_ny + j;
	m_Vehicle[k].m_grids[idx].m_concChem = concChem[k];
	m_Vehicle[k].m_grids[idx].m_Kw = partition_coef[k];
	m_Vehicle[k].m_grids[idx].m_D = diffu_coef[k];
      }
    }

  } // for k
}

void Skin::removeVehicle()
{
  double *concChem = new double[m_nChem];
  double *partition_coef = new double[m_nChem];
  double *diffu_coef = new double[m_nChem];

  for (int i=0; i<m_nChem; i++) {
    concChem[i] = 0;
    partition_coef[i] = 1e10; // very large partition coefficient
    diffu_coef[i] = 1e-100;   //   and very small diffusivity effectively remove the vehicle
  }
  resetVehicle(concChem, partition_coef, diffu_coef);

  delete [] concChem;
  delete [] partition_coef;
  delete [] diffu_coef;
}

// Compute flux into stratum corneum
void Skin::compFlux_2sc(double *flux)
{
  int i, j, idx;
  double f, mass_transfer_rate, total_area;
  double *concBdy = NULL;

  for ( i=0; i<m_nChem; i++ ) {

    m_Vehicle[i].setBdyMassInOutZero();

    concBdy = new double[m_StraCorn[i].m_ny];
    for (j=0; j<m_StraCorn[i].m_ny; j++)
      concBdy[j] = m_StraCorn[i].m_grids[j].m_concChem;
    m_Vehicle[i].setBoundaryConc(NULL, concBdy);
    delete [] concBdy;

    idx = (m_Vehicle[i].m_nx-1)*m_Vehicle[i].m_ny;
    mass_transfer_rate = 0;

    for ( j=0; j<m_Vehicle[i].m_ny; j++ ) { // y direction left to right    
      mass_transfer_rate += - m_Vehicle[i].compMassIrregGridsDown( m_Vehicle[i].m_grids[idx+j], m_Vehicle[i].m_grids[idx+j].m_concChem );
    }

    total_area = m_Vehicle[i].compTotalArea(3);
    flux[i] = mass_transfer_rate / total_area;
  } // for i, each chemical

}

// Compute flux from stratum corneum to the layer down (e.g. viable epidermis)
void Skin::compFlux_sc2down(double *flux)
{
  int i, j, idx;
  double f, mass_transfer_rate, area, total_area, conc_this, deriv_this, deriv_other;
  Grid *gridThiis;
  double *concBdy = NULL;

  for ( i=0; i<m_nChem; i++ ) {

    if (m_nViaEpd>0) {
      m_StraCorn[i].setBdyMassInOutZero();

      concBdy = new double[m_ViaEpd[i].m_ny];
      for (j=0; j<m_ViaEpd[i].m_ny; j++)
	concBdy[j] = m_ViaEpd[i].m_grids[j].m_concChem;
      m_StraCorn[i].setBoundaryConc(NULL, concBdy);
      delete [] concBdy;
    }

    idx = (m_StraCorn[i].m_nx-1)*m_StraCorn[i].m_ny;
    mass_transfer_rate = 0;

    for ( j=0; j<m_StraCorn[i].m_ny; j++ ) { // y direction left to right

      gridThiis = &m_StraCorn[i].m_grids[idx+j];
      conc_this = gridThiis->m_concChem;      

      if (m_nViaEpd>0)
	mass_transfer_rate += - m_StraCorn[i].compMassIrregGridsDown( *gridThiis, conc_this );
      else { // sc down is sink
	f = gridThiis->compFlux( &m_StraCorn[i].m_gridSink, conc_this, 0, gridThiis->m_dx/2, 0, &deriv_this, &deriv_other);
	area = m_StraCorn[i].compInterArea(*gridThiis, 3);
	mass_transfer_rate += - f*area;
      }
    }

    total_area = m_StraCorn[i].compTotalArea(3);
    flux[i] = mass_transfer_rate / total_area;

  } // for i, each chemical

}

double Skin::getSCYlen()
{
  double g, d, s, t;
  g=.075e-6; d=40e-6; s=0.075e-6; t=0.8e-6;
  return (d+s);
}

/*!
  Compute the mass (or mol, depending on concentration unit used) of solute in 
   each compartment
*/
double Skin::compCompartAmount()
{
  assert(m_nChem == 1);

  int i, j, idx, chemIdx=0;
  double tmp, total, frac;
  Comp *pCompThis = NULL;

  total = m_concVehicleInit[chemIdx] * m_Vehicle[chemIdx].compTotalVolume();
  m_amount[0] = total;
  
  idx = 1;
  frac = 0;
  for ( i = 0; i < m_nxComp; i++ ) {
    for ( j = 0; j < m_nyComp; j++ ) {
      pCompThis = m_CompIdx[i][j].pComp[chemIdx];
      tmp = pCompThis->getAmount() / total;
      m_amount[ idx + i*m_nyComp+j ] = tmp;
      frac += tmp;
    }
  }
  
  idx = 1 + m_nxComp*m_nyComp;
  if ( m_nBlood > 0 ) { // has blood
    tmp = m_Blood[chemIdx].getAmount() / total;
    m_amount[idx] = tmp;
    frac += tmp;

    tmp = m_Blood[chemIdx].getClearedAmount() / total;
    m_amount[idx+1] = tmp;
    frac += tmp;
  }
  else{ // has sink a.k.a. receptor fluid
    m_amount[idx] = 1-frac; // this is not ideal but has been verified
                            //  when switching off sink, i.e. using zero flux at the bottom
  }

  return total;

  /*
  fLayersAmount[0] = m_concVehicleInit[i] * m_Vehicle[i].compTotalVolume();
  fLayersAmount[1] = m_Vehicle[i].getAmount();

  m_StraCorn[i].getAmount(fLayersAmount+2, fLayersAmount+3, fLayersAmount+4);
  fLayersAmount[5] = m_ViaEpd[i].getAmount();
  fLayersAmount[6] = m_Dermis[i].getAmount();
  fLayersAmount[7] = m_Blood[i].getAmount();
  fLayersAmount[8] = m_Blood[i].getClearedAmount();
  */
}


/*!
  Compute the mass (or mol, depending on concentration unit used) of solute in 
   each layer
*/
void Skin::getLayersAmount(double *fLayersAmount, int dim, int idx_chem)
{
  // order of values (9 items):
  //    [initial amount in vehicle] [current amount in vehicle]
  //    [current amount in stratum corneum] [current amount in lipid of stratum corneum] [current amount in corneocyte of stratum corneum]
  //    [current amount in viable epidermis]
  //    [current amount in dermis]
  //    [current amount in blood]
  //    [current amount in skin, i.e. cumulative amount cleared from blood]
  assert( fLayersAmount && dim==9 );

  int i = idx_chem;

  // Compute the amount in vehicle
  fLayersAmount[0] = m_concVehicleInit[i] * m_Vehicle[i].compTotalVolume();
  fLayersAmount[1] = m_Vehicle[i].getAmount();

  m_StraCorn[i].getAmount(fLayersAmount+2, fLayersAmount+3, fLayersAmount+4);
  fLayersAmount[5] = m_ViaEpd[i].getAmount();
  fLayersAmount[6] = m_Dermis[i].getAmount();
  fLayersAmount[7] = m_Blood[i].getAmount();
  fLayersAmount[8] = m_Blood[i].getClearedAmount();
}

	
/*  ++++++++++++++++++++++++++++++++++
	I/O functions
	++++++++++++++++++++++++++++++++++ */

void Skin::getGridsConc(double *fGridsConc, int dim, int idx_chem)
{
  assert( fGridsConc && dim==m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de+m_dim_bd-1 );

  int i = idx_chem;

  // extracting concentration values from the various compartments
  m_Vehicle[i].getGridsConc( fGridsConc, m_dim_vh );
  m_StraCorn[i].getGridsConc( fGridsConc+m_dim_vh, m_dim_sc );
  m_ViaEpd[i].getGridsConc( fGridsConc+m_dim_vh+m_dim_sc, m_dim_ve );
  m_Dermis[i].getGridsConc( fGridsConc+m_dim_vh+m_dim_sc+m_dim_ve, m_dim_de );
  fGridsConc[m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de] = m_Blood[i].getConcChem();

}

void Skin::get1DConcSC(double *ret, int dim_ret, int idx_chem)
{
  assert( dim_ret == m_StraCorn[idx_chem].m_n_layer_x );
  m_StraCorn[idx_chem].comp1DConc();
  for ( int i = 0; i < dim_ret; i ++ )
    ret[i] = m_StraCorn[idx_chem].m_conc1D[i];
}


void Skin::get1DCoordSC(double *ret, int dim_ret, int idx_chem)
{
  assert( dim_ret == m_StraCorn[idx_chem].m_n_layer_x );
  for ( int i = 0; i < dim_ret; i ++ )
    ret[i] = m_StraCorn[idx_chem].m_coord1D[i];
}


void Skin::displayGrids()
{
  int j, idx;
  SayBye("Todo: improve this function");

  m_Vehicle[0].displayGrids();
  if (m_nStraCorn>0)
    m_StraCorn[0].displayGrids();
  if (m_nViaEpd>0)
    m_ViaEpd[0].displayGrids();
  if (m_nDermis>0)
    m_Dermis[0].displayGrids();
  if (m_b_has_blood)
    m_Blood[0].displayGrids();

  fflush(stdout);
}


void Skin::saveGrids(bool b_1st_time, const char fn[])
{
  int i, j, idx;
  char fn_tmp[1024];

#ifdef _DEBUG_
  printf("  Mass in components:\n");
#endif
  
  for (i=0; i<m_nChem; i++) {

    for (j=0; j<m_nVehicle; j++) {
      idx = i*m_nVehicle+j;
      sprintf(fn_tmp, "%s_vh_chem%d_comp%d.txt", fn, i, j);
      m_Vehicle[idx].saveGrids(b_1st_time, fn_tmp);
#ifdef _DEBUG_
      printf( "\tVehicle #%d \t %.3e\n", idx, m_Vehicle[idx].getAmount() );
#endif
    }

    for (j=0; j<m_nStraCorn; j++) {
      idx = i*m_nStraCorn+j;
      sprintf(fn_tmp, "%s_sc_chem%d_comp%d.txt", fn, i, j);
      m_StraCorn[idx].saveGrids(b_1st_time, fn_tmp);
#ifdef _DEBUG_
      printf( "\tStraCorn #%d \t %.3e\n", idx, m_StraCorn[idx].Comp::getAmount() );
#endif
    }

    for (j=0; j<m_nSurSebum; j++) {
      idx = i*m_nSurSebum+j;
      sprintf(fn_tmp, "%s_sursb_chem%d_comp%d.txt", fn, i, j);
      m_SurSebum[idx].saveGrids(b_1st_time, fn_tmp);
#ifdef _DEBUG_
      printf( "\tSurSebum #%d \t %.3e\n", idx, m_SurSebum[idx].getAmount() );
#endif
    }

    for (j=0; j<m_nSebum; j++) {
      idx = i*m_nSebum+j;
      sprintf(fn_tmp, "%s_sb_chem%d_comp%d.txt", fn, i, j);
      m_Sebum[idx].saveGrids(b_1st_time, fn_tmp);
#ifdef _DEBUG_
      printf( "\tSebum #%d \t %.3e\n", idx, m_Sebum[idx].getAmount() );
#endif
    }

    for (j=0; j<m_nViaEpd; j++) {
      idx = i*m_nViaEpd+j;
      sprintf(fn_tmp, "%s_ve_chem%d_comp%d.txt", fn, i, j);
      m_ViaEpd[idx].saveGrids(b_1st_time, fn_tmp);
#ifdef _DEBUG_
      printf( "\tViaEpd #%d \t %.3e\n", idx, m_ViaEpd[idx].getAmount() );
#endif
    }

    for (j=0; j<m_nDermis; j++) {
      idx = i*m_nDermis+j;
      sprintf(fn_tmp, "%s_de_chem%d_comp%d.txt", fn, i, j);
      m_Dermis[idx].saveGrids(b_1st_time, fn_tmp);
#ifdef _DEBUG_
      printf( "\tDermis #%d \t %.3e\n", idx, m_Dermis[idx].getAmount() );
#endif
    }

    if (m_b_has_blood) {
      sprintf(fn_tmp, "%s_bd_chem%d.txt", fn, i);
      m_Blood[i].saveConc(b_1st_time, fn_tmp);
    }
  }
}

/*! Save the total amount and percentage in individual compartments
 */
void Skin::saveAmount(bool b_1st_time, const char fn[])
{
  int i;
  FILE *file = NULL;

  if ( b_1st_time )
    file = fopen(fn, "w");
  else 
    file = fopen(fn, "a");
	
  for ( i = 0; i < m_n_amount; i++ ){
    fprintf(file, "%.5e\t", m_amount[i]);
  }
  fprintf(file, "\n");

  fclose(file);
}


void Skin::getXCoord(double *coord_x, int dim)
{
  assert( coord_x && dim==m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de+m_dim_bd-1 );

  // extracting X coordinates from the various compartments

  coord_x[0] = -1e10; // vehicle
  m_StraCorn[0].getXCoord( coord_x+m_dim_vh, m_dim_sc );
  m_ViaEpd[0].getXCoord( coord_x+m_dim_vh+m_dim_sc, m_dim_ve );
  m_Dermis[0].getXCoord( coord_x+m_dim_vh+m_dim_sc+m_dim_ve, m_dim_de );
  coord_x[m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de] = 1e10; // blood
}

void Skin::getYCoord(double *coord_y, int dim)
{
  assert( coord_y && dim==m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de+m_dim_bd-1 );

  // extracting Y coordinates from the various compartments

  coord_y[0] = 0; // vehicle
  m_StraCorn[0].getYCoord( coord_y+m_dim_vh, m_dim_sc );
  m_ViaEpd[0].getYCoord( coord_y+m_dim_vh+m_dim_sc, m_dim_ve );
  m_Dermis[0].getYCoord( coord_y+m_dim_vh+m_dim_sc+m_dim_ve, m_dim_de );
  coord_y[m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de] = 0; // blood
}

void Skin::saveCoord(const char fn_x[], const char fn_y[])
{
  if (m_nStraCorn>0)
    m_StraCorn[0].saveCoord(fn_x, fn_y);
  if (m_nViaEpd>0)
    m_ViaEpd[0].saveCoord(fn_x, fn_y); 
  if (m_nDermis>0)
    m_Dermis[0].saveCoord(fn_x, fn_y); 
}


void Skin::setScProperties(double lip_Kw, double lip_D, double cc_Kw, double cc_D)
{
  int i;
  assert (m_nChem == 1);
  if(m_nStraCorn>0) {
    for (i=0; i<m_nChem*m_nStraCorn; i++) 
      m_StraCorn[i].setGridsProperties(lip_Kw, lip_D, cc_Kw, cc_D);
  }
}
/*  END <I/O functions>
	------------------------------ */
