#include "stdafx.h"
#include "Skin_Setup.h"

/*! setup the arrangement of compartments
 */
void Skin_Setup::InitConfig(Chemical* chemSolute, Config& conf)
{
  strcpy(m_sComps, conf.m_sComps);

  if ( !strcmp(m_sComps, "VS") ||       //vehicle, stratum corneum
       !strcmp(m_sComps, "VSV") ||      //vehicle, stratum corneum, viable epidermis
       !strcmp(m_sComps, "VSVD") ||     //vehicle, stratum cornrum, viable epidermis, dermis
       !strcmp(m_sComps, "VSVDB") ) {   //vehicle, stratum cornrum, viable epidermis, dermis, blood
    InitVecCompart(chemSolute, conf.m_nChem,
		   &conf.m_conc_vehicle, &conf.m_partition_vehicle, &conf.m_diffu_vehicle,       
		   conf.m_x_len_vehicle, conf.m_area_vehicle,
		   conf.m_n_layer_x_sc, conf.m_n_layer_y_sc, conf.m_offset_y_sc,
		   conf.m_x_len_ve, conf.m_n_grids_x_ve, conf.m_x_len_de, conf.m_n_grids_x_de,
		   &conf.m_partition_dermis2blood, &conf.m_Kclear_blood, conf.m_bInfSrc);
  }
  else {
    SayBye("Options not implemented");
  }
    
}

void Skin_Setup::Release(void)
{
  Skin::Release();
}

//
void Skin_Setup::InitVecCompart(Chemical *chemSolute, int nChem,
				double *conc_vehicle, double *partition_vehicle, double *diffu_vehicle, 
				double dx_vehicle, double area_vehicle, 
				int n_layer_x_sc, int n_layer_y_sc, double offset_y_sc,
				double x_len_ve, int n_grids_x_ve, double x_len_de, int n_grids_x_de,
				double *par_dermis2blood, double *blood_k_clear, bool bInfSrc)
{
  if ( nChem > 1 )
    SayBye("Not implemented: diffusion of multiple compounds through 1D vector setup of compartments");
  m_nChem = nChem;
  
  Skin::Init();
  
  int i, nxComp, nyComp;
  double coord_x_start, coord_y_start, coord_x_end, coord_y_end;
  double x_len_sc, y_len_sc;

  m_dz_dtheta = 0.01; // fixing dz, the dimension perpendicular to x-y domain

  /* setup compartment matrix, actually a vector here */
  
  nyComp = 1;
  nxComp = strlen(m_sComps);

  if ( nxComp < 2 )
    SayBye("Error! At least two compartments are needed");
  m_nVehicle = m_nStraCorn = 1;
  
  if ( nxComp > 4 ) { // VSVDB
    m_b_has_blood = true;
    nxComp = 4; // blood is not a 2D diffusion compartment, thus doesn't count here
  }
  createCompMatrix(nxComp, nyComp);  

  /* define possible boundary condition (up, left, right, down) */
  BdyCondStr bdys_top = {ZeroFlux,  Periodic, Periodic, FromOther}; // top compartment
  BdyCondStr bdys_mid = {FromOther, Periodic, Periodic, FromOther}; // middle compartment
  BdyCondStr bdys_bot = {FromOther, Periodic, Periodic, ZeroFlux};  // bottom compartment

  
  /*  set up the compartments */

  y_len_sc = getSCYlen() * n_layer_y_sc;

  // Vehicle
  coord_x_start = 0; coord_y_start = 0;
  createVH(chemSolute, conc_vehicle, partition_vehicle, diffu_vehicle,
	   coord_x_start, coord_y_start, dx_vehicle, y_len_sc, area_vehicle,
	   bInfSrc, bdys_top, &coord_x_end, &coord_y_end);
  m_CompIdx[0][0].type = emVH;
  m_CompIdx[0][0].pComp = new Comp*[1];
  m_CompIdx[0][0].pComp[0] = &m_Vehicle[0];


  // SC
  coord_x_start = dx_vehicle; coord_y_start = 0;
  
  if ( nxComp == 2 )
    createSC(chemSolute, coord_x_start, coord_y_start, n_layer_x_sc, n_layer_y_sc, offset_y_sc, 
	     bdys_bot, &coord_x_end, &coord_y_end);
  else
    createSC(chemSolute, coord_x_start, coord_y_start, n_layer_x_sc, n_layer_y_sc, offset_y_sc, 
	     bdys_mid, &coord_x_end, &coord_y_end);
  
  m_CompIdx[1][0].type = emSC;
  m_CompIdx[1][0].pComp = new Comp*[1];
  m_CompIdx[1][0].pComp[0] = &m_StraCorn[0];

  m_dim_all = m_dim_vh + m_dim_sc;
    
  // VE
  if ( nxComp > 2 ) {
    m_nViaEpd = 1;
    m_dim_all += m_dim_ve;
    coord_x_start = coord_x_end; coord_y_start = 0;
    
    if ( nxComp == 3 )
      createVE(chemSolute, coord_x_start, coord_y_start, x_len_ve, y_len_sc, n_grids_x_ve, 1,
	       bdys_bot, &coord_x_end, &coord_y_end);
    else
      createVE(chemSolute, coord_x_start, coord_y_start, x_len_ve, y_len_sc, n_grids_x_ve, 1,
	       bdys_mid, &coord_x_end, &coord_y_end);
    
    m_CompIdx[2][0].type = emVE;
    m_CompIdx[2][0].pComp = new Comp*[1];
    m_CompIdx[2][0].pComp[0] = &m_ViaEpd[0];
  }

  // DE
  if ( nxComp > 3 ) {
    m_nDermis = 1;
    m_dim_all += m_dim_de;
    coord_x_start = coord_x_end; coord_y_start = 0;
    createDE(chemSolute, coord_x_start, coord_y_start, x_len_de, y_len_sc, n_grids_x_de, 1, m_b_has_blood,
	     bdys_bot, &coord_x_end, &coord_y_end);
    m_CompIdx[3][0].type = emDE;
    m_CompIdx[3][0].pComp = new Comp*[1];
    m_CompIdx[3][0].pComp[0] = &m_Dermis[0];
  }

  // BD
  if ( strlen(m_sComps) == 5 ) {
    m_dim_all += m_dim_bd;
    createBD(par_dermis2blood, blood_k_clear);
  }

  
  /* link the compartments through boundary setting
     bdy conditions: up/left/right/down */
  
  for (i=0; i<m_nChem; i++) {

    // Vehicle
    m_Vehicle[i].createBoundary(0, m_StraCorn[i].m_ny);    
    m_Vehicle[i].setBoundaryGrids(NULL, m_StraCorn[i].m_grids);

    // Stratum corneum
    if ( nxComp == 2 ) {
          m_StraCorn[i].createBoundary(0, 0);
	  m_StraCorn[i].setBoundaryGrids(NULL, NULL);
    }
    else {
      m_StraCorn[i].createBoundary(0, m_ViaEpd[i].m_ny);
      m_StraCorn[i].setBoundaryGrids(NULL, m_ViaEpd[i].m_grids);
    }

    // Viable epidermis
    if ( nxComp == 3 ) {
      m_ViaEpd[i].createBoundary(0, 0);
      m_ViaEpd[i].setBoundaryGrids(NULL, NULL);
    }
    else {
      m_ViaEpd[i].createBoundary(0, m_Dermis[i].m_ny);
      m_ViaEpd[i].setBoundaryGrids(NULL, m_Dermis[i].m_grids);
    }

    // Dermis
    if ( nxComp > 3) {
      m_Dermis[i].createBoundary(0, 0);    
      m_Dermis[i].setBoundaryGrids(NULL, NULL);
    }

  }


  // If InitReaction() is not called, set m_React.idx_substrate to -1 to indicate no reaction
  m_React.idx_substrate = -1;
}

/*!
 */
void Skin_Setup::InitMtxCompart(Chemical *chemSolute, int nChem,
				double *conc_vehicle, double *partition_vehicle, double *diffu_vehicle, 
				double dx_vehicle, double area_vehicle, 
				int n_layer_x_sc, int n_layer_y_sc, double offset_y_sc,
				double x_len_ve, int n_grids_x_ve, double x_len_de, int n_grids_x_de,
				double *par_dermis2blood, double *blood_k_clear, bool bInfSrc)
{
  if ( nChem > 1 )
    SayBye("Not implemented: diffusion of multiple compounds through 2D matrix setup of compartments");
}

/*
void Skin_Setup::InitConfig(Chemical *chemSolute, Config &conf)
{
  Init(chemSolute, conf.m_nChem, &conf.m_conc_vehicle, &conf.m_partition_vehicle, &conf.m_diffu_vehicle,       
       conf.m_x_len_vehicle, conf.m_area_vehicle,
       conf.m_n_layer_x_sc, conf.m_n_layer_y_sc, conf.m_offset_y_sc,
       conf.m_x_len_ve, conf.m_n_grids_x_ve, conf.m_x_len_de, conf.m_n_grids_x_de,
       &conf.m_partition_dermis2blood, &conf.m_Kclear_blood, conf.m_bInfSrc);
}
*/



