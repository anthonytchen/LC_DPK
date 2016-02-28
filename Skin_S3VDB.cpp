#include "stdafx.h"
#include "Skin_S3VDB.h"

/*void Skin_S3VDB::Init(Chemical *chemSolute, int nChem, 
		      double *conc_vehicle, double *partition_vehicle, double *diffu_vehicle, 
		      double x_len_sb, double y_len_sb, 
		      int n_layer_x_sc, int n_layer_y_sc, double offset_y_sc,
		      double x_len_ve, int n_grids_x_ve, double x_len_de, int n_grids_x_de,
		      double *par_dermis2blood, double *blood_k_clear,
		      bool bInfSrc)*/
void Skin_S3VDB::Init(Chemical *chemSolute, int nChem, 
		      double *conc_vehicle, double *partition_vehicle, double *diffu_vehicle, 
		      double x_len_sb_sur, int n_grids_x_sb_sur, double y_len_sb_har, int n_grids_y_sb_har,
		      int n_layer_x_sc, int n_layer_y_sc, double offset_y_sc,
		      int x_len_ve, int n_grids_x_ve, double x_len_de, int n_grids_x_de,
		      double *par_dermis2blood=NULL, double *blood_k_clear=NULL)
{
  int i, nxComp, nyComp, n_grids_y_sb_sur, dim;
  double coord_x_start, coord_y_start, coord_x_end, coord_y_end;
  double x_len_sc, y_len_sc;

  m_dz_dtheta = 0.01; // fixing dz, the dimension perpendicular to x-y domain

  // setup compartment matrix
  //  [ SurSB, SurSB ]
  //  [ SC   , HarSB ]
  nxComp = 2; nyComp = 2;
  m_nStraCorn = m_nSebum = 1;
  m_nSurSebum = 2;
  m_nVehicle = m_nViaEpd = m_nDermis = m_nBlood = 0;
  createCompMatrix(nxComp, nyComp);

  n_grids_y_sb_sur = 5;

  // boundary condition: up, left, right, down
  BdyCondStr bdys_sb_sur1 = {ZeroFlux,  ZeroFlux, FromOther, FromOther};
  BdyCondStr bdys_sb_sur2 = {ZeroFlux,  FromOther, ZeroFlux, FromOther};
  BdyCondStr bdys_sb_har = {FromOther,  FromOther, ZeroFlux, ZeroFlux};
  BdyCondStr bdys_sc = {FromOther, ZeroFlux, FromOther, ZeroFlux};
  //BdyCondStr bdys_ve = bdys_sc;
  //BdyCondStr bdys_de = {FromOther, Periodic, Periodic, ZeroFlux};

  m_nChem = nChem;

  /*  set up the compartments */

  y_len_sc = getSCYlen() * n_layer_y_sc;
  dim = 0;

  m_SurSebum = new SurSebum[m_nChem*2];
  m_Sebum = new Sebum[m_nChem*1];

  // SB, surface sebum
  coord_x_start = 0; coord_y_start = 0;
  createSurSB(chemSolute, coord_x_start, coord_y_start, x_len_sb_sur, y_len_sc, n_grids_x_sb_sur, n_grids_y_sb_sur,
	      bdys_sb_sur1, &coord_x_end, &coord_y_end, 0);
  m_SurSebum[0].setGridConc(1.0, 0, 0); // give the first grid some initial concentration
  m_CompIdx[0][0].type = emSurSB;
  m_CompIdx[0][0].pComp = new Comp*[1];
  m_CompIdx[0][0].pComp[0] = &m_SurSebum[0];
  dim += m_SurSebum[0].m_dim;

  coord_x_start = 0; coord_y_start = y_len_sc;
  createSurSB(chemSolute, coord_x_start, coord_y_start, x_len_sb_sur, y_len_sb_har, n_grids_x_sb_sur, 1,
	      bdys_sb_sur2, &coord_x_end, &coord_y_end, 1);
  m_CompIdx[0][1].type = emSurSB;
  m_CompIdx[0][1].pComp = new Comp*[1];
  m_CompIdx[0][1].pComp[0] = &m_SurSebum[1];
  dim += m_SurSebum[1].m_dim;

  // SC
  coord_x_start = coord_x_end; coord_y_start = 0;
  createSC(chemSolute, coord_x_start, coord_y_start, n_layer_x_sc, n_layer_y_sc, offset_y_sc, 
	   bdys_sc, &coord_x_end, &coord_y_end);
  m_CompIdx[1][0].type = emSC;
  m_CompIdx[1][0].pComp = new Comp*[1];
  m_CompIdx[1][0].pComp[0] = &m_StraCorn[0]; 
  dim += m_StraCorn[0].m_dim;

  // SB, hair sebum
  coord_x_start = coord_x_start; coord_y_start = y_len_sc;
  createSB(chemSolute, coord_x_start, coord_y_start, x_len_sb_sur, y_len_sb_har, n_grids_x_sb_sur, 1,
	   bdys_sb_har, &coord_x_end, &coord_y_end, 0);
  m_CompIdx[1][1].type = emSB;
  m_CompIdx[1][1].pComp = new Comp*[1];
  m_CompIdx[1][1].pComp[0] = &m_Sebum[0];
  dim += m_Sebum[0].m_dim;

  /*
  // VE
  coord_x_start = coord_x_end; coord_y_start = 0;
  createVE(chemSolute, coord_x_start, coord_y_start, x_len_ve, y_len_sc, n_grids_x_ve, 1,
	   bdys_ve, &coord_x_end, &coord_y_end);

  // DE
  coord_x_start = coord_x_end; coord_y_start = 0;
  createDE(chemSolute, coord_x_start, coord_y_start, x_len_de, y_len_sc, n_grids_x_de, 1, true,
	   bdys_de, &coord_x_end, &coord_y_end);

  // BD
  createBD(par_dermis2blood, blood_k_clear);
  */

  /* link the compartments through boundary setting
     bdy conditions: up/left/right/down */

  for (i=0; i<m_nChem; i++) {

    m_SurSebum[i*2+0].createBoundary(m_SurSebum[i*2+1].m_nx, m_StraCorn[i].m_ny);
    m_SurSebum[i*2+0].setBoundaryGrids(m_SurSebum[i*2+1].m_grids, m_StraCorn[i].m_grids);

    m_SurSebum[i*2+1].createBoundary(0, m_Sebum[i].m_ny);
    m_SurSebum[i*2+1].setBoundaryGrids(NULL, m_Sebum[i].m_grids);

    m_StraCorn[i].createBoundary(m_Sebum[i].m_nx, 0);
    m_StraCorn[i].setBoundaryGrids(m_Sebum[i].m_grids, NULL);

    m_Sebum[i].createBoundary(0, 0);
    m_Sebum[i].setBoundaryGrids(NULL, NULL);

  }

  // overall dimension
  m_dim_all =  dim;

  // If InitReaction() is not called, set m_React.idx_substrate to -1 to indicate no reaction
  m_React.idx_substrate = -1;
}

void Skin_S3VDB::Release(void)
{
  Skin::Release();
}

