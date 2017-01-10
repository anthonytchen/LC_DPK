/*!
  The header file for class Config, which reads from a text configuration
  file the properties of the chemical(s), the setup of the skin compartments and their dimensions, and the simulation parameters,
  and then feeds them to the appropriate classes (e.g. Chemical, Skin etc.)
 */
#ifndef _H_CONFIG_
#define _H_CONFIG_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
//#include "Chemical.h"

using namespace std;

class Config
{
 public:
  /* configuration of compartments */
  char m_sComps[1024];
  
  /* parameters relating to chemical (s) */
  int m_nChem; // number of compounds
  double m_mw, m_K_ow, m_pKa, m_frac_non_ion, m_frac_unbound, m_acid_base,
    m_partition_vehicle, m_diffu_vehicle;
  
  /* parameters relating to compartments (vehicle, skin, etc.) */
  double m_conc_vehicle, m_x_len_vehicle, m_area_vehicle;
  int m_bInfSrc;
  
  int m_n_layer_x_sc, m_n_layer_y_sc;
  double m_offset_y_sc;

  int m_n_grids_x_ve, m_n_grids_x_de;
  double m_x_len_ve, m_x_len_de;

  double m_partition_dermis2blood, m_Kclear_blood;

  int m_n_grids_x_sb_sur, m_n_grids_y_sb_sur, m_n_grids_x_sb_har, m_n_grids_y_sb_har;
  double m_x_len_sb_sur, m_y_len_sb_har;
  
  /* parameters relating to simulation */
  
 public:
  void ReadConfigFile(const char []);
  // void InitChemical(Chemical&);
  //void InitSkin();
  
  void Tokenize(const string& str, vector<string>& tokens, const string& delimiters = " ");
};

#endif
