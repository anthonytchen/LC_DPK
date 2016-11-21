#include "stdafx.h"
#include "Config.h"

void Config::ReadConfigFile(const char fn[])
{
  int i, value_i;
  bool blank_line;
  string line, name, value, msg;
  char value_c;
  double value_d;
  vector<string> tokens;
  
  ifstream file (fn, ifstream::in);
  if (file.is_open()) {
    while ( getline(file, line) ) {

      i = line.find_last_not_of(" \n\r\t");
      if ( i != line.npos )
	line.erase(i+1); // trim last "white" spaces
      
      if ( line.at(0) == '#' ) // this line is comment
	continue;
      
      blank_line = true;      
      for ( i=0; i<line.size(); i++ ) {
	if ( !isspace(line.at(i)) ) {
	  blank_line = false;
	  break;
	}
      }
      if ( blank_line ) // this line is blank
	continue;

      tokens.clear();
      Tokenize(line, tokens);
      if ( tokens.size() != 2 )
	continue;
      name = tokens.at(0);     
      value = tokens.at(1);
      value_c = value.at( value.find_first_not_of(" \n\r\t") );
      value_d = strtod( value.c_str(), NULL );
      value_i = strtol( value.c_str(), NULL, 10 );

      // parameters relating to chemical(s)
      if ( !name.compare("CHEM_NO") ) // number of compounds
	m_nChem = value_i;
      else if ( !name.compare("CHEM_MW") ) // molecular weight
	m_mw = value_d;      
      else if ( !name.compare("CHEM_KOW") ) // partition coefficient between octanol and water
	m_K_ow = value_d;	
      else if ( !name.compare("CHEM_PKA") ) // pKa -- acide dissociation constant
	m_pKa = value_d ; 
      else if ( !name.compare("CHEM_NONION") ) // fraction of solute non-ionised at pH 7.4
	m_frac_non_ion = value_d; 
      else if ( !name.compare("CHEM_UNBND") )  // fraction of solute unbound in a 2.7% albumin solution at pH 7.4    
	m_frac_unbound = value_d; 
      else if ( !name.compare("CHEM_ACIDBASE") )
	m_acid_base = value_c;
      else if ( !name.compare("CHEM_PAR_VEH") )
	m_partition_vehicle = value_d;
      else if ( !name.compare("CHEM_DIF_VEH") )
	m_diffu_vehicle = value_d;
      
      // parameters relating to the compartments

      else if ( !name.compare("VEH_INIT_CONC") )
	m_conc_vehicle = value_d;
      else if ( !name.compare("VEH_DX") )
	m_x_len_vehicle = value_d;
      else if ( !name.compare("VEH_AREA") )
	m_area_vehicle = value_d;
      else if ( !name.compare("VEH_INFINITE") )
	m_bInfSrc = value_i;
      
      else if ( !name.compare("SKIN_N_LAYER_X_SC") )
	m_n_layer_x_sc = value_i;
      else if ( !name.compare("SKIN_N_LAYER_Y_SC") )
	m_n_layer_y_sc = value_i;
      else if ( !name.compare("SKIN_OFFSET_Y_SC") )
	m_offset_y_sc = value_d;

      else if ( !name.compare("SKIN_N_GRIDS_X_VE") )
	m_n_grids_x_ve = value_i;
      else if ( !name.compare("SKIN_N_GRIDS_X_DE") )
	m_n_grids_x_de = value_i;
      else if ( !name.compare("SKIN_LEN_X_VE") )
	m_x_len_ve = value_d;
      else if ( !name.compare("SKIN_LEN_X_DE") )
	m_x_len_de = value_d;
      else if ( !name.compare("SKIN_PARTITION_DE2BD") ) // dermis to blood partition
	m_partition_dermis2blood = value_d;
      else if ( !name.compare("SKIN_CLEAR_BD") ) // blood clearance rate
	m_Kclear_blood = value_d;

      else if ( !name.compare("SKIN_N_GRIDS_X_SB_SUR") )
	m_n_grids_x_sb_sur = value_i;
      else if ( !name.compare("SKIN_LEN_X_SB_SUR") )
	m_x_len_sb_sur = value_d;
      else if ( !name.compare("SKIN_N_GRIDS_Y_SB_SUR") )
	m_n_grids_y_sb_sur = value_i;
      
      else if ( !name.compare("SKIN_N_GRIDS_X_SB_HAR") )
	m_n_grids_x_sb_har = value_i;      
      else if ( !name.compare("SKIN_N_GRIDS_Y_SB_HAR") )
	m_n_grids_y_sb_har = value_i;
      else if ( !name.compare("SKIN_LEN_Y_SB_HAR") )
	m_y_len_sb_har = value_d;


      // name not found
      else {
	msg = "Config name <" + name + "> not known";
	SayBye(msg.c_str());
      }
      
		
    }
    file.close();
  }
  else {
    msg = "Config file <" + string(fn) + "> doesn't exist";
    SayBye(msg.c_str());
  }
}

/*!
  Initialise objective <Chemical>

void Config::InitChemical(Chemical& target)
{
  target.Init(m_chem.m_mw, m_chem.m_K_ow, m_chem.m_pKa, m_chem.m_frac_non_ion, m_chem.m_frac_unbound, m_chem.m_acid_base);
}
*/

/*!
  Tokenize copied from (as at 30 June 2016)
  http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html
*/
void Config::Tokenize(const string& str, vector<string>& tokens, const string& delimiters)
{
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos     = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}
