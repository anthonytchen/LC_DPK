#include "stdafx.h"
#include "Config.h"

void Config::ReadConfigFile(const char fn[])
{
  int i;
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
		
      if ( !name.compare("CHEM_MW") ) // molecular weight
	m_chem.m_mw = value_d;      
      else if ( !name.compare("CHEM_KOW") ) // partition coefficient between octanol and water
	m_chem.m_K_ow = value_d;	
      else if ( !name.compare("CHEM_PKA") ) // pKa -- acide dissociation constant
	m_chem.m_pKa = value_d ; 
      else if ( !name.compare("CHEM_NONION") ) // fraction of solute non-ionised at pH 7.4
	m_chem.m_frac_non_ion = value_d; 
      else if ( !name.compare("CHEM_UNBND") )  // fraction of solute unbound in a 2.7% albumin solution at pH 7.4    
	m_chem.m_frac_unbound = value_d; 
      else if ( !name.compare("CHEM_ACIDBASE") )
	m_chem.m_acid_base = value_c;
      // else if ( !name.compare("CHEM_ACIDBASE") )
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
*/
void Config::InitChemical(Chemical& target)
{
  target.Init(m_chem.m_mw, m_chem.m_K_ow, m_chem.m_pKa, m_chem.m_frac_non_ion, m_chem.m_frac_unbound, m_chem.m_acid_base);
}

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
