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
#include "Chemical.h"

using namespace std;

class Config
{
 private:
  /* parameters relating to chemical (s) */
  Chemical m_chem;
  int m_nChem; // number of compounds
  
  /* parameters relating to skin */
  /* parameters relating to simulation */
  
 public:
  void ReadConfigFile(const char []);
  void InitChemical(Chemical&);
  void InitSkin();
  
  void Tokenize(const string& str, vector<string>& tokens, const string& delimiters = " ");
};

#endif
