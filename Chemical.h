/* The header file for chemical being applied to skin
   It should include all the properties (descriptors) of the chemical */
#ifndef _H_CHEMICAL_
#define _H_CHEMICAL_

class Chemical
{
 public:
  double m_mw, // molecular weight
    m_K_ow,	//partition coefficient between octanol and water
    m_pKa, // ionisation
    m_frac_non_ion, // fraction of solute non-ionised at pH 7.4
    m_frac_unbound; // fraction of solute unbound in a 2.7% albumin solution at pH 7.4
  char m_acid_base; // whether it's acid ('A') or base ('B')

 public:
  Chemical(void) {};
  ~Chemical(void) {};

  void Init(double, double, double, double, double, char);
};

#endif
