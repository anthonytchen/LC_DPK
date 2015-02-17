/* The header file for chemical being applied to skin
   It should include all the properties (descriptors) of the chemical
 */
#ifndef _H_CHEMICAL_
#define _H_CHEMICAL_

class Chemical
{
 public:
  dobule m_mw, // molecular weight
    m_K_ow;	//partition coefficient between octanol and water

 public:
  Chemical(void) {};
  ~Chemical(void) {};

  void Init(double, double);
};

#endif
