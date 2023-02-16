#ifndef SYS_TOOLS_HPP
#define SYS_TOOLS_HPP

#include <iostream>
#include <cmath>

using namespace std;

namespace SYS_T
{
  void print(const double * const &array, const int &size);

  int get_tn0(const int &n_t, const int &num_DL);

  void get_Lk( double * const &Lk, const double &Lt, const double &Lz,
    const double * const &alpha, const int &size );



}


#endif
