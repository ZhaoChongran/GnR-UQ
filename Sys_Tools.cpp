#include "Sys_Tools.hpp"

void SYS_T::print(const double * const &array, const int &size)
{
    for(int ii=0; ii<size; ++ii) cout<<array[ii]<<'\n';
    cout<<endl;
}


int SYS_T::get_tn0(const int &n_t, const int &num_DL)
{
  if(n_t <= num_DL - 1) return 0;
  else return n_t - num_DL + 1;
}


void SYS_T::get_Lk( double * const &Lk, const double &Lt, const double &Lz,
    const double * const &alpha, const int &size )
{
  for(int ii=0; ii<size; ++ii)
  {
    const double a = sin(alpha[ii]);
    const double b = cos(alpha[ii]);
    Lk[ii] = sqrt( Lt * Lt * a * a + Lz * Lz * b * b );
  }
}



// EOF
