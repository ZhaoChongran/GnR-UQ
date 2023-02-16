#include "Time_solver.hpp"

Time_solver::Time_solver( const int &in_ntpd, 
    const int &in_agemax, const int &in_simdays  )
  : num_step_per_day(in_ntpd), dt(1.0/(double) num_step_per_day),
    age_max(in_agemax), num_DL(num_step_per_day * age_max + 1),
    comp_t(in_simdays), num_t(comp_t * num_step_per_day + 1)
{}


Time_solver::~Time_solver()
{}


void Time_solver::print_timeinfo() const
{
  cout<<"Time stepping information: \n";
  cout<<" -- number of steps per day : "<<num_step_per_day<<endl;
  cout<<" -- time step size dt : "<<dt<<" day \n";
  cout<<" -- lifespan of collagen : "<<age_max<<" days \n";
  cout<<" -- number of time stpes for collagen/SM lifespan storage: ";
  cout<<num_DL<<endl;
  cout<<" -- number of days for the whole simulation : "<<comp_t<<" days \n";
  cout<<" -- number of steps for the whole simulation : "<<num_t<<endl;
}


// EOF
