#ifndef TIME_SOLVER_HPP
#define TIME_SOLVER_HPP

#include <iostream>

using namespace std;

class Time_solver
{
  public:
    Time_solver( const int &in_ntpd, const int &in_agemax,
      const int &in_simdays  );

    virtual ~Time_solver();
    
    int get_num_DL() const {return num_DL;}

    int get_num_t() const {return num_t;}

    double get_dt() const {return dt;}

    void print_timeinfo() const;

  private:
    const int num_step_per_day; // number of time integration per day
    const double dt;  // time step size = 1.0 / num_step_per_day
    const double age_max; // max lifespan for collagen and smooth muscle
    const int num_DL; // number of data for collagen and SM lifespan
                      // num_DL = num_step_per_day * age_max + 1
    const int comp_t; // number of days for the simulation
    const int num_t; //  number of time steps for the simulation
                     // num_t = comp_t * num_step_per_day + 1



};


#endif
