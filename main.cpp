// ==================================================================
// main.cpp
// This is my new main driver for the 1D growth & remodeling simulation
//
// Author: Ju Liu
// Date: June 3rd 2016
// ==================================================================
#include <fstream>

#include "Model_wall.hpp"
#include "Time_solver.hpp"

double test(const double &K_c1, const double &K_c2, const double &K_m1, const double &K_m2,
    const double &G_ch, const double &G_mh, const double &G_et, const double &G_ez);

int main()
{
  double K_c1 = 1.0;
  double K_c2 = 1.0;
  double K_m1 = 1.0;
  double K_m2 = 1.0;
  double G_ch = 1.07;
  double G_mh = 1.25;
  double G_et = 1.5;
  double G_ez = 1.5;
  test(K_c1, K_c2, K_m1, K_m2, G_ch, G_mh, G_et, G_ez);
}

double test(const double &K_c1, const double &K_c2, const double &K_m1, const double &K_m2,
    const double &G_ch, const double &G_mh, const double &G_et, const double &G_ez )
{
  const double pi = atan(1) * 4;

  // ----------- Time Solver ---------------
  const int steps_pday = 10;
  const int lifespan = 5000;
  const int simlength = 1000;
  const int ref_days = 300; 
  Time_solver * tsolver = new Time_solver(steps_pday, lifespan, simlength);

  tsolver->print_timeinfo();
  // ---------------------------------------

  // ----------- Wall Object ---------------
  const double dP = 1.0;
  const double dQ = 1.3;
  Model_wall * wall = new Model_wall(pi, dP, dQ, tsolver->get_num_t(),
      tsolver->get_num_DL(), tsolver->get_dt(), K_c1, K_c2, K_m1, K_m2, G_ch, G_mh, G_et, G_ez);

  wall->print_fluid_properties();

  wall->print_solid_properties();

  wall->check_initial_parameters();

  const double alpha_ckh[4] = {0.0, 0.5*pi, 0.25*pi, 0.75*pi};
  wall->check_initial_angle(alpha_ckh);

  wall->check_initial_stress();
  // ---------------------------------------


  // ------------ Nonlinear Solver ---------
  const double Max_error_a = 2.0e-9, Max_error_m = 1.0e-9;
  const int Max_it = 90;
  int num_it1 = 0, num_it2 = 0; 
  double beta = 0.3, tol_a = 100.0, tol_m = 100.0; 
  // ---------------------------------------


  // ----- Working variable in solver -----
  double L_z = 1.0, L_t;
  double a_act_p  = wall->get_a_M();
  double da_act_p = 0.0;
  const double k_act = 1.0 / 20.0; 

  double dwdLt_c, dwdLz_c, dwdLt_m, dwdLt_e;
  double ddwddLt_c, ddwddLt_m, ddwddLt_e;
  double M_ck[4] = {0.0, 0.0, 0.0, 0.0};
  double M_m = 0.0, M_c = 0.0;
  double Lc_k[4] = {0.0, 0.0, 0.0, 0.0};
  double Lc_kn[4] = {0.0, 0.0, 0.0, 0.0};
  double Lc_k_tau[4] = {0.0, 0.0, 0.0, 0.0};
  double alpha_tau[4] = {0.0, 0.5*pi, 0.0, 0.0};
  int tn0;
  double wt;
  double Lt_tau;
  double Lz_tau = 1.0;
  double Lm_n;
  double C_t, dC_t, T_act, dT_act;
  double Fa, dFa_da;
  double h_h[tsolver->get_num_t()];
  double total_M[tsolver->get_num_t()];
  double tau_w[tsolver->get_num_t()];
  double radius_t[tsolver->get_num_t()];
  h_h[0] = 0.0176066;
  total_M[0] = 0.0184869;
  tau_w[0] = 50.6;
  radius_t[0] = 0.142;
  // --------------------------------------

  // ----- Prepare file for recording -----
  //ofstream outfile( "results", ofstream::out | ofstream::trunc );

  //if(!outfile)
  //{
  //cerr<<"Error: unable to open file to record results. \n";
  //exit(EXIT_FAILURE);
  //}

  // --------------------------------------
  for( int n_t = 1; n_t < tsolver->get_num_t(); ++n_t )
  {
    double t = n_t * tsolver->get_dt();
   
    double P = wall->get_P(t,ref_days); 
    double Q = wall->get_Q(t,ref_days);

    // ! Warning : This predictor is not a standard one
    wall->predictor(n_t, 0.1);

    double a_t = wall->get_Da(n_t);

    // ! Warning : This predictor is not a standard one 
    double a_act = a_act_p + 0.5 * tsolver->get_dt() *
      (da_act_p + k_act * (a_t - (a_act_p + tsolver->get_dt()*da_act_p)));

    tol_m = 100.0; num_it1 = 0;

    tn0 = SYS_T::get_tn0(n_t, tsolver->get_num_DL()); 

    //outfile<<t<<'\t'<<P<<'\t'<<Q<<'\t';

    while( (tol_m > Max_error_m) && (num_it1 < Max_it) )
    {
      num_it1 += 1; 
      tol_a = 100.0; num_it2 = 0;
      while( (tol_a > Max_error_a) && (num_it2 < Max_it) )
      {
        num_it2 += 1;

        double tau_w = 4.0 * wall->get_mu() * Q / (pi*a_t*a_t*a_t);

        L_t = a_t / wall->get_a_M();

        // Update the angle based on L_t
        wall->set_Dalpha(n_t, L_t, L_z); 

        // calculate the stress and d_stress 
        // -- stress
        dwdLt_c = 0.0;
        dwdLz_c = 0.0;
        dwdLt_m = 0.0;
        ddwddLt_c = 0.0;
        ddwddLt_m = 0.0;

        // -- mass initialization
        for(int ii=0; ii<4; ++ii) M_ck[ii] = 0.0;

        M_m = 0.0;

        // Calculate initial mass/energy with degradation
        if( n_t <= tsolver->get_num_DL() )
        {
          wall->get_Lk(Lc_k, L_t, L_z, alpha_ckh);

          for(int ii=0; ii<4; ++ii)
          {
            M_ck[ii]   = wall->get_M_ck(ii, n_t);
            dwdLt_c   += wall->get_dwdLt_c(M_ck[ii], L_t, L_z, 
                alpha_ckh[ii], Lc_k[ii], 1.0);
            dwdLz_c   += wall->get_dwdLz_c(M_ck[ii], L_t, L_z, 
                alpha_ckh[ii], Lc_k[ii], 1.0);
            ddwddLt_c += wall->get_ddwddLt_c(M_ck[ii], L_t, L_z, 
                alpha_ckh[ii], Lc_k[ii], 1.0);
          }

          M_m        = wall->get_M_m(n_t);
          dwdLt_m   += wall->get_dwdLt_m(M_m, L_t, 1.0);
          ddwddLt_m += wall->get_ddwddLt_m(M_m, L_t, 1.0);
        }


        // Calculate viscoelasticity
        for(int n_tau = tn0; n_tau <= n_t; ++n_tau)
        {
          if(n_tau == tn0 || n_tau == n_t) wt = 0.5 * tsolver->get_dt();
          else wt = tsolver->get_dt();

          alpha_tau[2] = wall->get_Dalpha(n_tau);
          alpha_tau[3] = 2.0 * pi - alpha_tau[2];

          Lt_tau = wall->get_Da(n_tau) / wall->get_a_M();      

          wall->get_Lk(Lc_k_tau, Lt_tau, Lz_tau, alpha_tau);

          // This following is from the old code !!!
          wall->get_Lk(Lc_k, L_t, L_z, alpha_tau); 

          for(int ii=0; ii<4; ++ii)
          {
            Lc_kn[ii] = wall->get_Gch() * Lc_k[ii] / Lc_k_tau[ii];
            if(Lc_kn[ii] <= wall->get_y_Lkn())
            {
              const double new_cmass = wall->get_mc_tau(n_t, n_tau, 
                  ii, tsolver->get_dt(), wt);
              M_ck[ii] += new_cmass;
              dwdLt_c  += wall->get_dwdLt_c(new_cmass, L_t, L_z, 
                  alpha_tau[ii], Lc_k[ii], Lc_k_tau[ii]);
              dwdLz_c  += wall->get_dwdLz_c(new_cmass, L_t, L_z, 
                  alpha_tau[ii], Lc_k[ii], Lc_k_tau[ii]);
              ddwddLt_c += wall->get_ddwddLt_c(new_cmass, L_t, L_z,
                  alpha_tau[ii], Lc_k[ii], Lc_k_tau[ii]);
            }
          }

          Lm_n = wall->get_Gmh() * L_t / Lt_tau;

          if(Lm_n <= wall->get_y_Lmn())
          {
            const double new_mmass = wall->get_mm_tau(n_t, n_tau,
                tsolver->get_dt(), wt);
            M_m += new_mmass;

            dwdLt_m += wall->get_dwdLt_m(new_mmass, L_t, Lt_tau);

            ddwddLt_m += wall->get_ddwddLt_m(new_mmass, L_t, Lt_tau);
          }
        }

        dwdLt_e = wall->get_dwdLt_e( L_t * wall->get_Get(), 
            L_z * wall->get_Gez() );
        ddwddLt_e = wall->get_ddwddLt_e( L_t * wall->get_Get(), 
            L_z * wall->get_Gez() );

        // calculate active stress
        const double L_m_act = a_t / a_act;

        C_t    = wall->get_C_t( tau_w );
        dC_t   = wall->get_dC_t( Q, a_t );
        T_act  = wall->get_T_act( M_m, L_m_act, L_t*L_z, C_t );
        dT_act = wall->get_dT_act( M_m, L_m_act, a_t, L_z, a_act,
            L_t * L_z, C_t, dC_t );

        // calculate a_t and related data
        Fa = ( (dwdLt_c + dwdLt_m + dwdLt_e) / L_z ) + T_act - P * a_t;
        dFa_da = ( (ddwddLt_c + ddwddLt_m + ddwddLt_e) / (L_z*wall->get_a_M()) ) 
          + dT_act - P;

        a_t -= beta * Fa / dFa_da;

        a_act = (a_act_p + 0.5*tsolver->get_dt()*(da_act_p + k_act*a_t))
          /(1.0 + 0.5*tsolver->get_dt()*k_act);

        L_t = a_t / wall->get_a_M();

        tol_a = wall->l2error_a(a_t, n_t);

        // set in wall object
        wall->set_Da(n_t, a_t);
        wall->set_Dalpha(n_t, L_t, L_z); 
      } // end while tol_a > Max_error && num_it2 < Max_it

      if(num_it2 == Max_it) beta = 0.1;

      M_c = 0.0;

      for(int ii=0; ii<4; ++ii) M_c += M_ck[ii];

      double error_c, error_bottom_c, error_m, error_bottom_m;

      // calculate the new mass for collagen and muscle
      wall->update_m_c( n_t, L_t, L_z, dwdLt_c, dwdLz_c, M_c, C_t,
          error_c, error_bottom_c );

      wall->update_m_m( n_t, L_t, L_z, dwdLt_m, T_act, M_m, C_t,
          error_m, error_bottom_m );

      tol_m = sqrt((error_c + error_m) / (error_bottom_c + error_bottom_m));
    } // end while tol_m > Max_error && num_it1 < Max_it

    a_act_p = a_act;
    da_act_p = k_act * (a_t - a_act);

    wall->set_Dalpha(n_t, L_t, L_z);

    double M_e = wall->get_M_eh();
    total_M[n_t] = M_c + M_e + M_m;
    h_h[n_t] = total_M[n_t] / (wall->get_rho_s() * L_t * L_z);
    tau_w[n_t] = 4.0 * wall->get_mu() * Q / (pi*a_t*a_t*a_t);
    radius_t[n_t] = a_t;

    //outfile<<a_t<<'\t'<<h_h[n_t]<<'\t'<<M_c<<'\t'<<M_m<<'\t'<<M_e<<'\t'<<total_M<<'\t';
    //outfile<<wall->get_Dalpha(n_t);
    //outfile<<endl;

    // update DQ2 from tn0 to the current time step 
    for(int ii=tn0; ii<=n_t; ++ii)
    {
      wall->update_DQ2_c(ii, L_t, L_z, tsolver->get_dt());
      wall->update_DQ2_m(ii, L_t, L_z, tsolver->get_dt());
    }
    //cout<<"Time t= "<<t<<'\t';
    //cout<<"num_it1 = "<<num_it1<<'\t'<<"tol_a = "<<tol_a<<'\t';
    //cout<<"L_t = "<<L_t<<'\t'<<"h_h = "<<h_h<<'\t';
    //cout<<"total_M = "<<total_M<<'\t';
    //cout<<endl; 
  }

  //outfile.close();
  double new_a_h = radius_t[tsolver->get_num_t()-1];
  double new_h_h = h_h[tsolver->get_num_t()-1];
  double new_tau_w = tau_w[tsolver->get_num_t()-1];
  double new_M = total_M[tsolver->get_num_t()-1];
  const double tol_a_t = abs(radius_t[ref_days*steps_pday-1] - radius_t[0]) / radius_t[0];
  const double tol_h = abs(h_h[ref_days*steps_pday-1] - h_h[0]) / h_h[0];
  const double tol_M = abs(total_M[ref_days*steps_pday-1] - total_M[0]) / total_M[0];
  const double tol_tau_w = abs(tau_w[ref_days*steps_pday-1] - tau_w[0]) / tau_w[0];

  for (int n_t = 1; n_t < tsolver->get_num_t(); n_t++)
  {
    double t = n_t * tsolver->get_dt();
    if( (abs(radius_t[n_t] / new_a_h - 1.0) <= tol_a_t) && t > ref_days)
    {
      cout << "Time reaching homeostatic inner radius is " << t << endl; 
      break;
    }
  }
  for (int n_t = 1; n_t < tsolver->get_num_t(); n_t++)
  {
    double t = n_t * tsolver->get_dt();
    if( (abs(h_h[n_t] / new_h_h - 1.0) <= tol_h) && t > ref_days)
    {
      cout << "Time reaching homeostatic thickness is " << t << endl;
      break;
    }
  }
  for (int n_t = 1; n_t < tsolver->get_num_t(); n_t++)
  {
    double t = n_t * tsolver->get_dt(); 
    if( (abs(tau_w[n_t] / new_tau_w - 1.0) <= tol_tau_w) && t > ref_days)
    {
      cout << "Time reaching homeostatic WSS is " << t << endl; 
      break;
    }
  }
  for (int n_t = 1; n_t < tsolver->get_num_t(); n_t++)
  {
    double t = n_t * tsolver->get_dt();
    if( (abs(total_M[n_t] / new_M - 1.0) <= tol_M) && t > ref_days)
    {
      cout << "Time reaching homeostatic total mass is " << t << endl;
      break;
    }
  }

  delete wall; delete tsolver;
}
// EOF
