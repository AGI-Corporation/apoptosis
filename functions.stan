/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

// U = 0 if t < tau else 1
        // R   = Ro*np.exp(sigmu*t)
        // Qa  = kq*Ro/sigmu*(np.exp(sigmu*t)-1) - kq*Ro/sigmu*(np.exp(sigmu*(t-tau))-1) * U
        // Qc  = ( kq*Ro/(sigmu+kd)*(np.exp((sigmu+kd)*(t-tau))*np.exp(kd*tau) - np.exp(kd*tau))*np.exp(-kd*t) ) * U

real yt(real t, real R0, real mu, real kq, real td, real kd) {
  real sm = mu - kq;
  real exp_sm_t = exp(sm * t);
  real val;
  if (t < td) {
    val = R0 * (exp_sm_t + kq / sm * (exp_sm_t - 1));
  } else {
    real exp_sm_tmtd = exp(sm * (t - td));
    val = R0 * (exp_sm_t + kq / sm * (exp_sm_t - exp_sm_tmtd)
          + kq / (sm + kd) * (exp_sm_tmtd - exp(-kd * (t - td))));
  }
  return val > 0 ? val : 1e-9;
}

/* 
   Functions for solving the system numerically - use these to check that the
   analytic solution works.

*/


vector dsdt(real t, vector y, real R0, real sm, real kq, real td, real kd){
  vector[4] flux = [(sm + kq) * y[1],
                    kq * y[1],
                    t < td ? 0 : kq * R0 * exp(sm * (t - td)),
                    kd * y[3]]';
  return [flux[1]-flux[2], flux[2]-flux[3], flux[3]-flux[4]]';
}
real yt_num(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  real out = sum(ode_rk45(dsdt, [R0, 0, 0]', 0, {t}, R0, sm, kq, td, kd)[1]);
  return out > 0 ? out : 0.00001;
}
