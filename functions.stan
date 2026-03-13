/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

// U = 0 if t < tau else 1
        // R   = Ro*np.exp(sigmu*t)
        // Qa  = kq*Ro/sigmu*(np.exp(sigmu*t)-1) - kq*Ro/sigmu*(np.exp(sigmu*(t-tau))-1) * U
        // Qc  = ( kq*Ro/(sigmu+kd)*(np.exp((sigmu+kd)*(t-tau))*np.exp(kd*tau) - np.exp(kd*tau))*np.exp(-kd*t) ) * U

/**
 * Analytically solves the ODE system for total cell density.
 *
 * Performance optimizations:
 * - Simplified the analytic solution to reduce the number of exp() calls from 6 to at most 3.
 * - Used expm1() for better numerical stability and speed for small values.
 * - Avoided redundant recalculations of common terms.
 */
real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  real R0_kq = R0 * kq;
  if (t < td) {
    real sm_t = sm * t;
    return R0 * exp(sm_t) + R0_kq * (expm1(sm_t) / sm);
  } else {
    real sm_td = sm * td;
    real t_minus_td = t - td;
    real exp_sm_td = exp(sm_td);
    real y_td = R0 * exp_sm_td + R0_kq * (expm1(sm_td) / sm);
    real exp_sm_t_minus_td = exp(sm * t_minus_td);
    return y_td * exp_sm_t_minus_td
           + R0_kq * (exp_sm_t_minus_td - exp(-kd * t_minus_td)) / (sm + kd);
  }
}

/* 
   Functions for solving the system numerically - use these to check that the
   analytic solution works.

*/


vector dsdt(real t, vector y, real R0, real sm, real kq, real td, real kd){
  vector[4] flux = [(sm + kq) * y[1],
                    kq * y[1],
                    t < td ? 0 : kq * (R0 * exp(sm * (t - td))),
                    kd * y[3]]';
  return [flux[1]-flux[2], flux[2]-flux[3], flux[3]-flux[4]]';
}
real yt_num(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  real out = sum(ode_rk45(dsdt, [R0, 0, 0]', 0, {t}, R0, sm, kq, td, kd)[1]);
  return out > 0 ? out : 0.00001;
}
