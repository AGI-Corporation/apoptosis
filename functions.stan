/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

// U = 0 if t < tau else 1
        // R   = Ro*np.exp(sigmu*t)
        // Qa  = kq*Ro/sigmu*(np.exp(sigmu*t)-1) - kq*Ro/sigmu*(np.exp(sigmu*(t-tau))-1) * U
        // Qc  = ( kq*Ro/(sigmu+kd)*(np.exp((sigmu+kd)*(t-tau))*np.exp(kd*tau) - np.exp(kd*tau))*np.exp(-kd*t) ) * U

/**
 * Total cell density at time t.
 *
 * This function provides an analytically simplified solution to the ODE system,
 * reducing the number of exp() calls from 6 to at most 3.
 */
real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  if (t < td) {
    return R0 * (mu / sm * exp(sm * t) - kq / sm);
  } else {
    real mu_over_sm = mu / sm;
    real kq_over_sm_plus_kd = kq / (sm + kd);
    real term2_coeff = (kq * kd) / (sm * (sm + kd));
    real inv_exp_sm_td = exp(-sm * td);
    real combined_coeff = mu_over_sm - term2_coeff * inv_exp_sm_td;
    return R0 * (combined_coeff * exp(sm * t) - kq_over_sm_plus_kd * exp(-kd * (t - td)));
  }
}

// These are kept for backward compatibility or individual component access if needed,
// but yt() now uses a more efficient direct calculation.
real Rt(real t, real R0, real sm){
  return R0 * exp(sm * t);
}

real Qat(real t, real R0, real sm, real kq, real td){
  real U = t < td ? 0 : 1;
  return kq * R0 / sm * (exp(sm * t) - 1)
    - kq * R0 / sm * (exp(sm * (t - td)) - 1) * U;
}

real Qct(real t, real R0, real sm, real kq, real td, real kd){
  real U = t < td ? 0 : 1;
  return U
    * kq * R0 / (sm + kd)
    * (exp((sm + kd) * (t - td)) * exp(kd * td) - exp(kd * td))
    * exp(-kd * t);
}

/* 
   Functions for solving the system numerically - use these to check that the
   analytic solution works.

*/


vector dsdt(real t, vector y, real R0, real sm, real kq, real td, real kd){
  vector[4] flux = [(sm + kq) * y[1],
                    kq * y[1],
                    t < td ? 0 : kq * Rt(t - td, R0, sm),
                    kd * y[3]]';
  return [flux[1]-flux[2], flux[2]-flux[3], flux[3]-flux[4]]';
}
real yt_num(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  real out = sum(ode_rk45(dsdt, [R0, 0, 0]', 0, {t}, R0, sm, kq, td, kd)[1]);
  return out > 0 ? out : 0.00001;
}
