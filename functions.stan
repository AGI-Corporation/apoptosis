/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

// U = 0 if t < tau else 1
        // R   = Ro*np.exp(sigmu*t)
        // Qa  = kq*Ro/sigmu*(np.exp(sigmu*t)-1) - kq*Ro/sigmu*(np.exp(sigmu*(t-tau))-1) * U
        // Qc  = ( kq*Ro/(sigmu+kd)*(np.exp((sigmu+kd)*(t-tau))*np.exp(kd*tau) - np.exp(kd*tau))*np.exp(-kd*t) ) * U

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

/**
 * Analytic solution for total cell density yt.
 * Optimized to minimize exp() calls and simplified algebraically.
 */
real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  // Use a small epsilon for sm near zero to avoid division by zero
  // although in practice mu < kq usually.
  real sm_eff = (abs(sm) < 1e-9) ? (sm > 0 ? 1e-9 : -1e-9) : sm;

  real term1 = (mu * R0 / sm_eff) * exp(sm_eff * t);
  real res;

  if (t < td) {
    res = term1 - (kq * R0 / sm_eff);
  } else {
    real smkd = sm_eff + kd;
    real smkd_eff = (abs(smkd) < 1e-9) ? (smkd > 0 ? 1e-9 : -1e-9) : smkd;
    real common = kq * R0 / smkd_eff;
    res = term1 - (kd / sm_eff) * common * exp(sm_eff * (t - td)) - common * exp(kd * (td - t));
  }
  // Ensure numerical stability and positivity for log-likelihood
  return res > 1e-9 ? res : 1e-9;
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
