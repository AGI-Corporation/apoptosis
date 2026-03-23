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

/**
 * Total cell density yt at time t.
 *
 * This is an algebraically simplified version of the analytic solution:
 * yt = Rt + Qat + Qct
 *
 * Performance optimization:
 * - Reduces exp() calls from up to 6 to 3 (for t >= td) or 1 (for t < td).
 * - Adds a 0.00001 floor for numerical stability in log-likelihoods.
 */
real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  real out;

  if (t < td) {
    // For t < td:
    // Rt = R0 * exp(sm * t)
    // Qat = kq * R0 / sm * (exp(sm * t) - 1)
    // Qct = 0
    // yt = R0 * exp(sm * t) * (1 + kq/sm) - kq * R0 / sm
    // since 1 + kq/sm = (sm + kq)/sm = mu/sm
    out = (R0 / sm) * (mu * exp(sm * t) - kq);
  } else {
    // For t >= td:
    // Rt = R0 * exp(sm * t)
    // Qat = kq * R0 / sm * (exp(sm * t) - exp(sm * (t - td)))
    // Qct = kq * R0 / (sm + kd) * (exp(sm * (t - td)) - exp(kd * td - kd * t))
    // yt = R0 * exp(sm * t) + (kq * R0 / sm) * (exp(sm * t) - exp(sm * (t - td)))
    //      + (kq * R0 / (sm + kd)) * (exp(sm * (t - td)) - exp(-kd * (t - td)))

    real e_smt = exp(sm * t);
    real e_sm_t_td = exp(sm * (t - td));
    real e_kd_t_td = exp(-kd * (t - td));

    out = R0 * e_smt + (kq * R0 / sm) * (e_smt - e_sm_t_td)
          + (kq * R0 / (sm + kd)) * (e_sm_t_td - e_kd_t_td);
  }

  return out > 0.00001 ? out : 0.00001;
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
