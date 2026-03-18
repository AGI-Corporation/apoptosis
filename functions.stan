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
  if (t < td) {
    return (kq * R0 / sm) * (exp(sm * t) - 1);
  } else {
    return (kq * R0 / sm) * (exp(sm * t) - exp(sm * (t - td)));
  }
}

real Qct(real t, real R0, real sm, real kq, real td, real kd){
  if (t < td) {
    return 0;
  } else {
    return (kq * R0 / (sm + kd)) * (exp(sm * (t - td)) - exp(-kd * (t - td)));
  }
}

/**
 * Total cell density at time t.
 *
 * Optimized to reduce the number of exp() calls and simplify algebraic terms.
 * Reduces exp() calls from ~7 to 1 (for t < td) or 3 (for t >= td).
 */
real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  if (t < td) {
    return (R0 / sm) * (mu * exp(sm * t) - kq);
  } else {
    real est = exp(sm * t);
    real est_td = exp(sm * (t - td));
    real ekdt_td = exp(-kd * (t - td));
    // Simplified: Rt + Qat + Qct
    return (mu * R0 / sm) * est
           - (kq * kd * R0 / (sm * (sm + kd))) * est_td
           - (kq * R0 / (sm + kd)) * ekdt_td;
  }
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
