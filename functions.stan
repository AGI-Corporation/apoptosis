/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

/**
 * Total cell density at time t.
 *
 * @param t Time
 * @param R0 Initial density
 * @param mu Growth rate of resistant cells
 * @param kq Rate of transition to Classical Apoptosis Resistance
 * @param td Delay time before onset of death
 * @param kd Death rate
 * @return Total cell density
 */
real yt(real t, real R0, real mu, real kq, real td, real kd) {
  real sm = mu - kq;
  if (t < td) {
    return fmax(R0 / sm * (mu * expm1(sm * t) + sm), 1e-9);
  } else {
    real sm_td = sm * (t - td);
    return fmax(R0 * (
      (mu / sm) * exp(sm * t)
      - (kq * kd / (sm * (sm + kd))) * exp(sm_td)
      - (kq / (sm + kd)) * exp(-kd * (t - td))
    ), 1e-9);
  }
}

/* 
   Functions for solving the system numerically - use these to check that the
   analytic solution works.

*/


vector dsdt(real t, vector y, real R0, real sm, real kq, real td, real kd){
  vector[3] out;
  out[1] = sm * y[1];
  out[2] = kq * y[1] - (t < td ? 0 : kq * R0 * exp(sm * (t - td)));
  out[3] = (t < td ? 0 : kq * R0 * exp(sm * (t - td))) - kd * y[3];
  return out;
}

real yt_num(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  real out = sum(ode_rk45(dsdt, [R0, 0, 0]', 0, {t}, R0, sm, kq, td, kd)[1]);
  return out > 0 ? out : 1e-9;
}
