/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

/**
 * Total cell density at time t.
 *
 * This is an algebraically simplified analytic solution to the ODE system
 * described in the dsdt function below.
 *
 * @param t Time
 * @param R0 Initial cell density
 * @param mu Growth rate
 * @param kq Rate of transition to apoptosis
 * @param td Time delay before death
 * @param kd Rate of death
 * @return Total cell density (sum of live and apoptotic cells)
 */
real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  real out;

  if (t < td) {
    // Simplified: Rt(t) + Qat(t)
    // = R0 * exp(sm * t) + (kq * R0 / sm) * (exp(sm * t) - 1)
    // = (R0 / sm) * (mu * expm1(sm * t) + sm)
    out = (R0 / sm) * (mu * expm1(sm * t) + sm);
  } else {
    // Simplified: Rt(t) + Qat(t) + Qct(t)
    // = (R0 * mu / sm) * exp(sm * t)
    //   - (R0 * kq * kd / (sm * (sm + kd))) * exp(sm * (t - td))
    //   - (R0 * kq / (sm + kd)) * exp(-kd * (t - td))
    real sm_plus_kd = sm + kd;
    out = R0 * ( (mu / sm) * exp(sm * t)
                - (kq * kd / (sm * sm_plus_kd)) * exp(sm * (t - td))
                - (kq / sm_plus_kd) * exp(-kd * (t - td)) );
  }

  // Ensure numerical stability and avoid log(0)
  return out > 1e-9 ? out : 1e-9;
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
