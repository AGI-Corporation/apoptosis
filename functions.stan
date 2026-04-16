/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

/**
 * Total cell density at time t.
 *
 * Optimized analytic solution to minimize exp() calls and improve numerical stability.
 */
real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  real val;

  if (t < td){
    // Simplified: R0 * exp(sm * t) + (kq * R0 / sm) * (exp(sm * t) - 1)
    // Uses expm1 for better precision when sm * t is small.
    val = R0 * (exp(sm * t) + (kq / sm) * expm1(sm * t));
  } else {
    // Simplified combination of Rt, Qat, and Qct:
    // val = R0 * [ (mu/sm)*exp(sm*t) - (kq*kd/(sm*(sm+kd)))*exp(sm*(t-td)) - (kq/(sm+kd))*exp(-kd*(t-td)) ]
    real sm_plus_kd = sm + kd;
    real E1 = exp(sm * t);
    real E2 = exp(sm * (t - td));
    real E3 = exp(-kd * (t - td));

    val = R0 * ( (mu / sm) * E1 - (kq * kd / (sm * sm_plus_kd)) * E2 - (kq / sm_plus_kd) * E3 );
  }

  // Ensure positive return value for log-likelihood stability
  return val > 1e-9 ? val : 1e-9;
}

/* 
   Functions for solving the system numerically - use these to check that the
   analytic solution works.

*/

vector dsdt(real t, vector y, real R0, real sm, real kq, real td, real kd){
  // Rt(t) = R0 * exp(sm * t)
  // Fluxes between compartments:
  // R -> Qa : kq * R(t)
  // Qa -> Qc : kq * R(t-td) if t > td else 0
  // Qc -> ... : kd * Qc

  real Rt_val = R0 * exp(sm * t);
  vector[4] flux;
  flux[1] = (sm + kq) * y[1];
  flux[2] = kq * y[1];
  flux[3] = t < td ? 0 : kq * R0 * exp(sm * (t - td));
  flux[4] = kd * y[3];

  return [flux[1]-flux[2], flux[2]-flux[3], flux[3]-flux[4]]';
}

real yt_num(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  array[1] vector[3] sol = ode_rk45(dsdt, [R0, 0, 0]', 0, {t}, R0, sm, kq, td, kd);
  real out = sum(sol[1]);
  return out > 0 ? out : 1e-9;
}
