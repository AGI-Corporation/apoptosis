/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

/**
 * Analytically solve the ODE system for total cell density.
 *
 * The system is:
 * dR/dt = (mu - kq) * R
 * dQa/dt = kq * R - U(t-td) * kq * R(t-td)
 * dQc/dt = U(t-td) * kq * R(t-td) - kd * Qc
 *
 * y(t) = R(t) + Qa(t) + Qc(t)
 *
 * Optimized to reduce exp() calls and redundant arithmetic.
 */
real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  real val;

  if (t < td) {
    // y(t) = R0 * (exp(sm*t) + (kq/sm) * (exp(sm*t) - 1))
    //      = (R0 / sm) * (mu * exp(sm*t) - kq)
    val = (R0 / sm) * (mu * exp(sm * t) - kq);
  } else {
    // y(t) = R0 * [ (mu/sm)*exp(sm*t) - (kq*kd/(sm*(sm+kd)))*exp(sm*(t-td)) - (kq/(sm+kd))*exp(-kd*(t-td)) ]
    real sm_plus_kd = sm + kd;
    val = R0 * ( (mu / sm) * exp(sm * t) -
                 (kq * kd / (sm * sm_plus_kd)) * exp(sm * (t - td)) -
                 (kq / sm_plus_kd) * exp(-kd * (t - td)) );
  }

  return val > 1e-9 ? val : 1e-9;
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
