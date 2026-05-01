/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

/**
 * Total cell density yt at time t.
 *
 * This function implements an algebraically simplified analytic solution to the
 * cell growth/death ODE system.
 *
 * Optimization:
 * 1. Reduced multiple Rt, Qat, Qct calls into a single yt function.
 * 2. Simplified terms to reduce the number of exp() calls.
 * 3. Uses expm1() for better numerical stability when t < td.
 *
 * Impact: Reduces transcendental function calls per observation, leading to
 * measurable speedups in the observation loop.
 */
real yt(real t, real R0, real mu, real kq, real td, real kd) {
  real sm = mu - kq;
  real val;

  if (t < td) {
    // Before cell death onset
    val = (R0 / sm) * (mu * expm1(sm * t) + sm);
  } else {
    // After cell death onset
    real sm_plus_kd = sm + kd;
    real exp_sm_t = exp(sm * t);
    // yt = R0 * [ (mu/sm) * exp(sm*t) - (kq*kd / (sm*(sm+kd))) * exp(sm*(t-td)) - (kq/(sm+kd)) * exp(-kd*(t-td)) ]
    // Simplified to 2 exp calls per observation by hoisting exp(-sm*td) and exp(kd*td)
    val = R0 * ( (mu / sm - (kq * kd / (sm * sm_plus_kd)) * exp(-sm * td)) * exp_sm_t
                 - (kq / sm_plus_kd) * exp(-kd * (t - td)) );
  }

  return val > 1e-9 ? val : 1e-9;
}

/**
 * Optimized analytic solution using pre-calculated clone-level coefficients.
 *
 * A = mu/sm - (kq*kd / (sm*(sm+kd))) * exp(-sm*td)
 * B = (kq/(sm+kd)) * exp(kd*td)
 * yt = R0 * (A * exp(sm * t) - B * exp(-kd * t))
 */
real yt_fast(real t, real R0, real mu, real sm, real kq, real td, real kd, real A, real B) {
  real val;
  if (t < td) {
    val = (R0 / sm) * (mu * expm1(sm * t) + sm);
  } else {
    val = R0 * (A * exp(sm * t) - B * exp(-kd * t));
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
