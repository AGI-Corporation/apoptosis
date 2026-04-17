/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

/**
 * Total cell density at time t.
 *
 * This function provides an optimized analytic solution for the ODE system:
 * dR/dt = (mu - kq) * R
 * dQa/dt = kq * R - kq * R(t - td) * U(t - td)
 * dQc/dt = kq * R(t - td) * U(t - td) - kd * Qc
 * where U is the Heaviside step function.
 *
 * The solution is simplified to minimize expensive exp() calls and improve numerical stability.
 */
real yt(real t, real R0, real mu, real kq, real td, real kd) {
  real sm = mu - kq;
  real R = R0 * exp(sm * t);

  if (t < td) {
    // For t < td, only R and Qa are present (and Qa = kq * R0 * (exp(sm*t) - 1) / sm)
    // Simplified: yt = R + Qa = R0 * exp(sm*t) * (1 + kq/sm) - kq*R0/sm
    // Using expm1 for better precision when sm*t is small.
    return R0 + (mu / sm) * R0 * expm1(sm * t);
  } else {
    // For t >= td, R, Qa, and Qc are present.
    // Rt = R0 * exp(sm * t)
    // Qat = (kq * R0 / sm) * (exp(sm * t) - exp(sm * (t - td)))
    // Qct = (kq * R0 / (sm + kd)) * (exp(sm * (t - td)) - exp(-kd * t + kd * td - sm * td))

    // Original Qct:
    // real U = t < td ? 0 : 1;
    // return U * kq * R0 / (sm + kd) * (exp((sm + kd) * (t - td)) * exp(kd * td) - exp(kd * td)) * exp(-kd * t);
    // = kq * R0 / (sm + kd) * (exp(sm*t - sm*td + kd*t - kd*td + kd*td - kd*t) - exp(kd*td - kd*t))
    // = kq * R0 / (sm + kd) * (exp(sm*(t - td)) - exp(-kd*(t - td)))

    real exp_sm_t = exp(sm * t);
    real exp_sm_t_minus_td = exp(sm * (t - td));
    real exp_neg_kd_t_minus_td = exp(-kd * (t - td));

    real Rt = R0 * exp_sm_t;
    real Qat = (kq * R0 / sm) * (exp_sm_t - exp_sm_t_minus_td);
    real Qct = (kq * R0 / (sm + kd)) * (exp_sm_t_minus_td - exp_neg_kd_t_minus_td);

    return Rt + Qat + Qct;
  }
}

/* 
   Functions for solving the system numerically - use these to check that the
   analytic solution works.

*/

vector dsdt(real t, vector y, real R0, real sm, real kq, real td, real kd){
  real Rt_lag = R0 * exp(sm * (t - td));
  vector[4] flux = [(sm + kq) * y[1],
                    kq * y[1],
                    t < td ? 0 : kq * Rt_lag,
                    kd * y[3]]';
  return [flux[1]-flux[2], flux[2]-flux[3], flux[3]-flux[4]]';
}

real yt_num(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  real out = sum(ode_rk45(dsdt, [R0, 0, 0]', 0, {t}, R0, sm, kq, td, kd)[1]);
  return out > 0 ? out : 1e-9;
}
