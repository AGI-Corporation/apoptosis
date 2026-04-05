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
 * Total cell density at time t, analytically solved.
 *
 * Optimization:
 * 1. Algebraic simplification to minimize expensive exp() calls.
 * 2. Implements a floor of 1e-9 to ensure numerical stability for log-likelihood.
 */
real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  real val;

  if (t < td) {
    // simplified Rt + Qat (Qct is 0)
    // R0 * exp(sm*t) + kq*R0/sm * (exp(sm*t) - 1)
    // = R0 * (exp(sm*t) + kq/sm*exp(sm*t) - kq/sm)
    // = R0 * (exp(sm*t) * (1 + kq/sm) - kq/sm)
    // 1 + kq/sm = (sm + kq)/sm = mu/sm
    val = R0 * (mu / sm * exp(sm * t) - kq / sm);
  } else {
    real exp_sm_t = exp(sm * t);
    real exp_sm_t_minus_td = exp(sm * (t - td));
    real exp_neg_kd_t_minus_td = exp(-kd * (t - td));

    // Rt + Qat (for t >= td)
    // R0 * exp(sm*t) + kq*R0/sm * (exp(sm*t) - 1) - kq*R0/sm * (exp(sm*(t-td)) - 1)
    // = R0 * exp(sm*t) + kq*R0/sm * (exp(sm*t) - exp(sm*(t-td)))
    // = R0 * (exp(sm*t) * (1 + kq/sm) - kq/sm * exp(sm*(t-td)))
    // = R0 * (mu/sm * exp(sm*t) - kq/sm * exp(sm*(t-td)))

    // Qct (for t >= td)
    // kq * R0 / (sm + kd) * (exp((sm + kd) * (t - td)) * exp(kd * td) - exp(kd * td)) * exp(-kd * t)
    // = kq * R0 / (sm + kd) * (exp(sm*(t-td) + kd*(t-td) + kd*td) - exp(kd*td)) * exp(-kd*t)
    // = kq * R0 / (sm + kd) * (exp(sm*(t-td) + kd*t) - exp(kd*td)) * exp(-kd*t)
    // = kq * R0 / (sm + kd) * (exp(sm*(t-td)) - exp(kd*td - kd*t))
    // = kq * R0 / (sm + kd) * (exp(sm*(t-td)) - exp(-kd*(t-td)))

    val = R0 * (mu / sm * exp_sm_t
                - kq / sm * exp_sm_t_minus_td
                + kq / (sm + kd) * (exp_sm_t_minus_td - exp_neg_kd_t_minus_td));
  }

  return fmax(1e-9, val);
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
