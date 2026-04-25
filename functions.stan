/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

// U = 0 if t < tau else 1
        // R   = Ro*np.exp(sigmu*t)
        // Qa  = kq*Ro/sigmu*(np.exp(sigmu*t)-1) - kq*Ro/sigmu*(np.exp(sigmu*(t-tau))-1) * U
        // Qc  = ( kq*Ro/(sigmu+kd)*(np.exp((sigmu+kd)*(t-tau))*np.exp(kd*tau) - np.exp(kd*tau))*np.exp(-kd*t) ) * U

/**
 * Total cell density yt at time t.
 *
 * This is the analytic solution to the ODE system:
 * dR/dt = (mu - kq) * R
 * dQa/dt = kq * R - kq * R(t - td) * U(t - td)
 * dQc/dt = kq * R(t - td) * U(t - td) - kd * Qc
 *
 * where sm = mu - kq.
 *
 * yt = R(t) + Qa(t) + Qc(t)
 *
 * Simplification:
 * t < td:
 * yt = R0 * exp(sm * t) + kq * R0 / sm * (exp(sm * t) - 1)
 *    = R0 * (exp(sm * t) + kq/sm * exp(sm * t) - kq/sm)
 *    = R0 * ( (1 + kq/sm) * exp(sm * t) - kq/sm )
 * Since 1 + kq/sm = (sm + kq)/sm = mu/sm:
 * yt = R0/sm * (mu * exp(sm * t) - kq)
 *    = R0/sm * (mu * (exp(sm * t) - 1) + mu - kq)
 *    = R0/sm * (mu * expm1(sm * t) + sm)
 *
 * t >= td:
 * yt = R(t) + Qa(t) + Qc(t)
 * R(t) = R0 * exp(sm * t)
 * Qa(t) = kq * R0 / sm * (exp(sm * t) - 1) - kq * R0 / sm * (exp(sm * (t - td)) - 1)
 *      = kq * R0 / sm * (exp(sm * t) - exp(sm * (t - td)))
 * Qc(t) = kq * R0 / (sm + kd) * (exp(sm * (t - td)) - exp(-kd * (t - td)))
 *
 * yt = R0 * exp(sm * t) + kq * R0 / sm * (exp(sm * t) - exp(sm * (t - td))) + kq * R0 / (sm + kd) * (exp(sm * (t - td)) - exp(-kd * (t - td)))
 *    = R0 * [ exp(sm * t) * (1 + kq/sm) - exp(sm * (t - td)) * (kq/sm - kq/(sm + kd)) - kq/(sm + kd) * exp(-kd * (t - td)) ]
 *    = R0 * [ mu/sm * exp(sm * t) - exp(sm * (t - td)) * (kq * (sm + kd) - kq * sm) / (sm * (sm + kd)) - kq/(sm + kd) * exp(-kd * (t - td)) ]
 *    = R0 * [ mu/sm * exp(sm * t) - exp(sm * (t - td)) * (kq * kd) / (sm * (sm + kd)) - kq/(sm + kd) * exp(-kd * (t - td)) ]
 */
real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  if (t < td) {
    return R0 / sm * (mu * expm1(sm * t) + sm);
  } else {
    real term1 = mu / sm * exp(sm * t);
    real term2 = (kq * kd) / (sm * (sm + kd)) * exp(sm * (t - td));
    real term3 = kq / (sm + kd) * exp(-kd * (t - td));
    return R0 * (term1 - term2 - term3);
  }
}

/* 
   Functions for solving the system numerically - use these to check that the
   analytic solution works.

*/

real Rt(real t, real R0, real sm){
  return R0 * exp(sm * t);
}

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
  return out > 0 ? out : 1e-9;
}
