/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

// U = 0 if t < tau else 1
        // R   = Ro*np.exp(sigmu*t)
        // Qa  = kq*Ro/sigmu*(np.exp(sigmu*t)-1) - kq*Ro/sigmu*(np.exp(sigmu*(t-tau))-1) * U
        // Qc  = ( kq*Ro/(sigmu+kd)*(np.exp((sigmu+kd)*(t-tau))*np.exp(kd*tau) - np.exp(kd*tau))*np.exp(-kd*t) ) * U

/**
 * Total cell density yt = Rt + Qat + Qct
 * sm = mu - kq
 *
 * Case 1: t < td
 * yt = R0 * exp(sm*t) + (kq * R0 / sm) * (exp(sm*t) - 1)
 *    = R0 * exp(sm*t) + (kq * R0 / sm) * exp(sm*t) - (kq * R0 / sm)
 *    = R0 * exp(sm*t) * (1 + kq/sm) - (kq * R0 / sm)
 * Since 1 + kq/sm = (sm + kq)/sm = mu/sm:
 * yt = R0/sm * (mu * exp(sm*t) - kq)
 *    = R0/sm * (mu * (exp(sm*t) - 1) + mu - kq)
 *    = R0/sm * (mu * expm1(sm*t) + sm)
 *
 * Case 2: t >= td
 * yt = Rt(t) + Qat(t) + Qct(t)
 * Rt(t) = R0 * exp(sm*t)
 * Qat(t) = (kq*R0/sm) * (exp(sm*t) - 1) - (kq*R0/sm) * (exp(sm*(t-td)) - 1)
 *        = (kq*R0/sm) * (exp(sm*t) - exp(sm*(t-td)))
 * Qct(t) = (kq*R0/(sm+kd)) * (exp((sm+kd)*(t-td)) * exp(kd*td) - exp(kd*td)) * exp(-kd*t)
 *        = (kq*R0/(sm+kd)) * (exp(sm*t - sm*td + kd*t - kd*td + kd*td - kd*t) - exp(kd*td - kd*t))
 *        = (kq*R0/(sm+kd)) * (exp(sm*(t-td)) - exp(-kd*(t-td)))
 *
 * yt = R0*exp(sm*t) + (kq*R0/sm)*(exp(sm*t) - exp(sm*(t-td))) + (kq*R0/(sm+kd))*(exp(sm*(t-td)) - exp(-kd*(t-td)))
 *    = R0*exp(sm*t)*(1 + kq/sm) - (kq*R0/sm)*exp(sm*(t-td)) + (kq*R0/(sm+kd))*exp(sm*(t-td)) - (kq*R0/(sm+kd))*exp(-kd*(t-td))
 *    = R0*exp(sm*t)*(mu/sm) + exp(sm*(t-td)) * R0*kq * (1/(sm+kd) - 1/sm) - (kq*R0/(sm+kd))*exp(-kd*(t-td))
 *    1/(sm+kd) - 1/sm = (sm - (sm+kd)) / (sm*(sm+kd)) = -kd / (sm*(sm+kd))
 * yt = R0*(mu/sm)*exp(sm*t) - R0*kq*kd/(sm*(sm+kd)) * exp(sm*(t-td)) - (kq*R0/(sm+kd))*exp(-kd*(t-td))
 */

real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  if (t < td) {
    return R0 / sm * (mu * expm1(sm * t) + sm);
  } else {
    real term1 = R0 * (mu / sm) * exp(sm * t);
    real term2 = (R0 * kq * kd) / (sm * (sm + kd)) * exp(sm * (t - td));
    real term3 = (R0 * kq) / (sm + kd) * exp(-kd * (t - td));
    real out = term1 - term2 - term3;
    return out > 0 ? out : 1e-9;
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
  real R = y[1];
  real Qa = y[2];
  real Qc = y[3];

  real dR = sm * R;
  real dQa = kq * R - (t < td ? 0 : kq * Rt(t - td, R0, sm));
  real dQc = (t < td ? 0 : kq * Rt(t - td, R0, sm)) - kd * Qc;

  return [dR, dQa, dQc]';
}

real yt_num(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  real out = sum(ode_rk45(dsdt, [R0, 0, 0]', 0, {t}, R0, sm, kq, td, kd)[1]);
  return out > 0 ? out : 1e-9;
}
