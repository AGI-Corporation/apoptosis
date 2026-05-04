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
 * Algebraically simplified analytic solution for cell density.
 * Reduces the number of exp() calls from up to 6 down to 3,
 * and uses expm1() for better numerical stability when t < td.
 */
real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  if (t < td) {
    // simplified from R(t) + Qa(t)
    // R0 * exp(sm*t) + kq*R0/sm * (exp(sm*t) - 1)
    // = R0/sm * (sm*exp(sm*t) + kq*exp(sm*t) - kq)
    // = R0/sm * (mu*exp(sm*t) - kq)
    // = R0/sm * (mu*(exp(sm*t)-1) + mu - kq)
    // = R0/sm * (mu*expm1(sm*t) + sm)
    return (R0 / sm) * (mu * expm1(sm * t) + sm);
  } else {
    // simplified from R(t) + Qa(t) + Qc(t)
    // = R0 * mu/sm * exp(sm*t) - (kq*kd*R0)/(sm*(sm+kd)) * exp(sm*(t-td)) - (kq*R0)/(sm+kd) * exp(-kd*(t-td))
    real sm_kd = sm + kd;
    return R0 * ( (mu / sm) * exp(sm * t)
                - (kq * kd / (sm * sm_kd)) * exp(sm * (t - td))
                - (kq / sm_kd) * exp(-kd * (t - td)) );
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
  return out > 0 ? out : 1e-9;
}
