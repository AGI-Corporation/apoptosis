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
  if (t < td) {
    return (kq * R0 / sm) * expm1(sm * t);
  } else {
    return (kq * R0 / sm) * (exp(sm * t) - exp(sm * (t - td)));
  }
}

real Qct(real t, real R0, real sm, real kq, real td, real kd){
  if (t < td) {
    return 0;
  } else {
    return (kq * R0 / (sm + kd)) * (exp(sm * (t - td)) - exp(-kd * (t - td)));
  }
}

real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  real val = Rt(t, R0, sm) + Qat(t, R0, sm, kq, td) + Qct(t, R0, sm, kq, td, kd);
  return val > 1e-9 ? val : 1e-9;
}

/* 
   Functions for solving the system numerically - use these to check that the
   analytic solution works.

*/


vector dsdt(real t, vector y, real R0, real sm, real kq, real td, real kd){
  vector[4] flux;
  flux[1] = (sm + kq) * y[1];
  flux[2] = kq * y[1];
  flux[3] = t < td ? 0 : kq * Rt(t - td, R0, sm);
  flux[4] = kd * y[3];
  return [flux[1]-flux[2], flux[2]-flux[3], flux[3]-flux[4]]';
}

real yt_num(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  real out = sum(ode_rk45(dsdt, [R0, 0, 0]', 0, {t}, R0, sm, kq, td, kd)[1]);
  return out > 0 ? out : 1e-9;
}
