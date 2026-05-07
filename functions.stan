/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  if (t < td) {
    // Simplified: R(t) + Qa(t)
    // yt = R0 * (exp(sm * t) + kq * expm1(sm * t) / sm)
    return fmax(1e-9, R0 * (exp(sm * t) + kq * (expm1(sm * t) / sm)));
  } else {
    // Simplified: R(t) + Qa(t) + Qc(t)
    // yt = R0 * [ (mu / sm) * exp(sm * t)
    //            - (kq * kd / (sm * (sm + kd))) * exp(sm * (t - td))
    //            - (kq / (sm + kd)) * exp(-kd * (t - td)) ]

    real exp_sm_t = exp(sm * t);
    real exp_sm_t_td = exp(sm * (t - td));
    real exp_kd_t_td = exp(-kd * (t - td));

    real val = R0 * ( (mu / sm) * exp_sm_t
                 - (kq * kd / (sm * (sm + kd))) * exp_sm_t_td
                 - (kq / (sm + kd)) * exp_kd_t_td );

    return fmax(1e-9, val);
  }
}

/* 
   Functions for solving the system numerically - use these to check that the
   analytic solution works.

*/


vector dsdt(real t, vector y, real R0, real sm, real kq, real td, real kd){
  vector[3] out;
  out[1] = sm * y[1];
  out[2] = kq * y[1] - (t < td ? 0 : kq * R0 * exp(sm * (t - td)));
  out[3] = (t < td ? 0 : kq * R0 * exp(sm * (t - td))) - kd * y[3];
  return out;
}

real yt_num(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  real out = sum(ode_rk45(dsdt, [R0, 0, 0]', 0, {t}, R0, sm, kq, td, kd)[1]);
  return out > 0 ? out : 1e-9;
}
