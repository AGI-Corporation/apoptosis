/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  if (t < td) {
    return (R0 / sm) * (mu * expm1(sm * t) + sm);
  } else {
    real sm_plus_kd = sm + kd;
    real A = kq * R0 / sm_plus_kd;
    real B = (mu / sm) * exp(sm * td) - (kq / sm);

    real t_minus_td = t - td;
    // Corrected analytic solution for t >= td
    // yt = Rt(t) + Qat(t) + Qct(t)
    // Rt(t) = R0 * exp(sm * t)
    // Qat(t) = (kq * R0 / sm) * (exp(sm * t) - 1) - (kq * R0 / sm) * (exp(sm * (t - td)) - 1)
    //        = (kq * R0 / sm) * (exp(sm * t) - exp(sm * (t - td)))
    // Qct(t) = [kq * R0 / (sm + kd)] * [exp(sm*(t-td)) * exp(kd*td) - exp(kd*td)] * exp(-kd*t)
    //        = [kq * R0 / (sm + kd)] * [exp(sm*(t-td) + kd*td - kd*t) - exp(kd*td - kd*t)]
    //        = [kq * R0 / (sm + kd)] * [exp((sm+kd)*(t-td)) - 1] * exp(kd*td - kd*t) // No, this is getting complicated.
    // Let's use the original formulas but more efficiently.

    real exp_sm_t = exp(sm * t);
    real exp_sm_t_minus_td = exp(sm * t_minus_td);

    real Rt = R0 * exp_sm_t;
    real Qat = (kq * R0 / sm) * (exp_sm_t - exp_sm_t_minus_td);
    real Qct = (kq * R0 / sm_plus_kd) * (exp_sm_t_minus_td - exp(-kd * t_minus_td));

    real out = Rt + Qat + Qct;
    return out > 0 ? out : 1e-9;
  }
}

/* 
   Functions for solving the system numerically - use these to check that the
   analytic solution works.

*/


vector dsdt(real t, vector y, real R0, real sm, real kq, real td, real kd){
  // Rt(t) is cell density of type R.
  // dR/dt = (mu-kq)R = sm*R
  // R(t) = R0 * exp(sm*t)
  // But wait, the original code had:
  // flux = [(sm + kq) * y[1], kq * y[1], t < td ? 0 : kq * Rt(t - td, R0, sm), kd * y[3]]
  // dy[1]/dt = flux[1] - flux[2] = (sm + kq)R - kq*R = sm*R. Correct.
  // dy[2]/dt = flux[2] - flux[3] = kq*R(t) - (t < td ? 0 : kq*R(t-td))
  // dy[3]/dt = flux[3] - flux[4] = (t < td ? 0 : kq*R(t-td)) - kd*Qc

  real Rt_current = R0 * exp(sm * t);
  real Rt_delayed = t < td ? 0 : R0 * exp(sm * (t - td));

  vector[3] dydt;
  dydt[1] = sm * y[1];
  dydt[2] = kq * y[1] - kq * Rt_delayed;
  dydt[3] = kq * Rt_delayed - kd * y[3];

  return dydt;
}

real yt_num(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  if (t == 0) return R0;
  array[1] vector[3] sol = ode_rk45(dsdt, [R0, 0.0, 0.0]', 0.0, {t}, R0, sm, kq, td, kd);
  real out = sum(sol[1]);
  return out > 0 ? out : 1e-9;
}
