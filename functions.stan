/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

real yt(real t, real R0, real mu, real kq, real td, real kd) {
  real sm = mu - kq;
  real out;

  if (t < td) {
    // Simplified: R0 * exp(sm * t) + kq * R0 / sm * (exp(sm * t) - 1)
    // = R0 * (exp(sm * t) + kq / sm * exp(sm * t) - kq / sm)
    // = R0 * ( (1 + kq/sm) * exp(sm * t) - kq / sm )
    // = R0 * ( (sm + kq)/sm * exp(sm * t) - kq / sm )
    // = R0 * ( mu / sm * exp(sm * t) - kq / sm )
    out = R0 * (mu / sm * exp(sm * t) - kq / sm);
  } else {
    real sm_kd = sm + kd;
    // Simplified:
    // Rt = R0 * exp(sm * t)
    // Qat = kq * R0 / sm * (exp(sm * t) - 1) - kq * R0 / sm * (exp(sm * (t - td)) - 1)
    // Qct = kq * R0 / sm_kd * (exp(sm * (t - td)) - exp(-kd * (t - td)))
    // yt = R0 * mu / sm * exp(sm * t) - kq * R0 / sm - kq * R0 / sm * exp(sm * (t - td)) + kq * R0 / sm + Qct
    // yt = R0 * mu / sm * exp(sm * t) - kq * R0 / sm * exp(sm * (t - td)) + kq * R0 / sm_kd * exp(sm * (t - td)) - kq * R0 / sm_kd * exp(-kd * (t - td))
    // yt = R0 * mu / sm * exp(sm * t) - kq * R0 * (1/sm - 1/sm_kd) * exp(sm * (t - td)) - kq * R0 / sm_kd * exp(-kd * (t - td))
    // 1/sm - 1/sm_kd = (sm + kd - sm) / (sm * sm_kd) = kd / (sm * sm_kd)
    out = R0 * mu / sm * exp(sm * t)
          - kq * kd * R0 / (sm * sm_kd) * exp(sm * (t - td))
          - kq * R0 / sm_kd * exp(-kd * (t - td));
  }

  return out > 0 ? out : 0.00001;
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
  return out > 0 ? out : 0.00001;
}
