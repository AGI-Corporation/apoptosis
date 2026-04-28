/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

real Rt(real t, real R0, real sm){
  return R0 * exp(sm * t);
}

real Qat(real t, real R0, real sm, real kq, real td){
  if (t < td) return 0;
  return kq * R0 / sm * (exp(sm * t) - 1)
    - kq * R0 / sm * (exp(sm * (t - td)) - 1);
}

real Qct(real t, real R0, real sm, real kq, real td, real kd){
  if (t < td) return 0;
  return kq * R0 / (sm + kd)
    * (exp((sm + kd) * (t - td)) * exp(kd * td) - exp(kd * td))
    * exp(-kd * t);
}

real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  real out;
  if (t < td){
    out = R0 / sm * (mu * expm1(sm * t) + sm);
  } else {
    out = R0 * mu / sm * exp(sm * t)
          - R0 * kq * kd / (sm * (sm + kd)) * exp(sm * (t - td))
          - R0 * kq / (sm + kd) * exp(-kd * (t - td));
  }
  return out > 1e-9 ? out : 1e-9;
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
