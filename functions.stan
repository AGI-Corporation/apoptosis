/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  if (t < td) {
    // return Rt(t, R0, sm) + Qat(t, R0, sm, kq, td);
    // simplified: (R0/sm) * (mu * exp(sm*t) - kq)
    return (R0 / sm) * (mu * exp(sm * t) - kq);
  } else {
    // return Rt(t, R0, sm) + Qat(t, R0, sm, kq, td) + Qct(t, R0, sm, kq, td, kd);
    // simplified: (R0*mu/sm)*exp(sm*t) - (kq*kd*R0/(sm*(sm+kd)))*exp(sm*(t-td)) - (kq*R0/(sm+kd))*exp(-kd*(t-td))
    real term1 = (R0 * mu / sm) * exp(sm * t);
    real term2 = (kq * kd * R0 / (sm * (sm + kd))) * exp(sm * (t - td));
    real term3 = (kq * R0 / (sm + kd)) * exp(-kd * (t - td));
    real res = term1 - term2 - term3;
    // Add a small floor to ensure numerical stability as in original yt_num
    return res > 1e-9 ? res : 1e-9;
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
  return out > 0 ? out : 0.00001;
}
