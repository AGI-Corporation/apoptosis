/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

real yt(real t, real R0, real mu, real kq, real td, real kd) {
  real sm = mu - kq;
  if (t < td) {
    // Rt + Qat = R0 * exp(sm*t) + kq*R0/sm * (exp(sm*t) - 1)
    // = R0/sm * (sm*exp(sm*t) + kq*exp(sm*t) - kq)
    // = R0/sm * (mu*exp(sm*t) - kq)
    // = R0/sm * (mu*expm1(sm*t) + sm)
    return (R0 / sm) * (mu * expm1(sm * t) + sm);
  } else {
    // Simplified analytic solution for t >= td
    // Reduces number of exp() calls and arithmetic operations
    real term1 = (mu / sm) * exp(sm * t);
    real term2 = (kq * kd / (sm * (sm + kd))) * exp(sm * (t - td));
    real term3 = (kq / (sm + kd)) * exp(-kd * (t - td));
    return R0 * (term1 - term2 - term3);
  }
}

/* 
   Functions for solving the system numerically - use these to check that the
   analytic solution works.

*/

vector dsdt(real t, vector y, real R0, real sm, real kq, real td, real kd){
  real Rt_val = R0 * exp(sm * t);
  vector[4] flux = [(sm + kq) * y[1],
                    kq * y[1],
                    t < td ? 0 : kq * R0 * exp(sm * (t - td)),
                    kd * y[3]]';
  return [flux[1]-flux[2], flux[2]-flux[3], flux[3]-flux[4]]';
}

real yt_num(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  real out = sum(ode_rk45(dsdt, [R0, 0, 0]', 0, {t}, R0, sm, kq, td, kd)[1]);
  return out > 0 ? out : 0.00001;
}
