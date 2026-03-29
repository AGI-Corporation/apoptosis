/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  real val;
  if (t < td){
    // Simplified analytic solution for t < td (1 exp call)
    val = (R0 / sm) * (mu * exp(sm * t) - kq);
  } else {
    // Simplified analytic solution for t >= td (max 3 exp calls)
    real sm_kd = sm + kd;
    real term1 = (mu * R0 / sm) * exp(sm * t);
    real term2 = (kq * kd * R0 / (sm * sm_kd)) * exp(sm * (t - td));
    real term3 = (kq * R0 / sm_kd) * exp(-kd * (t - td));
    val = term1 - term2 - term3;
  }
  return fmax(val, 1e-9); // Floor for numerical stability
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
