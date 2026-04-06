/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

// U = 0 if t < tau else 1
        // R   = Ro*math.exp(sm * t)
        // Qa  = kq*Ro/sm*(math.exp(sm * t)-1) - kq*Ro/sm*(math.exp(sm * (t-td)) - 1) * U
        // Qc  = ( kq*Ro/(sm + kd) * (math.exp((sm + kd) * (t - td)) * math.exp(kd * td) - math.exp(kd * td)) * math.exp(-kd * t) ) * U

real yt(real t, real R0, real mu, real kq, real td, real kd) {
  real sm = mu - kq;
  if (t < td) {
    return R0 * exp(sm * t) + (kq * R0 / sm) * (exp(sm * t) - 1);
  } else {
    real term1 = R0 * exp(sm * t);
    real term2 = (kq * R0 / sm) * (exp(sm * t) - exp(sm * (t - td)));
    real term3 = (kq * R0 / (sm + kd)) * (exp(sm * (t - td)) - exp(-kd * (t - td)));
    return term1 + term2 + term3 + 1e-9;
  }
}

/* 
   Functions for solving the system numerically - use these to check that the
   analytic solution works.

*/


vector dsdt(real t, vector y, real R0, real sm, real kq, real td, real kd){
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
