/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

// U = 0 if t < tau else 1
// R   = Ro*np.exp(sigmu*t)
// Qa  = kq*Ro/sigmu*(np.exp(sigmu*t)-1) - kq*Ro/sigmu*(np.exp(sigmu*(t-tau))-1) * U
// Qc  = ( kq*Ro/(sigmu+kd)*(np.exp((sigmu+kd)*(t-tau))*np.exp(kd*tau) - np.exp(kd*tau))*np.exp(-kd*t) ) * U

real Rt(real t, real R0, real sm) {
  return R0 * exp(sm * t);
}

real Qat(real t, real R0, real sm, real kq, real td) {
  real U = t < td ? 0 : 1;
  return (kq * R0 / sm) * (exp(sm * t) - 1)
         - (kq * R0 / sm) * (exp(sm * (t - td)) - 1) * U;
}

real Qct(real t, real R0, real sm, real kq, real td, real kd) {
  real U = t < td ? 0 : 1;
  return U * (kq * R0 / (sm + kd))
         * (exp((sm + kd) * (t - td)) * exp(kd * td) - exp(kd * td))
         * exp(-kd * t);
}

/**
 * Total cell density yt at time t.
 *
 * Optimized to reduce exp() calls from 7 down to 1 (t < td) or 3 (t >= td).
 */
real yt(real t, real R0, real mu, real kq, real td, real kd) {
  real sm = mu - kq;
  if (t < td) {
    // For t < td: yt = Rt(t) + Qat(t)
    // Rt = R0 * exp(sm * t)
    // Qat = (kq * R0 / sm) * (exp(sm * t) - 1)
    // yt = R0 * exp(sm * t) + (kq * R0 / sm) * exp(sm * t) - (kq * R0 / sm)
    // yt = R0 * exp(sm * t) * (1 + kq / sm) - (kq * R0 / sm)
    // yt = R0 * exp(sm * t) * (sm + kq) / sm - (kq * R0 / sm)
    // yt = R0 * exp(sm * t) * mu / sm - (kq * R0 / sm)
    real val = (R0 / sm) * (mu * exp(sm * t) - kq);
    return val > 0 ? val : 0.00001;
  } else {
    // For t >= td: yt = Rt + Qat + Qct
    // Rt = R0 * exp(sm * t)
    // Qat = (kq * R0 / sm) * (exp(sm * t) - exp(sm * (t - td)))
    // Qct = (kq * R0 / (sm + kd)) * (exp(sm * (t - td)) - exp(kd * (td - t)))
    // Total yt = R0 * exp(sm * t) * (1 + kq / sm)
    //            + exp(sm * (t - td)) * (kq * R0 / (sm + kd) - kq * R0 / sm)
    //            - (kq * R0 / (sm + kd)) * exp(kd * (td - t))
    // yt = (mu * R0 / sm) * exp(sm * t)
    //      - (kq * R0 * kd / (sm * (sm + kd))) * exp(sm * (t - td))
    //      - (kq * R0 / (sm + kd)) * exp(kd * (td - t))
    real term1 = (mu * R0 / sm) * exp(sm * t);
    real term2 = (kq * R0 * kd / (sm * (sm + kd))) * exp(sm * (t - td));
    real term3 = (kq * R0 / (sm + kd)) * exp(kd * (td - t));
    real val = term1 - term2 - term3;
    return val > 0 ? val : 0.00001;
  }
}

/* 
   Functions for solving the system numerically - use these to check that the
   analytic solution works.

*/

vector dsdt(real t, vector y, real R0, real sm, real kq, real td, real kd) {
  vector[4] flux = [(sm + kq) * y[1], kq * y[1],
                    t < td ? 0 : kq * Rt(t - td, R0, sm), kd * y[3]]';
  return [flux[1] - flux[2], flux[2] - flux[3], flux[3] - flux[4]]';
}
real yt_num(real t, real R0, real mu, real kq, real td, real kd) {
  real sm = mu - kq;
  real out = sum(ode_rk45(dsdt, [R0, 0, 0]', 0, {t}, R0, sm, kq, td, kd)[1]);
  return out > 0 ? out : 0.00001;
}
