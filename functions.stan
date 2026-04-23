/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:
*/

/**
 * Total cell density at time t.
 *
 * @param t Time
 * @param R0 Initial cell density
 * @param mu Growth rate
 * @param kq Death rate constant
 * @param td Delay time
 * @param kd Rate constant for transition to dead cells
 * @return Total cell density
 */
real yt(real t, real R0, real mu, real kq, real td, real kd) {
  real sm = mu - kq;
  real val;
  if (t < td) {
    // Algebraically simplified analytic solution for t < td
    val = (R0 / sm) * (mu * expm1(sm * t) + sm);
  } else {
    // Algebraically simplified analytic solution for t >= td
    // Reduces number of exp() calls from 6 to 3 compared to the modular version.
    val = R0 * ( (mu / sm) * exp(sm * t)
                 - (kq * kd) / (sm * (sm + kd)) * exp(sm * (t - td))
                 - (kq / (sm + kd)) * exp(-kd * (t - td)) );
  }
  // Small floor to ensure numerical stability and prevent log(0)
  return fmax(val, 1e-9);
}

/*
   Functions for solving the system numerically - use these to check that the
   analytic solution works.
*/

vector dsdt(real t, vector y, real R0, real sm, real kq, real td, real kd) {
  vector[4] flux;
  flux[1] = (sm + kq) * y[1];
  flux[2] = kq * y[1];
  flux[3] = t < td ? 0 : kq * R0 * exp(sm * (t - td));
  flux[4] = kd * y[3];
  return [flux[1] - flux[2], flux[2] - flux[3], flux[3] - flux[4]]';
}

real yt_num(real t, real R0, real mu, real kq, real td, real kd) {
  real sm = mu - kq;
  real out = sum(ode_rk45(dsdt, [R0, 0, 0]', 0, {t}, R0, sm, kq, td, kd)[1]);
  return fmax(out, 1e-9);
}
