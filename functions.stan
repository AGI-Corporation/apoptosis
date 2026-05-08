/*
  Functions for analytically solving the ODE system in order to find the total
  cell density at a given time, given initial density R0 and parameters p where:

*/

real yt(real t, real R0, real mu, real kq, real td, real kd){
  real sm = mu - kq;
  if (t < td) {
    // R0 * exp(sm * t) + kq * R0 / sm * (exp(sm * t) - 1)
    // = R0 * exp(sm * t) + kq * R0 / sm * exp(sm * t) - kq * R0 / sm
    // = R0 * exp(sm * t) * (1 + kq / sm) - kq * R0 / sm
    // = R0 * exp(sm * t) * (sm + kq) / sm - kq * R0 / sm
    // = R0 / sm * (mu * exp(sm * t) - kq)

    // Using expm1 for numerical stability when sm*t is small
    return R0 / sm * (mu * expm1(sm * t) + sm);
  } else {
    // Rt = R0 * exp(sm * t)
    // Qat = kq * R0 / sm * (exp(sm * t) - 1) - kq * R0 / sm * (exp(sm * (t - td)) - 1)
    //     = kq * R0 / sm * (exp(sm * t) - exp(sm * (t - td)))
    // Qct = kq * R0 / (sm + kd) * (exp((sm + kd) * (t - td)) * exp(kd * td) - exp(kd * td)) * exp(-kd * t)
    //     = kq * R0 / (sm + kd) * (exp(sm * (t - td) + kd * (t - td) + kd * td - kd * t) - exp(kd * td - kd * t))
    //     = kq * R0 / (sm + kd) * (exp(sm * (t - td)) - exp(-kd * (t - td)))

    // yt = Rt + Qat + Qct
    // yt = R0 * exp(sm * t) + kq * R0 / sm * (exp(sm * t) - exp(sm * (t - td))) + kq * R0 / (sm + kd) * (exp(sm * (t - td)) - exp(-kd * (t - td)))
    // yt = R0 * exp(sm * t) * (1 + kq / sm) - kq * R0 / sm * exp(sm * (t - td)) + kq * R0 / (sm + kd) * exp(sm * (t - td)) - kq * R0 / (sm + kd) * exp(-kd * (t - td))
    // yt = R0 * exp(sm * t) * (mu / sm) + R0 * kq * exp(sm * (t - td)) * (1 / (sm + kd) - 1 / sm) - R0 * kq / (sm + kd) * exp(-kd * (t - td))
    // 1/(sm+kd) - 1/sm = (sm - (sm + kd)) / (sm * (sm + kd)) = -kd / (sm * (sm + kd))

    // yt = R0 * (mu / sm * exp(sm * t) - (kq * kd) / (sm * (sm + kd)) * exp(sm * (t - td)) - kq / (sm + kd) * exp(-kd * (t - td)))

    real term1 = mu / sm * exp(sm * t);
    real common_q = kq / (sm + kd);
    real term2 = common_q * kd / sm * exp(sm * (t - td));
    real term3 = common_q * exp(-kd * (t - td));

    real out = R0 * (term1 - term2 - term3);
    return out > 1e-9 ? out : 1e-9;
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
