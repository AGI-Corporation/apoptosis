import math

def Rt_orig(t, R0, sm):
    return R0 * math.exp(sm * t)

def Qat_orig(t, R0, sm, kq, td):
    U = 0 if t < td else 1
    return kq * R0 / sm * (math.exp(sm * t) - 1) - kq * R0 / sm * (math.exp(sm * (t - td)) - 1) * U

def Qct_orig(t, R0, sm, kq, td, kd):
    U = 0 if t < td else 1
    return U * kq * R0 / (sm + kd) * (math.exp((sm + kd) * (t - td)) * math.exp(kd * td) - math.exp(kd * td)) * math.exp(-kd * t)

def yt_orig(t, R0, mu, kq, td, kd):
    sm = mu - kq
    return Rt_orig(t, R0, sm) + Qat_orig(t, R0, sm, kq, td) + Qct_orig(t, R0, sm, kq, td, kd)

def yt_simp(t, R0, mu, kq, td, kd):
    sm = mu - kq
    if t < td:
        return R0 / sm * (mu * math.exp(sm * t) - kq)
    else:
        return (R0 * mu / sm * math.exp(sm * t)
                - kq * kd * R0 / (sm * (sm + kd)) * math.exp(sm * (t - td))
                - kq * R0 / (sm + kd) * math.exp(-kd * (t - td)))

params = [
    (1.0, 2.5, 0.7, 3.6, 2.8, 0.5), # t < td
    (4.0, 2.5, 0.7, 3.6, 2.8, 0.5), # t > td
    (2.8, 2.5, 0.7, 3.6, 2.8, 0.5), # t == td
]

for p in params:
    orig = yt_orig(*p)
    simp = yt_simp(*p)
    print(f"Params: {p}")
    print(f"Orig: {orig}, Simp: {simp}, Diff: {orig - simp}")
