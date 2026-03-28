import numpy as np
import pytest
from scipy.stats import lognorm, norm

from util import (
    get_99_pct_params_ln,
    get_99_pct_params_n,
    get_lognormal_params_from_qs,
    get_normal_params_from_qs,
)


class TestGetLognormalParamsFromQs:
    def test_basic_round_trip(self):
        """Parameters should reproduce the specified quantiles."""
        x1, x2, p1, p2 = 1.0, 2.0, 0.1, 0.9
        mu, sigma = get_lognormal_params_from_qs(x1, x2, p1, p2)
        assert lognorm.ppf(p1, s=sigma, scale=np.exp(mu)) == pytest.approx(x1, rel=1e-6)
        assert lognorm.ppf(p2, s=sigma, scale=np.exp(mu)) == pytest.approx(x2, rel=1e-6)

    def test_sigma_is_positive(self):
        _, sigma = get_lognormal_params_from_qs(0.5, 2.0, 0.05, 0.95)
        assert sigma > 0

    def test_tight_range_gives_small_sigma(self):
        """A narrow quantile interval should yield a small sigma."""
        _, sigma = get_lognormal_params_from_qs(0.99, 1.01, 0.01, 0.99)
        assert sigma < 0.1

    def test_wide_range_gives_large_sigma(self):
        _, sigma_wide = get_lognormal_params_from_qs(0.01, 100.0, 0.01, 0.99)
        _, sigma_narrow = get_lognormal_params_from_qs(0.9, 1.1, 0.01, 0.99)
        assert sigma_wide > sigma_narrow

    def test_returns_two_floats(self):
        result = get_lognormal_params_from_qs(1.0, 2.0, 0.25, 0.75)
        assert len(result) == 2
        assert all(isinstance(v, float) for v in result)

    def test_known_case_median(self):
        """When p1=0.5 and x1=x2, mu should equal log(x1)."""
        x = 3.0
        mu, _ = get_lognormal_params_from_qs(x, x * np.e, 0.5, norm.cdf(1))
        assert mu == pytest.approx(np.log(x), rel=1e-6)

    def test_priors_from_fit_models(self):
        """Values used in fit_models.py should produce valid parameters."""
        cases = [
            (0.65, 0.73), (1, 5), (0.4, 7.5), (0.05, 2.5), (2, 3), (0.05, 0.13)
        ]
        for x1, x2 in cases:
            mu, sigma = get_lognormal_params_from_qs(x1, x2, 0.01, 0.99)
            assert sigma > 0
            assert lognorm.ppf(0.01, s=sigma, scale=np.exp(mu)) == pytest.approx(x1, rel=1e-5)
            assert lognorm.ppf(0.99, s=sigma, scale=np.exp(mu)) == pytest.approx(x2, rel=1e-5)


class TestGetNormalParamsFromQs:
    def test_basic_round_trip(self):
        """Parameters should reproduce the specified quantiles."""
        x1, x2, p1, p2 = -1.0, 1.0, 0.1, 0.9
        mu, sigma = get_normal_params_from_qs(x1, x2, p1, p2)
        assert norm.ppf(p1, loc=mu, scale=sigma) == pytest.approx(x1, rel=1e-6)
        assert norm.ppf(p2, loc=mu, scale=sigma) == pytest.approx(x2, rel=1e-6)

    def test_symmetric_quantiles_give_zero_mu(self):
        """Symmetric quantile interval around zero should give mu = 0."""
        mu, _ = get_normal_params_from_qs(-1.0, 1.0, 0.05, 0.95)
        assert mu == pytest.approx(0.0, abs=1e-10)

    def test_sigma_is_positive(self):
        _, sigma = get_normal_params_from_qs(0.0, 1.0, 0.01, 0.99)
        assert sigma > 0

    def test_recovers_standard_normal(self):
        """Recover mu=0, sigma=1 from the ±1-sigma quantiles."""
        p1 = norm.cdf(-1.0)
        p2 = norm.cdf(1.0)
        mu, sigma = get_normal_params_from_qs(-1.0, 1.0, p1, p2)
        assert mu == pytest.approx(0.0, abs=1e-10)
        assert sigma == pytest.approx(1.0, rel=1e-6)

    def test_recovers_shifted_distribution(self):
        """Recover known mu and sigma for a shifted normal distribution."""
        mu_true, sigma_true = 5.0, 2.0
        x1 = norm.ppf(0.25, loc=mu_true, scale=sigma_true)
        x2 = norm.ppf(0.75, loc=mu_true, scale=sigma_true)
        mu, sigma = get_normal_params_from_qs(x1, x2, 0.25, 0.75)
        assert mu == pytest.approx(mu_true, rel=1e-6)
        assert sigma == pytest.approx(sigma_true, rel=1e-6)

    def test_returns_two_floats(self):
        result = get_normal_params_from_qs(-1.0, 1.0, 0.1, 0.9)
        assert len(result) == 2
        assert all(isinstance(v, float) for v in result)

    def test_mu_shifts_with_interval(self):
        """Shifting the quantile interval should shift mu accordingly."""
        mu_low, _ = get_normal_params_from_qs(0.0, 1.0, 0.25, 0.75)
        mu_high, _ = get_normal_params_from_qs(10.0, 11.0, 0.25, 0.75)
        assert mu_high > mu_low


class TestGet99PctParamsLn:
    def test_delegates_to_lognormal(self):
        """Should equal get_lognormal_params_from_qs with p1=0.01, p2=0.99."""
        x1, x2 = 1.0, 5.0
        mu1, s1 = get_99_pct_params_ln(x1, x2)
        mu2, s2 = get_lognormal_params_from_qs(x1, x2, 0.01, 0.99)
        assert mu1 == pytest.approx(mu2)
        assert s1 == pytest.approx(s2)

    def test_round_trip_1st_99th_percentiles(self):
        """1st and 99th percentiles should match the input values."""
        x1, x2 = 0.65, 0.73
        mu, sigma = get_99_pct_params_ln(x1, x2)
        assert lognorm.ppf(0.01, s=sigma, scale=np.exp(mu)) == pytest.approx(x1, rel=1e-6)
        assert lognorm.ppf(0.99, s=sigma, scale=np.exp(mu)) == pytest.approx(x2, rel=1e-6)

    def test_sigma_is_positive(self):
        _, sigma = get_99_pct_params_ln(2.0, 3.0)
        assert sigma > 0


class TestGet99PctParamsN:
    def test_delegates_to_normal(self):
        """Should equal get_normal_params_from_qs with p1=0.01, p2=0.99."""
        x1, x2 = -2.0, 2.0
        mu1, s1 = get_99_pct_params_n(x1, x2)
        mu2, s2 = get_normal_params_from_qs(x1, x2, 0.01, 0.99)
        assert mu1 == pytest.approx(mu2)
        assert s1 == pytest.approx(s2)

    def test_round_trip_1st_99th_percentiles(self):
        """1st and 99th percentiles should match the input values."""
        x1, x2 = -1.0, 1.0
        mu, sigma = get_99_pct_params_n(x1, x2)
        assert norm.ppf(0.01, loc=mu, scale=sigma) == pytest.approx(x1, rel=1e-6)
        assert norm.ppf(0.99, loc=mu, scale=sigma) == pytest.approx(x2, rel=1e-6)

    def test_sigma_is_positive(self):
        _, sigma = get_99_pct_params_n(-1.0, 1.0)
        assert sigma > 0
