import sys
import types
from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import pytest

# Stub out cmdstanpy and its utils before fit_models is imported, since
# cmdstanpy is not available in the test environment and the module-level
# imports in fit_models.py would otherwise fail.
_cmdstanpy_stub = types.ModuleType("cmdstanpy")
_cmdstanpy_stub.CmdStanModel = MagicMock()
_cmdstanpy_utils_stub = types.ModuleType("cmdstanpy.utils")
_cmdstanpy_utils_stub.get_logger = MagicMock(return_value=MagicMock())
_cmdstanpy_utils_stub.jsondump = MagicMock()
sys.modules.setdefault("cmdstanpy", _cmdstanpy_stub)
sys.modules.setdefault("cmdstanpy.utils", _cmdstanpy_utils_stub)

from fit_models import get_infd_kwargs, get_stan_input  # noqa: E402


def make_mock_msmts():
    """Minimal msmts DataFrame with two designs, two clones, two replicates."""
    data = {
        "design": ["Empty", "Empty", "Empty", "A", "A", "A"],
        "clone": ["C1", "C1", "C1", "C2", "C2", "C2"],
        "replicate": ["C1-1st", "C1-1st", "C1-1st", "C2-1st", "C2-1st", "C2-1st"],
        "day": [1.0, 2.0, 3.0, 1.0, 2.0, 3.0],
        "y": [1.0, 0.8, 0.6, 1.1, 0.9, 0.7],
        "design_ab": ["BASE", "BASE", "BASE", "A", "A", "A"],
        "design_abc": ["BASE", "BASE", "BASE", "A", "A", "A"],
        "design_ab_fct": [1, 1, 1, 2, 2, 2],
        "design_abc_fct": [1, 1, 1, 2, 2, 2],
        "design_fct": [1, 1, 1, 2, 2, 2],
        "clone_fct": [1, 1, 1, 2, 2, 2],
        "replicate_fct": [1, 1, 1, 2, 2, 2],
    }
    return pd.DataFrame(data)


MOCK_PRIORS = {
    "prior_mu": (0.0, 1.0),
    "prior_kq": (0.5, 0.3),
}


class TestGetStanInput:
    def setup_method(self):
        self.msmts = make_mock_msmts()
        self.priors = MOCK_PRIORS

    def test_priors_included_in_output(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        for key in self.priors:
            assert key in result

    def test_N_equals_number_of_rows(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        assert result["N"] == len(self.msmts)

    def test_N_is_int(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        assert isinstance(result["N"], int)

    def test_N_test_equals_N(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        assert result["N_test"] == result["N"]

    def test_R_equals_unique_replicate_count(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        assert result["R"] == self.msmts["replicate"].nunique()

    def test_R_is_int(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        assert isinstance(result["R"], int)

    def test_C_equals_unique_clone_count(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        assert result["C"] == self.msmts["clone"].nunique()

    def test_C_is_int(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        assert isinstance(result["C"], int)

    def test_t_matches_day_values(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        np.testing.assert_array_equal(result["t"], self.msmts["day"].values)

    def test_t_test_matches_day_values(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        np.testing.assert_array_equal(result["t_test"], self.msmts["day"].values)

    def test_y_matches_vcd_values(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        np.testing.assert_array_equal(result["y"], self.msmts["y"].values)

    def test_y_test_matches_y(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        np.testing.assert_array_equal(result["y_test"], result["y"])

    def test_replicate_length_equals_N(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        assert len(result["replicate"]) == result["N"]

    def test_clone_array_length_equals_R(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        assert len(result["clone"]) == result["R"]

    def test_design_included_when_not_null(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        assert "design" in result

    def test_D_included_when_not_null(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        assert "D" in result

    def test_D_equals_max_design_fct(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        assert result["D"] == int(self.msmts["design_ab_fct"].max())

    def test_D_is_int(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        assert isinstance(result["D"], int)

    def test_design_excluded_when_null_in_col_name(self):
        result = get_stan_input(self.msmts, self.priors, "design_null")
        assert "design" not in result

    def test_D_excluded_when_null_in_col_name(self):
        result = get_stan_input(self.msmts, self.priors, "design_null")
        assert "D" not in result

    def test_likelihood_is_int(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        assert isinstance(result["likelihood"], int)

    def test_likelihood_value_is_zero_or_one(self):
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        assert result["likelihood"] in (0, 1)

    def test_design_array_length_equals_C(self):
        """design array maps clone to design, so its length should equal C."""
        result = get_stan_input(self.msmts, self.priors, "design_ab")
        assert len(result["design"]) == result["C"]

    def test_design_abc_col_also_supported(self):
        result = get_stan_input(self.msmts, self.priors, "design_abc")
        assert "design" in result
        assert result["D"] == int(self.msmts["design_abc_fct"].max())


class TestGetInfdKwargs:
    def setup_method(self):
        self.msmts = make_mock_msmts()
        self.stan_input = {"y": self.msmts["y"].values}

    def test_returns_dict(self):
        result = get_infd_kwargs(self.msmts, "design_ab", self.stan_input)
        assert isinstance(result, dict)

    def test_log_likelihood_is_llik(self):
        result = get_infd_kwargs(self.msmts, "design_ab", self.stan_input)
        assert result["log_likelihood"] == "llik"

    def test_posterior_predictive_is_yrep(self):
        result = get_infd_kwargs(self.msmts, "design_ab", self.stan_input)
        assert result["posterior_predictive"] == "yrep"

    def test_observed_data_contains_y(self):
        result = get_infd_kwargs(self.msmts, "design_ab", self.stan_input)
        np.testing.assert_array_equal(result["observed_data"]["y"], self.stan_input["y"])

    def test_coords_has_clone(self):
        result = get_infd_kwargs(self.msmts, "design_ab", self.stan_input)
        assert "clone" in result["coords"]

    def test_coords_has_replicate(self):
        result = get_infd_kwargs(self.msmts, "design_ab", self.stan_input)
        assert "replicate" in result["coords"]

    def test_coords_has_cv_effects(self):
        result = get_infd_kwargs(self.msmts, "design_ab", self.stan_input)
        assert "cv_effects" in result["coords"]
        assert result["coords"]["cv_effects"] == ["q", "tau", "d"]

    def test_design_coord_present_when_not_null(self):
        result = get_infd_kwargs(self.msmts, "design_ab", self.stan_input)
        assert "design" in result["coords"]

    def test_design_coord_absent_when_null(self):
        result = get_infd_kwargs(self.msmts, "design_null", self.stan_input)
        assert "design" not in result["coords"]

    def test_base_dims_always_present(self):
        result = get_infd_kwargs(self.msmts, "design_ab", self.stan_input)
        dims = result["dims"]
        for key in ("R0", "cq", "cd", "ct", "sd_cv", "llik"):
            assert key in dims, f"Missing dim key: {key}"

    def test_design_dims_present_when_not_null(self):
        result = get_infd_kwargs(self.msmts, "design_ab", self.stan_input)
        dims = result["dims"]
        for key in ("dq", "dd", "dt", "avg_delay"):
            assert key in dims, f"Missing dim key: {key}"

    def test_design_dims_absent_when_null(self):
        result = get_infd_kwargs(self.msmts, "design_null", self.stan_input)
        dims = result["dims"]
        for key in ("dq", "dd", "dt", "avg_delay"):
            assert key not in dims, f"Unexpected dim key: {key}"

    def test_R0_dim_is_replicate(self):
        result = get_infd_kwargs(self.msmts, "design_ab", self.stan_input)
        assert result["dims"]["R0"] == ["replicate"]

    def test_cq_dim_is_clone(self):
        result = get_infd_kwargs(self.msmts, "design_ab", self.stan_input)
        assert result["dims"]["cq"] == ["clone"]

    def test_llik_dim_is_replicate(self):
        result = get_infd_kwargs(self.msmts, "design_ab", self.stan_input)
        assert result["dims"]["llik"] == ["replicate"]
