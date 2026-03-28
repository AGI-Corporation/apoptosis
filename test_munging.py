import pandas as pd
import pytest

from munging import DESIGN_AB, DESIGN_ABC, prepare_data, stan_factorize


def make_raw_data(**overrides):
    """Create a minimal raw DataFrame matching the expected CSV column layout."""
    data = {
        "Clone": ["AB1", "AB1", "B1", "B1", "E1", "E1"],
        "Plasmid": ["Bak-,Bax-", "Bak-,Bax-", "Bax-", "Bax-", "Empty", "Empty"],
        "Run": ["1st", "2nd", "1st", "2nd", "1st", "2nd"],
        "Treatment": ["15ug/mL Puromycin"] * 6,
        "Day": [1, 2, 1, 2, 1, 2],
        "VCD": [1.0, 0.8, 1.1, 0.9, 1.2, 1.0],
    }
    data.update(overrides)
    return pd.DataFrame(data)


class TestStanFactorize:
    def test_basic_returns_integer_codes(self):
        s = pd.Series(["a", "b", "a", "c"])
        result = stan_factorize(s)
        assert set(result.unique()) == {1, 2, 3}

    def test_codes_start_at_one(self):
        s = pd.Series(["x", "y", "z"])
        result = stan_factorize(s)
        assert result.min() == 1
        assert result.max() == 3

    def test_same_value_gets_same_code(self):
        s = pd.Series(["a", "b", "a"])
        result = stan_factorize(s)
        assert result.iloc[0] == result.iloc[2]
        assert result.iloc[0] != result.iloc[1]

    def test_first_param_gets_code_1(self):
        s = pd.Series(["a", "b", "c"])
        result = stan_factorize(s, first="b")
        b_code = result[s == "b"].iloc[0]
        assert b_code == 1

    def test_first_param_shifts_other_codes(self):
        s = pd.Series(["a", "b", "c"])
        result = stan_factorize(s, first="b")
        assert result[s == "b"].iloc[0] == 1
        assert result[s == "a"].iloc[0] != 1
        assert result[s == "c"].iloc[0] != 1

    def test_first_not_in_values_raises_value_error(self):
        s = pd.Series(["a", "b", "c"])
        with pytest.raises(ValueError):
            stan_factorize(s, first="z")

    def test_single_unique_value_gets_code_1(self):
        s = pd.Series(["a", "a", "a"])
        result = stan_factorize(s)
        assert (result == 1).all()

    def test_preserves_index_alignment(self):
        s = pd.Series(["b", "a", "b"], index=[10, 20, 30])
        result = stan_factorize(s, first="a")
        assert result.loc[20] == 1
        assert result.loc[10] == result.loc[30]

    def test_no_first_produces_n_unique_codes(self):
        s = pd.Series(["x", "y", "z", "x", "y"])
        result = stan_factorize(s)
        assert set(result.unique()) == {1, 2, 3}


class TestPrepareData:
    def test_filters_to_requested_treatment(self):
        raw = make_raw_data(
            Treatment=["15ug/mL Puromycin", "other", "15ug/mL Puromycin",
                       "other", "15ug/mL Puromycin", "other"]
        )
        result = prepare_data(raw, treatment="15ug/mL Puromycin")
        assert result["treatment"].eq("15ug/mL Puromycin").all()

    def test_excludes_day_zero_rows(self):
        raw = make_raw_data(Day=[0, 1, 0, 2, 1, 0])
        result = prepare_data(raw, treatment="15ug/mL Puromycin")
        assert result["day"].gt(0).all()

    def test_excludes_none_design_rows(self):
        raw = make_raw_data(
            Plasmid=["None", "Bak-,Bax-", "None", "Bax-", "Empty", "Empty"]
        )
        result = prepare_data(raw, treatment="15ug/mL Puromycin")
        assert "None" not in result["design"].values

    def test_expected_columns_present(self):
        result = prepare_data(make_raw_data(), treatment="15ug/mL Puromycin")
        expected = [
            "design", "clone", "treatment", "replicate", "day", "y",
            "design_ab", "design_abc", "design_fct", "design_ab_fct",
            "design_abc_fct", "clone_fct", "replicate_fct",
            "is_A", "is_B", "is_C", "is_AB", "is_AC", "is_BC", "is_ABC",
        ]
        for col in expected:
            assert col in result.columns, f"Missing column: {col}"

    def test_design_ab_mapping_bak_bax(self):
        result = prepare_data(make_raw_data(), treatment="15ug/mL Puromycin")
        rows = result[result["design"] == "Bak-,Bax-"]
        assert (rows["design_ab"] == "AB").all()

    def test_design_ab_mapping_bax_only(self):
        result = prepare_data(make_raw_data(), treatment="15ug/mL Puromycin")
        rows = result[result["design"] == "Bax-"]
        assert (rows["design_ab"] == "B").all()

    def test_design_ab_mapping_empty(self):
        result = prepare_data(make_raw_data(), treatment="15ug/mL Puromycin")
        rows = result[result["design"] == "Empty"]
        assert (rows["design_ab"] == "BASE").all()

    def test_design_abc_mapping_bok_only(self):
        raw = make_raw_data(
            Plasmid=["Bok-", "Bok-", "Bax-", "Bax-", "Empty", "Empty"]
        )
        result = prepare_data(raw, treatment="15ug/mL Puromycin")
        rows = result[result["design"] == "Bok-"]
        assert (rows["design_abc"] == "C").all()

    def test_design_abc_mapping_all_three(self):
        raw = make_raw_data(
            Plasmid=["Bak-,Bax-,Bok-", "Bak-,Bax-,Bok-", "Bax-", "Bax-",
                     "Empty", "Empty"]
        )
        result = prepare_data(raw, treatment="15ug/mL Puromycin")
        rows = result[result["design"] == "Bak-,Bax-,Bok-"]
        assert (rows["design_abc"] == "ABC").all()

    def test_is_A_flag_set_for_bak_designs(self):
        result = prepare_data(make_raw_data(), treatment="15ug/mL Puromycin")
        bak_rows = result[result["design"].str.contains("Bak")]
        assert bak_rows["is_A"].all()

    def test_is_A_flag_not_set_for_non_bak(self):
        result = prepare_data(make_raw_data(), treatment="15ug/mL Puromycin")
        non_bak = result[~result["design"].str.contains("Bak")]
        assert not non_bak["is_A"].any()

    def test_is_B_flag_set_for_bax_designs(self):
        result = prepare_data(make_raw_data(), treatment="15ug/mL Puromycin")
        bax_rows = result[result["design"].str.contains("Bax")]
        assert bax_rows["is_B"].all()

    def test_is_C_flag_set_for_bok_designs(self):
        raw = make_raw_data(
            Plasmid=["Bok-", "Bok-", "Bax-", "Bax-", "Empty", "Empty"]
        )
        result = prepare_data(raw, treatment="15ug/mL Puromycin")
        bok_rows = result[result["design"].str.contains("Bok")]
        assert bok_rows["is_C"].all()

    def test_is_AB_requires_both_bak_and_bax(self):
        result = prepare_data(make_raw_data(), treatment="15ug/mL Puromycin")
        bak_bax = result[result["design"] == "Bak-,Bax-"]
        assert bak_bax["is_AB"].all()
        bax_only = result[result["design"] == "Bax-"]
        assert not bax_only["is_AB"].any()

    def test_is_ABC_requires_all_three(self):
        raw = make_raw_data(
            Plasmid=["Bak-,Bax-,Bok-", "Bak-,Bax-,Bok-", "Bak-,Bax-",
                     "Bak-,Bax-", "Empty", "Empty"]
        )
        result = prepare_data(raw, treatment="15ug/mL Puromycin")
        abc_rows = result[result["design"] == "Bak-,Bax-,Bok-"]
        assert abc_rows["is_ABC"].all()
        ab_rows = result[result["design"] == "Bak-,Bax-"]
        assert not ab_rows["is_ABC"].any()

    def test_replicate_combines_clone_and_run(self):
        result = prepare_data(make_raw_data(), treatment="15ug/mL Puromycin")
        assert result["replicate"].str.contains("-").all()
        assert result["replicate"].str.contains("AB1").any()

    def test_factorize_cols_start_at_one(self):
        result = prepare_data(make_raw_data(), treatment="15ug/mL Puromycin")
        for col in ["clone_fct", "replicate_fct", "design_fct",
                    "design_ab_fct", "design_abc_fct"]:
            assert result[col].min() == 1, f"{col} minimum should be 1"

    def test_design_fct_first_is_empty(self):
        """'Empty' design should map to factor code 1."""
        result = prepare_data(make_raw_data(), treatment="15ug/mL Puromycin")
        empty_code = result.loc[result["design"] == "Empty", "design_fct"].iloc[0]
        assert empty_code == 1

    def test_design_ab_fct_first_is_base(self):
        """'BASE' design_ab should map to factor code 1."""
        result = prepare_data(make_raw_data(), treatment="15ug/mL Puromycin")
        base_code = result.loc[result["design_ab"] == "BASE", "design_ab_fct"].iloc[0]
        assert base_code == 1

    def test_raises_when_no_matching_treatment(self):
        """prepare_data raises ValueError when given a treatment not present in the data."""
        with pytest.raises(ValueError):
            prepare_data(make_raw_data(), treatment="nonexistent")

    def test_y_column_is_vcd(self):
        raw = make_raw_data(VCD=[2.5, 3.0, 1.5, 2.0, 4.0, 3.5])
        result = prepare_data(raw, treatment="15ug/mL Puromycin")
        assert set(result["y"]) == {2.5, 3.0, 1.5, 2.0, 4.0, 3.5}

    def test_day_column_is_float(self):
        result = prepare_data(make_raw_data(), treatment="15ug/mL Puromycin")
        assert result["day"].dtype == float


class TestDesignMappings:
    def test_design_ab_bak_bax_is_ab(self):
        assert DESIGN_AB["Bak-,Bax-"] == "AB"

    def test_design_ab_empty_is_base(self):
        assert DESIGN_AB["Empty"] == "BASE"

    def test_design_ab_none_is_base(self):
        assert DESIGN_AB["None"] == "BASE"

    def test_design_ab_bok_only_is_base(self):
        assert DESIGN_AB["Bok-"] == "BASE"

    def test_design_ab_bak_only_is_a(self):
        assert DESIGN_AB["Bak-"] == "A"

    def test_design_ab_bax_only_is_b(self):
        assert DESIGN_AB["Bax-"] == "B"

    def test_design_ab_bak_bok_is_a(self):
        assert DESIGN_AB["Bak-,Bok-"] == "A"

    def test_design_ab_bax_bok_is_b(self):
        assert DESIGN_AB["Bax-,Bok-"] == "B"

    def test_design_ab_all_three_is_ab(self):
        assert DESIGN_AB["Bak-,Bax-,Bok-"] == "AB"

    def test_design_abc_all_three_is_abc(self):
        assert DESIGN_ABC["Bak-,Bax-,Bok-"] == "ABC"

    def test_design_abc_bok_only_is_c(self):
        assert DESIGN_ABC["Bok-"] == "C"

    def test_design_abc_bak_bok_is_ac(self):
        assert DESIGN_ABC["Bak-,Bok-"] == "AC"

    def test_design_abc_bax_bok_is_bc(self):
        assert DESIGN_ABC["Bax-,Bok-"] == "BC"

    def test_design_abc_empty_is_base(self):
        assert DESIGN_ABC["Empty"] == "BASE"

    def test_design_abc_bak_only_is_a(self):
        assert DESIGN_ABC["Bak-"] == "A"

    def test_design_abc_bax_only_is_b(self):
        assert DESIGN_ABC["Bax-"] == "B"
