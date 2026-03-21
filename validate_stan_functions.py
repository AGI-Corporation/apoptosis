from cmdstanpy import CmdStanModel
import pandas as pd

STAN_FILE = "validation_model.stan"
DATA = {
    "N": 5,
    "t": [1.0, 2.0, 3.0, 4.0, 5.0],
    "R0": 2.5,
    "sm": 0.5, # Positive growth rate for easier visualization/debugging
    "kq": 0.2,
    "td": 2.5,
    "kd": 0.1,
}
SAMPLE_CONFIG = dict(
    fixed_param=True,
    iter_sampling=1,
    iter_warmup=0,
    adapt_engaged=False,
)

def main():
   model = CmdStanModel(stan_file=STAN_FILE)
   mcmc = model.sample(DATA, **SAMPLE_CONFIG)
   print("input data:\n", DATA)
   draws = mcmc.draws_pd()
   y_n = [draws[f"y_n[{i+1}]"].iloc[0] for i in range(DATA["N"])]
   y_a = [draws[f"y_a[{i+1}]"].iloc[0] for i in range(DATA["N"])]
   df = pd.DataFrame({"t": DATA["t"], "numerical": y_n, "analytic": y_a})
   df["diff"] = (df["numerical"] - df["analytic"]).abs()
   print(df)
   if (df["diff"] < 1e-5).all():
       print("\nSUCCESS: Analytic solution matches numerical solution.")
   else:
       print("\nFAILURE: Significant difference between solutions.")

if __name__ == "__main__":
    main()
