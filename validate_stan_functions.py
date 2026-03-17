from cmdstanpy import CmdStanModel
import os

STAN_FILE = "validation_model.stan"
DATA = {
    "N": 5,
    "t": [1, 2, 3, 4, 5],
    "R0": 2.5,
    "sm": -2.956,
    "kq": 3.686,
    "td": 2.839,
    "kd": 0.499,
}
SAMPLE_CONFIG = dict(
    fixed_param=True,
    iter_sampling=1,
    iter_warmup=0,
    stanc_options={"include-paths": ["."]}
)

def main():
   model = CmdStanModel(stan_file=STAN_FILE)
   mcmc = model.sample(DATA, **SAMPLE_CONFIG)
   print("input data:\n", DATA)
   pd_draws = mcmc.draws_pd().T
   print("results:\n", pd_draws)

   # Extract y_n and y_a
   y_n = [pd_draws.loc[f"y_n[{i+1}]", 0] for i in range(DATA["N"])]
   y_a = [pd_draws.loc[f"y_a[{i+1}]", 0] for i in range(DATA["N"])]

   print("\nVerification:")
   for i in range(DATA["N"]):
       diff = abs(y_n[i] - y_a[i])
       print(f"t={DATA['t'][i]}: numerical={y_n[i]:.6f}, analytic={y_a[i]:.6f}, diff={diff:.6e}")

if __name__ == "__main__":
    # Ensure CMDSTAN environment variable is set for the script if not already
    if "CMDSTAN" not in os.environ:
        os.environ["CMDSTAN"] = os.path.expanduser("~/.cmdstan/cmdstan-2.26.1")
    main()
