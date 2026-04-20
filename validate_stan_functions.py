from cmdstanpy import CmdStanModel

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
    adapt_engaged=False,
)

def main():
   model = CmdStanModel(stan_file=STAN_FILE)
   mcmc = model.sample(DATA, **SAMPLE_CONFIG)
   print("input data:\n", DATA)
   df = mcmc.draws_pd()
   y_n = [col for col in df.columns if col.startswith('y_n')]
   y_a = [col for col in df.columns if col.startswith('y_a')]
   print("Numerical:\n", df[y_n].values)
   print("Analytic:\n", df[y_a].values)


if __name__ == "__main__":
    main()
