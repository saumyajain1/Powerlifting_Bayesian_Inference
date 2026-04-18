# Powerlifting Bayesian Inference

This project studies how raw powerlifting performance scales with bodyweight, age, and sex using the 2025 OpenPowerlifting raw SBD results. The analysis is framed as a Bayesian workflow with motivated priors, prior predictive checks, posterior predictive checks, model comparison, and a direct comparison between Hamiltonian Monte Carlo and variational inference. The modeling strategy starts with a Gaussian generalized linear model (equivalently, Bayesian linear regression with identity link) and then extends beyond standard GLMs to a nonlinear Student-t model. The main substantive question is whether the relationship between bodyweight and total strength shows diminishing returns.

The project compares a Gaussian linear baseline with a richer nonlinear Student-t model. The baseline is validated against a semi-analytic reference posterior, while the richer model is evaluated through predictive checks and model comparison. Because the fitted models are generative for `TotalKg` conditional on age, sex, and bodyweight, they support prior and posterior predictive simulation rather than only fitted mean curves. The main result is that the nonlinear model finds a clear flattening of the bodyweight-strength relationship at heavier bodyweights and performs modestly better than the linear baseline, while the variational approximations perform substantially worse than HMC.

The final report is available at [reports/final/final_report.pdf](/Users/saumyajain/Desktop/Powerlifting_Bayesian_Inference/reports/final/final_report.pdf).

## Run

Before running the workflow, manually download the OpenPowerlifting bulk CSV and place it at `datasets/raw/openpowerlifting.csv`.

Then run the workflow in script order from the project root:

1. `Rscript scripts/01_prepare_openpowerlifting.R`
2. `Rscript scripts/02_eda.R`
3. `Rscript scripts/03_prior_likelihood_design.R`
4. `Rscript scripts/04_prepare_model_data.R`
5. `Rscript scripts/05_fit_baseline_model.R`
6. `Rscript scripts/06_check_baseline_model.R`
7. `Rscript scripts/07_fit_main_model.R`
8. `Rscript scripts/08_check_main_model.R`
9. `Rscript scripts/09_compare_models_and_summarize.R`

Large model artifacts are written to `artifacts/`, while lighter intermediate outputs and figures are written under `datasets/intermediate/` and `figures/`.
