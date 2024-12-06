# Tomographic reconstruction of a disease transmission landscape via GPS recorded random paths

Identifying areas in a landscape where individuals have a higher probability of becoming infected with a pathogen is a crucial step towards disease management. We perform a novel epidemiological tomography for the estimation of landscape propensity to disease infection, using GPS animal tracks in a manner analogous to  tomographic techniques in Positron Emission Tomography (PET). Our study data consists of individual tracks of white-tailed deer (Odocoileus virginianus) and three exotic Cervid species moving freely in a 172-ha high-fenced game preserve over given time periods. A serological test was performed on each individual to measure the antibody concentration of epizootic hemorrhagic disease viruses (EHDV) at the beginning and at the end of each tracking period. EHDV is a vector-borne viral disease indirectly transmitted between ruminant hosts by biting midges (Culicoides). We model the data as a binomial linear inverse problem, where spatial coherence is enforced  with a total variation regularization. The smoothness of the reconstructed propensity map is selected by the quantile universal threshold, which can also  test the null hypothesis that the propensity map is spatially constant. We apply our method to simulated and real data, showing good statistical properties during simulations and consistent results and interpretations compared to intensive field estimations.

## Usage instructions

1. Clone repo.
2. Install cvx https://cvxr.com/cvx/doc/install.html
3. Add repo to MATLAB path.
4. Main function is `tvQUT` which estimates propensity based on design matrix $X$, response $Y$ and total variation matrix $D$. It internally calculates regularization parameter using function `lambdaQUT`.
5. All desing matrices $X$, responses $Y$ and matrices $D$ used to estimate EHDV propensities are available in folder `data`.
6. To replicate paper results on deer data, run `deer/deermoves_bootstrap_plus.m`.
7. To replicate all single simulations run `simulations/simulate_main.m`.
8. To replicate simulations with uncertainty quantification run `simulations/simulate_main_bootstrap_parallel.m`.
9. To replicate simulation tables and plots use `simulations/results_single.m` `simulations/results_coverage.m` after step 8.
10. To replicate splines results, run R files in folder `splines`


