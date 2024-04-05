# Spatial estimation of virus infection propensity in hosts determined from GPS-based space-time locations

Identifying areas in a landscape where individuals have a higher probability of becoming infected with a pathogen is a crucial step towards disease management. Our study data consists of GPS-based tracks of individual white-tailed deer (\textit{Odocoileus virginianus}) and three exotic Cervid species moving freely in a 172-ha high-fenced game preserve over given time periods. A serological test was performed on each individual to measure the antibody concentration of epizootic hemorrhagic disease virus (EHDV) for each of three serotypes (EHDV-1, -2, and -6) at the beginning and at the end of each tracking period. EHDV is a vector-borne viral disease indirectly transmitted between ruminant hosts by biting midges (\textit{Culicoides} spp.). The purpose of this study is to estimate the spatial distribution of infection propensity by performing an epidemiological tomography of a region using tracers. We model the data as a binomial linear inverse problem, where spatial coherence is enforced  with a total variation regularization. The smoothness of the reconstructed propensity map is selected by the quantile universal threshold, which can also  test the null hypothesis that the propensity map is spatially constant. We apply our method to simulated and real data, showing good statistical properties during simulations and consistent results and interpretations compared to intensive field estimations.

## Usage instructions

1. Clone repo.
2. Install cvx https://cvxr.com/cvx/doc/install.html
3. Add repo to MATLAB path.
4. Main function is `tvQUT` which estimates propensity based on design matrix $X$, response $Y$ and total variation matrix $D$. It internally calculates regularization parameter using function `lambdaQUT`.
5. All desing matrices $X$, responses $Y$ and matrices $D$ used to estimate EHDV propensities are available in folder `data`.
6. To replicate paper results on deer data, run `deer/deermoves.m`.
7. To replicate all simulations run `simulations/simulate_main.m` and `simulations/simulate_main_text.m`.
8. To replicate simulation plots use `simulations/plot_effect_n.m`, `simulations/plot_effect_nx.m` and `simulations/plot_effect_tn.m` after step 7.


