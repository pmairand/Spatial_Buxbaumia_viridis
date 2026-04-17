# Role of dispersal and microhabitat on the spatial distribution of a coarse wood-associated moss at the local scale
 
**Master Thesis - Paul Mairand**  
Office National des Forêts (ONF) - 2025

 
## Study species
 
*Buxbaumia viridis* (green shield-moss) is one of the few bryophytes protected in France and Europe (Habitats Directive, Annex II 92/43/EEC ; Bern Convention, Annex I ; ECCB 1995 - VU). It is a saprolignicolous moss that colonizes coarse woody debris (CWD) in coniferous montane forests, requiring advanced stages of wood decomposition and dense canopy cover.

 
## Objectives
 
Most existing studies on *B. viridis* have focused on chorological distribution or habitat niche modelling at broad spatial scales.This study addresses a gap in the literature by investigating the **fine-scale spatial structure** of the species within individual forest plots in the Massif Central (France).
 
The core question is : **does *B. viridis* simply follow the deadwood resource, or does it exhibit independent spatial processes** (dispersal limitation, conspecific facilitation) that shape its colonisation pattern beyond what microhabitat quality alone can predict?

 
## Study area
 
Four forest plots in the **Massif Central, France** (sites M2, M3, M4, M5), covering montane coniferous forests with varying colonisation rates (6-22%).

 
## Methods
 
### Exploratory analysis - Second-order statistics
 
- **Pair Correlation Function (PCF)** : quantifies spatial aggregation of deadwood substrates relative to complete spatial randomness (CSR). Values g(r) > 1 indicate aggregation at scale r.
- **Mark Connection Function (MCF)** : estimates the probability that two substrates at distance r are both colonised by *B. viridis*, under a random labelling null hypothesis. Departures above the simulation envelope indicate clustering of colonised substrates independent of the underlying CWD distribution.
### Point Process Models - PPM
 
Inhomogeneous Poisson Point Process Models (PPM) fitted by maximum likelihood with three spatially continuous covariates:
 
$$\log \lambda(u) = \beta_0 + \beta_1 \cdot \text{dens10}(u) + \beta_2 \cdot \text{alt}(u) + \beta_3 \cdot \text{sap}(u)$$
 
- `dens10` : local deadwood density weighted by diameter (10 m bandwidth)
- `alt` : smoothed local altitude
- `sap` : substrate saproxylation stage
### Joint INLA model - Bayesian spatial modelling
 
A joint latent Gaussian model fitted via **INLA** combined with the **SPDE approach**, coupling two processes through a shared spatial random field:
 
**Process 1 - Log-Gaussian Cox Process (LGCP)** for deadwood distribution:
 
$$\log \lambda_1(s) = \mu_1 + w_1(s)$$
 
**Process 2 - Binomial model** for *B. viridis* colonisation probability:
 
$$\text{logit}(p_i) = \alpha + \beta_{\text{diam}} \cdot \text{diam}_i + \beta_{\text{sapmax}} \cdot SapMAX_i + \beta_{\text{ess}} \cdot [\text{conifer}] + \beta_{\text{sup}} \cdot [\text{stump}] + \beta_s \cdot w_1(s_i) + w_2(s_i)$$
 
The coupling term $\beta_s \cdot w_1(s_i)$ tests whether areas of high deadwood density independently favour colonisation once individual substrate covariates are accounted for.

 
## Main results
 
- **Deadwood** is strongly aggregated at the metric scale (r < 2-3 m) at all four sites, driven by localised tree fall events.
- **MCF** reveals independent spatial clustering of colonised substrates only at the most heavily colonised site (M3, 22%), consistent with short-distance spore dispersal or local source effects.
- **Saproxylation stage** is the dominant driver of colonisation intensity across all sites (β = 1.54-2.44, p < 0.001).
- **Joint INLA** confirms that local deadwood density has no independent effect on colonisation probability (β_s ≈ 0) once substrate-level covariates are controlled.

 
## Repository structure
 
```
01_Statistics/
├── 01_R_Scripts/
│   ├── 00_Bux_Articles_plot.R       # Bibliographic data visualisation
│   ├── 01_Data_management.R         # Data cleaning and preparation
│   ├── 02_PCF_MCF_PPM.R             # Exploratory point pattern analysis
│   ├── 03_INLA.R                    # Simple INLA model
│   ├── 04_joint_INLA.R              # Joint INLA model (LGCP + Binomial)
│   ├── 05_joint_INLA_validation.R   # Model validation
│   └── Functions/                   # Custom R functions
│       ├── 01_PCF_fun.R
│       ├── 02_MCF_fun.R
│       ├── 03_INLA_fun.R
│       ├── 04_joint_INLA_fun.R
│       └── 05_joint_INLA_valid_fun.R
└── 04_Quarto/
    └── Buxbaumia_viridis-Stage_P_Mairand.qmd   # Full analysis report
```
 
>  Raw data are not included in this repository.

 
## Main R packages
 
- [`spatstat`](https://spatstat.org/) - point pattern analysis (PCF, MCF, PPM)
- [`R-INLA`](https://www.r-inla.org/) - Bayesian spatial modelling

 
## Key references
 
- Rue, Martino & Chopin (2009) - INLA
- Lindgren, Rue & Lindström (2011) - SPDE approach
- Krainski et al. (2018) - *Advanced Spatial Modeling with SPDE*
- Simpson et al. (2016) - PC priors
- Kropik et al. (2020) - *B. viridis* ecology and dispersal
- Offerhaus et al. (2019) - *B. viridis* distribution in Europe

## Data availability

The data will be made available on Zenodo as soon as this study is published.

