# Tutorials

LASER provides a library of reusable epidemiological components as well as canonical transmission models: SI, SIS, SIR, SIRS, SEIR, SEIRS, spatial networks, and more. This flexible framework is customizable, able to be applied to a wide variety of epidemiological modeling questions. From simple transmission modes to complex spatial dynamics, there are numerous ways to combine and apply LASER modules.

The following code examples and interactive Jupyter notebooks are designed to provide implementation suggestions and to help familiarize users with core LASER functionality. We encourage new users to explore LASER through the notebooks, and challenge them to create new methods for answering similar questions.

The tutorials begin with a simple SIR model, with sample code demonstrating how to configure LASER to implement the `SIRModel` class. As the tutorial progresses, complexity is added to the model including spatial components and human migration.

The notebook section contains a suite of interactive notebooks that progress through transmission modes, starting from SI models with no births and sequentially adding additional dynamics to create SIR models with vital dynamics and migration. In addition to the series of notebooks on transmission dynamics, there are notebooks which explore specific LASER features such as births, mortality, and seasonality. We recommend working through the transmission component notebooks prior to the feature notebooks.

Suggested order:

1. [SI model with no demographics](notebooks/01_SI_nobirths_logistic_growth.ipynb)
1. [SI model with constant population demographics](notebooks/02_SI_wbirths_logistic_growth.ipynb)
1. [SIS model with no demographics](notebooks/03_SIS_nobirths_logistic_growth.ipynb)
1. [Outbreak size in the SIR model](notebooks/04_SIR_nobirths_outbreak_size.ipynb)
1. [Average age at infection in the SIR model](notebooks/05_SIR_wbirths_age_distribution.ipynb)
1. [Intrinsic periodicity of the SIR model](notebooks/06_SIR_wbirths_natural_periodicity.ipynb)
1. [Exploring the critical community size of an SIR model](notebooks/07_SIR_CCS.ipynb)
1. [The relationship between coupling and incidence correlation in a 2-patch model](notebooks/08_2patch_SIR_wbirths_correlation.ipynb)
1. [Modeling the spread of rabies in one dimension](notebooks/09_rabies_diffusion_1D.ipynb)
1. [Periodicity of measles in England and Wales](notebooks/10_EW_periodicity.ipynb)

Once you have completed the above notebooks, please explore the following notebooks in any desired order to learn how to implement and integrate custom components:

- [England and Wales measles analysis](notebooks/EW_analysis.ipynb)
- [SEI and SEIS model implementations](notebooks/SEI_and_SEIS_implementation.ipynb)
- [Explore the vital dynamics births components](notebooks/births.ipynb)
- [Test constant population components](notebooks/constant_pop.ipynb)
- [Numba compatible distributions](notebooks/distributions.ipynb)
- [Grid function examples](notebooks/grid_examples.ipynb)
- [Explore the vital dynamics mortality components](notebooks/mortality.ipynb)
- [Explore routine immunization](notebooks/ri_exploration.ipynb)
- [Seasonality in transmission](notebooks/seasonality.ipynb)
