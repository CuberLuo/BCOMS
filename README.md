# BCOMS

An Adaptive Bayesian Phase I/II Dose-Finding Design for Drug Combination Immunotherapy Trials Incorporating Survival Endpoints

## Overview

BCOMS is an adaptive Bayesian dose-finding design for combination immunotherapy trials in Phase I/II, incorporating survival endpoints. It is implemented in R and JAGS, supporting both simulation studies and practical trial design.

## Docker Deployment

### Prerequisites

- Install [Docker](https://docs.docker.com/get-docker/)

### Build Image

```bash
docker compose build
```

### Run Simulation

```bash
docker compose up
```

After completion, results will be saved in the `output/` directory.

## Local Run

### Requirements

- R >= 4.5.2
- JAGS (for Bayesian inference)

### Install Dependencies

```r
install.packages("renv")
renv::restore()
```

### Run Simulation

```bash
cd R
Rscript simulate.main.R
```
