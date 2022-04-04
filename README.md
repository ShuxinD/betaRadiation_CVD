# Read Me

This repository is the data analysis part in the paper of "Long-term Exposure to Particle Radioactivity (Gross Beta Activity) and PM2.5 and Cardiovascular Diseases Mortality" (may update name later)

## Data

### Health outcome of interest

The annual death rate per ZCTA from 2001/01 to 2015/12 in MA (from MA the Department of Public Health) \<-- death count per ZCTA/population count per ZCTA

-   [01convertMAdeath_ZCTA.R](https://github.com/ShuxinD/betaRadiation_CVD/blob/main/codes/01convertMAdeath_ZCTA.R): convert MA death record from one-row-per-person to one-row-per-zcta
-   population count per ZCTA is from Joel's server from census data

### PR and PM<sub>2.5</sub>

Annual average of PR and PM<sub>2.5</sub>: aggregate from daily 1km<sup>2</sup> grid to annual ZCTA.

-   Annual averages of PR per ZCTA are from Longxiang
-   Annual averages of PM<sub>2.5</sub> per ZIP code are from Joel's server

### Meterological data: temperature

Daily temperature from Daymet

-   [02seasonaltempDaymet_ZCTA.R](https://github.com/ShuxinD/betaRadiation_CVD/blob/main/codes/02seasonaltempDaymet_ZCTA.R): convert daily temperature to seasonal average per ZCTA for each year

### Other covariates

Socioeconomics status and general life habits covariates are from an existing dataset:

-   Census: "Median household income", "Median value of house", "Percentage of Hispanic", "Percentage of black", "Pencentage below poverty line", "Percentage without high school diploma", "Population density"
-   the Behavioral Risk Factor Surveillance System (BRFSS): "Mean BMI", "Smoking rate"

**We use [03mergeDatasetsALL.R](https://github.com/ShuxinD/betaRadiation_CVD/blob/main/codes/03mergeDatasetsALL.R) to construct the final dataset for analysis.**

## Statistics

### Descriptive

[04descriptiveStats.R](https://github.com/ShuxinD/betaRadiation_CVD/blob/main/codes/04descriptiveStats.R)

-   generate the correlation table among each variable
-   generate the distribution (mean, SD, quartiles) for all covariates

### Models

[05modelling.R](https://github.com/ShuxinD/betaRadiation_CVD/blob/main/codes/05modelling.R)

#### DID

Differences-in-difference:

``` {.r}
mod <- gam(deathcount ~ Beta + pm25 + summer_tmean + winter_tmean + as.factor(year) + as.factor(ZCTA5CE10) + offset(log(pcount)),
               data = dt,
               family = quasipoisson(link = "log"))
tb <- summary(mod)$p.table
```

#### GLMM

Generalized linear mix-effect model:

``` {.r}
mod <- gamm(deathcount ~ Beta + pm25 + summer_tmean + winter_tmean + medhouseholdincome + medianhousevalue + hispanic + pct_blk + poverty + education + popdensity + mean_bmi + smoke_rate + as.factor(year) + offset(log(pcount)),
              data = dt,
              random = list(ZCTA5CE10=~1),
              family = quasipoisson(link = "log"))
```

**The results for each pollutant were plotted via [06plotting.R](https://github.com/ShuxinD/betaRadiation_CVD/blob/main/codes/06plotting.R)**

### Model: Interaction between PR and PM

``` {.r}
 mod <- gam(deathcount ~ pm25 + Beta + I(pm25*Beta) + 
               summer_tmean + winter_tmean + as.factor(year) + as.factor(ZCTA5CE10) + offset(log(pcount)),
             data = dt,
             family = quasipoisson(link = "log"))

mod <- gamm(deathcount ~pm25 + Beta + I(pm25*Beta) + 
                summer_tmean + winter_tmean + medhouseholdincome + medianhousevalue + hispanic + pct_blk + poverty + education + popdensity + mean_bmi + smoke_rate + as.factor(year) + offset(log(pcount)),
              data = dt,
              random = list(ZCTA5CE10=~1),
              family = quasipoisson(link = "log"))
```
