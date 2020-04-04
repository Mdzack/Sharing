---
title: 'Chinook Salmon Abundance in North America'
author: "Maia Zack"
date: "01/11/2019"
html_document: default
---
##### Contact: mdzack@ucdavis.edu


#### <a href="#info"> 1. Current Information </a>

#### <a href="#hypothesese"> 2. Questions, Hypotheses, and Predictions </a>

#### <a href="#methods"> 3. Methods </a>

#### <a href="#results"> 4. Results </a>

#### <a href="#models"> 5. Extra Models - Methods </a>

#### <a href="#extras"> 6. Extra Modesl - Results </a>

#### <a href="#appendix"> 7. Appendix </a> 

![](Chinook_PIc.jpg)
_Image by Edmund Lee via Department of Commerce and Labor Bureau of Fisheries_

<a name="info"></a>

### Current Information
As seen in your fact sheet, Chinook salmon (or King Salmon) are important to local indigenous peoples, the salmon fishing industry, and both marine and freshwater ecosystems ( [WWF](https://www.worldwildlife.org/species/pacific-salmon)). Previous research has also shown a decrease in Chinook populations over recent years. They have been negatively affected by factors such as exploitation, climate change, and habitat degradation, loss, and impediments ( [NOAA](https://www.fisheries.noaa.gov/species/chinook-salmon-protected); [Ohlberger et al., 2018](https://onlinelibrary.wiley.com/doi/full/10.1111/faf.12272)).

<a name="hypothesese"></a>

### Research Question, Hypotheses, and Predictions
Question: How has Chinook salmon abundance changed over time in North America?

Hypotheses

- Hypothesis 1: Chinook salmon populations have declined over time in North America.

- Null hypothesis: Chinook salmon populations have not changed over time in North America.

- Hypothesis 2: Chinook salmon populations have increased over time in North America.

Predictions:

- I predict that chinook populations will continue to decline over the coming decades, if intervention is not taken.

(More information on planning of the report can be found in the [pre-registration](https://github.com/EdDataScienceEES/MaiaZack2019/blob/master/Challenges/Challenge3/Project_Planning/preregistration.md) and [workflow](https://github.com/EdDataScienceEES/MaiaZack2019/blob/master/Challenges/Challenge3/Project_Planning/Chinook_Workflow.md) documents.

<a name="methods"></a>

## Methods 
I evaluated chinook salmon population data from 33 studies from 1970 to 2012 in North America.  All data was taken from population data from the Living Planet Index ([The Living Planet Database](http://www.livingplanetindex.org/data_portal)).

I used the programming language R (version 3.6.1)  manipulate data and create graphs, tables, and models for this report. I also used the R package "MCMCglmm" to run the Bayesian model ([Hadfield, 2010](https://www.jstatsoft.org/article/view/v033i02/v33i02.pdf))The code I used can be found in Appendix 2.

``` r
chinook_mcmc <- MCMCglmm(scalepop ~ I(year-1970), random = ~year, data = chinook, nitt = 60000)
``` 
In the model, year acted as the predictor and explanatory variable, while scaled population was the dependent variable. I used a Gaussian distribution rather than a Poisson distribution because I used continuous, scaled population data, rather than the count data presented the population data (Figure 4, Appendix 1a).

I considered year a fixed and a random variable in the model.  Years are not truly independent because temporal autocorrelation exists.  For instance, if there is a disease that breaks out one year its affects on the population will reverberate into the next year.  I did not consider the location of population a random effect.  There did not appear to be any significant difference in sampling methods which might warrent its consideration as a random effect.  Also, the research question is to look at population abundance in North America, therefore, I did not want to ignore variance by location.

Finally, after I determined the model's structures I used trace plots of the bayesian model to the assess model's convergence before creating predictions (Figure 5, Appendix 1b).

<a name="results"></a>

## Results

From data taken from 33 locations between 1970 and 2012, I found that salmon populations are decreasing by a relative value of -0.0027 spawners per year (in scaled population data) in North America (slope = -0.0027, lower CI = -0.010, upper CI = 0.0040, Table 1, Figure 1). I determined hat the model did converge because the trace plot of the random and fixed effect, year, have a  "fuzzy caterpillar" appearance that is characteristic of a model that has converged (Figure 5, Appendix 1b). The confidence intervals from a Bayesian model did cover 0, indicating the results are not significant.  However, while the results are not siginificant and the null hypothesis cannot be rejected. It is worth noting that all models found a decrease in population.  Therefore, considering the threats to Chinook, I would recommend collecting more data and running further studies. This is the benefit of using a Bayesian model. It can be adjusted and re-run if more data is collected.

<img src="Chinook_Bayesian_Plot.png" width="700" />

__Figure 1__ The predictions of a Bayesian MCMCglmm model for Chinook salmon abundance plotted with [Living Planet Index](http://www.livingplanetindex.org/data_portal) data show a slight decrease over time. However, confidence intervals show this decrease is insignificant.

<br/><br/>

<table style="text-align:center"><caption><strong>Summary Results of a MCMCglmm Bayesian Model of Salmon Abundance between 1970 and 2012</strong></caption>
<tr><td colspan="8" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td>variable</td><td>post.mean</td><td>l.95..CI</td><td>u.95..CI</td><td>eff.samp</td><td>pMCMC</td><td>effect</td></tr>
<tr><td colspan="8" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">1</td><td>(Intercept)</td><td>-0.304</td><td>-0.457</td><td>-0.146</td><td>5,700</td><td>0.001</td><td>fixed</td></tr>
<tr><td style="text-align:left">2</td><td>I(year - 1970)</td><td>-0.003</td><td>-0.010</td><td>0.004</td><td>5,700</td><td>0.440</td><td>fixed</td></tr>
<tr><td style="text-align:left">3</td><td>year</td><td>0.050</td><td>0.022</td><td>0.082</td><td>5,942.019</td><td></td><td>random</td></tr>
<tr><td style="text-align:left">4</td><td>units</td><td>0.281</td><td>0.250</td><td>0.311</td><td>6,021.932</td><td></td><td>residual</td></tr>
<tr><td colspan="8" style="border-bottom: 1px solid black"></td></tr></table>
__Table 1__ MCMCglmm bayesian model outputs on the abundance of Chinook salmon.

<br/><br/>

<a name="models"></a>

## Extra Models - Methods 

I also created a frequentist general linear model. However, general linear models cannnot take into account for the autocorrelation that is present between years, so I did not use its resutls to evaluate the hypothesese.

``` r
chinook_lm <- lm(scalepop ~ I(year-1970), data = chinook)
```
<br/><br/>

A frequentist hierarchical linear model was also done.  Year was considered as both a random and a fixed effect, with the same reasoning as the Bayesian model.
``` r
chinook_lm4 <- lmer(scalepop ~ I(year-1970) + (1|year), data = chinook)
```
I decided to use the Bayesian model over the frequentist mixed model because the Bayesian model can be adjusted as the Living Planet Database is updated.  Also, the Bayesian model's confidence intervals are a clearer and more reliable indicator of significance over using p-values for the linear mixed model.

<br/><br/>

<a name="extras"></a>

## Extra Models - Results
#### General Linear Model
From data taken from 33 locations between 1970 and 2012, I found that salmon populations decreased relatively by -0.0027 spawners per year (in scaled population data) in North America (slope = -0.0022, standard error = +/- 0.0019, Table 2, Figure 2).  However, the magnitude of the decrease is quite low and the standard eror is quite high.

<br/><br/>


<img src="Chinook_Basic_Linear.png" width="700" />

__Figure 2__ The predictions of a general linear model for Chinook salmon abundance plotted with [Living Planet Index](http://www.livingplanetindex.org/data_portal) data show a slight decrease over time.  However, the magnitude of the decrease is quite low and the standard eror is quite high.

<br/><br/>

<table style="text-align:center"><caption><strong>Summary Results of General Linear Model of Salmon Abundance between 1970 and 2012</strong></caption>
<tr><td colspan="8" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Statistic</td><td>N</td><td>Mean</td><td>St. Dev.</td><td>Min</td><td>Pctl(25)</td><td>Pctl(75)</td><td>Max</td></tr>
<tr><td colspan="8" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">estimate</td><td>2</td><td>-0.157</td><td>0.220</td><td>-0.313</td><td>-0.235</td><td>-0.080</td><td>-0.002</td></tr>
<tr><td style="text-align:left">std.error</td><td>2</td><td>0.022</td><td>0.029</td><td>0.002</td><td>0.012</td><td>0.033</td><td>0.043</td></tr>
<tr><td style="text-align:left">statistic</td><td>2</td><td>-4.224</td><td>4.370</td><td>-7.314</td><td>-5.769</td><td>-2.679</td><td>-1.134</td></tr>
<tr><td style="text-align:left">p.value</td><td>2</td><td>0.129</td><td>0.182</td><td>0</td><td>0.1</td><td>0.2</td><td>0</td></tr>
<tr><td colspan="8" style="border-bottom: 1px solid black"></td></tr></table>

__Table 2__ Basic linear model outputs on the abundance of Chinook salmon.

<br/><br/>
<br/><br/>

#### Mixed Linear Model

From data taken from 33 locations between 1970 and 2012, we found that salmon populations are decreasing by -0.0027 spawners per year (in scaled population data) in North America (slope = -0.0027, standard error = +/- 0.0034, Table 3).  The error of the model is quite large, so I would consider the results insignificant.

<img src="Chinook_Mixed_Plot.png" width="700" />

__Figure 3__ The predictions of a linear mixed model for Chinook salmon abundance plotted with [Living Planet Index](http://www.livingplanetindex.org/data_portal) data show a slight decrease over time.  However, the magnitude of the decrease is quite low and the standard eror is quite high.

<br/><br/>

<table style="text-align:center"><caption><strong>Summary Results of Linear Mixed Model of Salmon Abundance between 1970 and 2012</strong></caption>
<tr><td colspan="8" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Statistic</td><td>N</td><td>Mean</td><td>St. Dev.</td><td>Min</td><td>Pctl(25)</td><td>Pctl(75)</td><td>Max</td></tr>
<tr><td colspan="8" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">estimate</td><td>4</td><td>0.111</td><td>0.351</td><td>-0.302</td><td>-0.078</td><td>0.297</td><td>0.529</td></tr>
<tr><td style="text-align:left">std.error</td><td>2</td><td>0.041</td><td>0.053</td><td>0.003</td><td>0.022</td><td>0.059</td><td>0.078</td></tr>
<tr><td style="text-align:left">statistic</td><td>2</td><td>-2.335</td><td>2.186</td><td>-3.881</td><td>-3.108</td><td>-1.562</td><td>-0.789</td></tr>
<tr><td colspan="8" style="border-bottom: 1px solid black"></td></tr></table>

__Table 3__ Mixed linear model outputs on the abundance of Chinook salmon.

<br/><br/>
<br/><br/>

<a name="appendix"></a> 

### Appendix

#### Appendix 1 - Methods

##### 1a. Scaled Population Histogram

<img src="Chinook_hist.png" width="700" />

__Figure 4__ Histogram plots of scaled popultation counts for each population location from 1970 to 2012.
<br/><br/>

##### 1b. Trace and Density Plots of MCMCglmm model of Chinook Salmon Abundance

<img src="Chinook_Trace_Plots.png" width="700" />

__Figure 5__ Trace and density plots of a MCMCglmm model of chinook salmon abundance used to assess model convergence.

<br/><br/>
<br/><br/>

#### Appendix 2 - Code

``` r
##%######################################################%##
#                                                          #
####          Challenge 3- Statistical Modeling         ####
####    Maia Zack 1/11/2019, University of Edinburgh    ####
####           Contact: mdzack@ucdavis.edu              ####
#                                                          #
##%######################################################%##

### Pre-registration: 
### Final Report: https://github.com/EdDataScienceEES/challenge3-statistical-modelling-2019-Mdzack/blob/master/Report/Chinook_Report.md

## Load Libraries ----
library(tidyverse)
library(lme4)
library(MCMCglmm)
library(ggplot2)
library(ggeffects)
library(broom)
library(stargazer)

## Load Living Planet Data ----
load("Data/LPI_species.Rdata")

## Data Manipulation and Exploration----

# Explore data
head(LPI_species)
str(LPI_species)

# Select for necessary data
chinook <- LPI_species %>% 
  # Filter for "Chinook Salmon"
  filter(Common.Name == "Chinook salmon") %>% 
  
  # Get rid of unnecessary columns
  select(id, Location.of.population, year, scalepop)

# Change year to a character
chinook$year <- as.character(chinook$year)

# Get rid of "x"'s before year #'s
chinook$year <- parse_number(chinook$year)

# Turn year back to a numeric 
chinook$year <- as.numeric(chinook$year)

# Make a new column for year for linear mixed model
chinook <- mutate(chinook, year2 = I(year-1970))
## I don't know why, but if I create a duplicated year column and use that in my linear model my mixed linear model predictions includes standard error and confidence levels, however, when I use only year ggpredict doesn't include these values


## Frequentist Linear model ----

# Explore distribution of data
(chinook_hist <- ggplot(chinook, aes(x = scalepop)) +
   facet_wrap(~Location.of.population, nrow = 10) +
   geom_histogram() + 
   theme_bw() + 
   theme(panel.grid = element_blank()) +
   labs(title = "Scaled Population Distributions",
        y = "\nCount", x = "Chinook Population (salmon count, scaled)\n") +
   ggsave(filename = "Report/Chinook_hist.png", device = "png", height = 30, width = 20))
# Data is normally dist within locations of pop, so use Gaussian dist

# Creating general linear model with no random effects 
chinook_lm <- lm(scalepop ~ I(year-1970), data = chinook)

# Summary of simple linear model
summary(chinook_lm)

# Create a table of the results from general linear model
lm.table <- as.data.frame(tidy(chinook_lm))

# Create report ready table
stargazer(lm.table, type = "html", out = "Report/lm_table.html",
          title = "Summary Results of General Linear Model of Salmon Abundance between 1970 and 2012")

# Checking the assumptions through plotting 
plot(chinook_lm) # There are some severe deviations in the QQplot about halfway through

## Visualize general linear model ----
(ggplot(chinook, aes(x = year, y = scalepop)) +
   geom_point(colour = "#3182bd", alpha = 0.5) +
   geom_smooth(method = lm, colour = "#3182bd", fill = "#9ecae1", alpha = 0.6) +
   theme_bw() + 
   theme(panel.grid = element_blank()) +
   annotate("text", x = 2000, y = 0.90, label = "Slope = -0.0022, Std. Error = ±0.0019") + 
   scale_x_continuous(breaks = c(1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010)) +
   labs(x = "\nYear", y = "Chinook Population (salmon count, scaled)\n",
        title = "General Linear Model Predictions of Salmon Abundance") +
  ggsave(filename = "Report/Chinook_Basic_Linear.png", device = "png", width = 8, height = 6))


## Frequentist Mixed Effects Linear Model ----
# use lmer because I have continuous data, rather than glmer is for other response types to create linear mixed model 
chinook_lm4 <- lmer(scalepop ~ year2 + (1|year), data = chinook)

# Review results of model
summary(chinook_lm4)

# Look at the results of model
plot(chinook_lm4)

# Tidy results in a table
chinook.lm4.table <- as.data.frame(tidy(chinook_lm4))

# Create report ready table
stargazer(chinook.lm4.table, type = "html", out = "Report/lm4_table.html",
          title = "Summary Results of Linear Mixed Model of Salmon Abundance between 1970 and 2012")

## Visualize mixed linear model ----
# Create model predictions 
lm4_predictions <- ggpredict(chinook_lm4, terms = c("year2"))

# Create plot 
(ggplot() +
    geom_line(data = lm4_predictions, aes(x = x, y = predicted), alpha = 1,
              size = 2, colour = "#D92323")+
    geom_ribbon(data = lm4_predictions, aes(ymin = conf.low, ymax = conf.high, x = x),
                alpha = 0.3, fill = "#D92323") +
    geom_point(data = chinook, aes(x = year2, y = scalepop), alpha = 0.3, 
               size = 2, colour = "#F7A6A6") +
    theme_bw() + 
    theme(panel.grid = element_blank()) +
    annotate("text", x = 35, y = 0.90, label = "Slope = -0.0027, Std. Error = ±0.0034") +  
    scale_y_continuous(limits = c (-1, 1)) +
    scale_x_continuous(limits = c (0, 45)) +
    #breaks = c(1970, 1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015)) +
    labs(x = "\nYear (-1970)", y = " Chinook Population (salmon count, scaled)\n ",
         title = "Linear Mixed Model Predictions of Salmon Abundance") +
    ggsave(filename = "Report/Chinook_Mixed_Plot.png", device = "png", width = 8, height = 6))


## Bayesian Linear model ----

# Create Bayesian model for abundance
chinook_mcmc <- MCMCglmm(scalepop ~ I(year-1970), random = ~year, data = chinook, nitt = 60000)

# Review results of the model
summary(chinook_mcmc)

# Look at the trace plots to assess convergence 
png("Report/Chinook_Trace_Plots.png")
plot(chinook_mcmc$VCV)
dev.off()

# Tidy summary data for table
bayesian.table <- as.data.frame(tidyMCMC(chinook_mcmc))

## Visualize Bayesian model ----
# Create predcitions for Bayesian model
mcmc_pred <- ggpredict(chinook_mcmc, terms = c("year"))

# Create predictions plot
(ggplot() +
    geom_line(data = mcmc_pred, aes(x = x, y = predicted), alpha = 1,
              size = 2, colour = "#07571B") +
    geom_ribbon(data = mcmc_pred, aes(ymin = conf.low, ymax = conf.high, x = x),
                alpha = 0.3, fill = "#07571B") +
    geom_point(data = chinook, aes(x = year, y = scalepop), alpha = 0.3, 
               size = 2, colour = "#3B7D46") +
    theme_bw() + 
    theme(panel.grid = element_blank()) +
    annotate("text", x = 2000, y = 0.90, label = "Slope = -0.0027, lower CI = -0.010, upper CI = 0.0040") +  
    scale_y_continuous(limits = c (-1, 1)) +
    scale_x_continuous(limits = c (1970, 2015),
                       breaks = c(1970, 1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015)) +
    labs(x = "\nYear", y = "Chinook Population (salmon count, scaled)\n",
         title = "Bayesian Model Predictions of Salmon Abundance") +
  ggsave(filename = "Report/Chinook_Bayesian_Plot.png", device = "png", width = 8, height = 6))

# Create function to clean model
clean.MCMC <- function(x) {
  sols <- summary(x)$solutions  ## pull out relevant info from model summary
  Gcovs <- summary(x)$Gcovariances
  Rcovs <- summary(x)$Rcovariances
  fixed <- data.frame(row.names(sols), sols, row.names = NULL)  ## convert to dataframes with the row.names as the first col
  random <- data.frame(row.names(Gcovs), Gcovs, row.names = NULL)
  residual <- data.frame(row.names(Rcovs), Rcovs, row.names = NULL)
  names(fixed)[names(fixed) == "row.names.sols."] <- "variable"  ## change the columns names to variable, so they all match
  names(random)[names(random) == "row.names.Gcovs."] <- "variable"
  names(residual)[names(residual) == "row.names.Rcovs."] <- "variable"
  fixed$effect <- "fixed"  ## add ID column for type of effect (fixed, random, residual)
  random$effect <- "random"
  residual$effect <- "residual"
  modelTerms <- as.data.frame(bind_rows(fixed, random, residual))  # merge it all together
}

# Clean MCMCglmm model
mcmc_tidied <- clean.MCMC(chinook_mcmc)

# Create report ready table
stargazer(mcmc_tidied, type = "html", summary = FALSE, out = "Report/MCMC_Table.html",
          title = "Summary Results of a MCMCglmm Bayesian Model of Salmon Abundance between 1970 and 2012")

```






