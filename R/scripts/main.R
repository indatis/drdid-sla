getwd()
setwd("~/PHD UNICAS/Project - DRDID extension/drdid-sla-main")
list.files()
rm(list=ls())

library(DRDID)
library(dplyr)
library(tidyr)
library(ggplot2)

source("R/functions/drdid_att_panel.R")
source("R/functions/sla_runner.R")
source("R/functions/sla_diagnostics.R")
source("R/functions/sla_plots.R")
source("R/functions/simulate_panel.R")
source("R/functions/sla_outputs.R")

source("R/scripts/01_lalonde_example.R")
source("R/scripts/02_simulation_strong.R")
source("R/scripts/03_simulation_null.R")
source("R/scripts/04_krueger_example.R")
source("R/scripts/05_adh_example.R")
source("R/scripts/06_callaway_example.R")

source("R/functions/make_sla_gifs.R")


