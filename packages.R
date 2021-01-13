
# Packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  stargazer,
  huge,
  future.apply,
  dplyr,
  stats,
  MASS,
  Matrix,
  polycor,
  plyr,
  # those beneath always load last
  tidyverse)
