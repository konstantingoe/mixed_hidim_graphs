
# Packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  stargazer,
  huge,
  future.apply,
  parallel,
  dplyr,
  stats,
  MASS,
  gdata,
  Matrix,
  polycor,
  corpcor,
  plyr)
