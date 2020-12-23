
# Packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  pcalg,
  mvtnorm,
  weights,
  ks,
  stargazer,
  e1071,
  DescTools,
  mipfp,
  foreign,
  glasso,
  corpcor,
  huge,
  glmnet,
  randomForest,
  igraph,
  mice,
  future.apply,
  VIM,
  dplyr,
  naniar,
  Hmisc,
  MASS,
  data.table,
  polycor,
  plyr,
  rlang,
  Rgraphviz,
  # those beneath always load last
  ggplot2,
  tidyverse)
