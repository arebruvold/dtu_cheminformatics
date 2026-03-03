# Usage

# setup
library(tidyverse)

# list of compounds

cmps <- c(
  "perfluorooctanoic acid",
  "PFOA perfluorooctanoic acid",
  "PFOA cpd",
  "pentadecafluorooctanoic acid"
) %>%
  as_tibble() %>%
  rename(compound = 1)

cmps_info <- cmps %>%
  mutate(pcprops = pcprop_from_name_vec(compound)) %>%
  unnest(pcprops)
