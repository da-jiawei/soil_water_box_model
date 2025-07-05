rm(list = ls())
pacman::p_load(tidyverse, readxl, patchwork)

# load data 
dat = read_xlsx("data/DanceBayou_all_carb_data.xlsx") |>
  select(ID, site, depth, d13C, d18O, `Sr, L1`, `Mg, L1`, `Mn, L1`)
names(dat) = c("ID", "site", "depth", "d13", "d18", "SrCa_mass", "MgCa_mass", "MnCa_mass")
dat = dat %>%
  filter(site == "1") %>%
  mutate(MgCa = 1e3 * MgCa_mass * 40.078 / 24.305,
         SrCa = 1e3 * SrCa_mass * 40.078 / 87.62,
         MnCa = 1e3 * MnCa_mass * 40.078 / 54.94,
         depth = as.numeric(depth))

ggplot(dat, aes(x = d13, y = d18, color = ID)) +
  geom_point(shape = 21) +
  scale_y_reverse()
