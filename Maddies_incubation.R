rm(list = ls())
pacman::p_load(tidyverse, readxl, patchwork)
theme = theme(axis.ticks = element_line(color = "black"),
              axis.title = element_text(size = 12), 
              axis.text = element_text(color = "black", size = 10),
              text = element_text(color = "black", size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              panel.grid = element_blank())

# LOAD DATA ----
incubation = read_csv("data/Incubation 2 Cleaned Results.csv")[c(14,18,22), 2:13] |>
  mutate(day = c(3, 7, 21))
names(incubation)[1:12] = c("Na", "Mg", "K", "Ca", "Li", "B", "Si", "Mn", "Zn", "Sr", "Ba", "pH")
# convert from mg/L to mol/L
incubation = incubation |>
  mutate(Na = Na / (1e3*23),
         Mg = Mg / (1e3*24),
         K = K / (1e3*39),
         Ca = Ca / (1e3*40),
         Li = Li / (1e3/6.94),
         B = B / (1e3*10.81),
         Si = Si / (1e3*28),
         Mn = Mn / (1e3*55),
         Zn = Zn / (1e3*65),
         Sr = Sr / (1e3*87.62),
         Ba = Ba / (1e3*137.327))
element = incubation |>
  pivot_longer(cols = 1:12,
               names_to = "name",
               values_to = "value") |>
  filter(name != "pH") |>
  group_by(name) |>
  summarise(mean = mean(value))
ggplot(element, aes(x = name, y = mean)) +
  geom_point()
element

# Model setup ----
source('model/box_model_ode_major_elements.R')
time = seq(0, 1e4, 1)
parms["F_in"] = 0 # 3e-5
parms["F_evap"] = 0
parms["k_degas"] = 0
parms["k_diss"] = 0
parms["Mn_rd"] = 5e-3
results = ode(y = state,
              times = time,
              parms = parms,
              func = soil_water_chemistry)
results = as.data.frame(results)
calcite = results |> filter(Jp > 0)
par(mfrow = c(3,2), mar = margin(3, 3, 1, 1))
plot(results$time, results$pH, type = "l", xlab = "time (s)", ylab = "pH", mgp = c(2, .8, 0))
plot(results$time, results$Alk, type = "l", xlab = "time (s)", ylab = "Alk (mol/L)", mgp = c(2, .8, 0))
plot(results$time, results$DIC_sw, type = "l", xlab = "time (s)", ylab = "DIC (mol/L)", mgp = c(2, .8, 0))
plot(calcite$time, 1e3*calcite$Mn_sw, type = "l", xlab = "time (s)", ylab = "Mn (mmol/L)", mgp = c(2, .8, 0))
plot(calcite$time, calcite$Jp, type = "l", xlab = "time (s)", ylab = "calcite (mol/s)", mgp = c(2, .8, 0))
plot(results$time, results$HCO3_rd, type = "l", xlab = "time (s)", ylab = "Mn reduction rate (mol/s)", mgp = c(2, .8, 0))

ggplot(results) +
  geom_point(aes(x = time, y = Mn_sw*1e3), size = 1) +
  theme_bw() + theme +
  theme(legend.title = element_text(size = 12, margin = margin(b = 15))) +
  labs(y = "Mn (mmol/L)",
       x = "time (h)") +
  scale_y_continuous(limits = c(4e-3, 6e-3))
  
