rm(list = ls())
pacman::p_load(tidyverse, readxl, ggpubr, scales, patchwork)
theme = theme(axis.ticks = element_line(color = "black"),
              axis.title = element_text(size = 12), 
              axis.text = element_text(color = "black", size = 10),
              text = element_text(color = "black", size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              panel.grid = element_blank())

# Load data and model ---- 
DB = read_xlsx("data/DanceBayou_all_carb_data.xlsx")
DB1 = data.frame(DB$ID, DB$site, DB$depth, DB$d13C, DB$d18O, DB$`Sr, L1`, DB$`Mg, L1`, DB$`Mn, L1`)
names(DB1) = c("ID", "site", "depth", "d13", "d18", "SrCa_mass", "MgCa_mass", "MnCa_mass")
DB1 = DB1 %>%
  filter(site == "1") %>%
  mutate(MgCa = 1e3 * MgCa_mass * 40.078 / 24.305,
         SrCa = 1e3 * SrCa_mass * 40.078 / 87.62,
         MnCa = 1e3 * MnCa_mass * 40.078 / 54.94)
DB1$depth = as.numeric(DB1$depth)
cat("\014")

# Simple version ----
# no Fin, no CO2 degassing nor dissolution, constant Mn reduction rate
source('model/box_model_ode.R')
time = seq(0, 1.9e4, 10)
parms["F_in"] = 0 
# parms["F_evap"] = 0
parms["k_degas"] = 0
parms["k_diss"] = 0
state["DIC_sw"] = 2.5e-3
state["Mn_sw"] = 0
parms["Mn_rd"] = 5e-8
# parms["lamda"] = 1
results = ode(y = state,
              times = time,
              parms = parms,
              func = soil_water_chemistry)
results = as.data.frame(results) |> filter(time > 100)
calcite = results |> filter(Jp > 0)

par(mfrow = c(3,2), mar = margin(3, 3, 1, 1))
plot(results$time, results$pH, type = "l", xlab = "time (s)", ylab = "pH", mgp = c(2, .8, 0))
plot(results$time, results$Alk * 1e3, type = "l", xlab = "time (s)", ylab = "Alk (mmol/L)", mgp = c(2, .8, 0))
plot(results$time, results$DIC_sw * 1e3, type = "l", xlab = "time (s)", ylab = "DIC (mmol/L)", mgp = c(2, .8, 0))
plot(calcite$time, calcite$MnCa_c * 1e3, type = "l", xlab = "time (s)", ylab = expression("Mn/Ca"[c]*" (mmol/mol)"), mgp = c(2, .8, 0))
plot(calcite$time, calcite$Jp, type = "l", xlab = "time (s)", ylab = "calcite (mol/s)", mgp = c(2, .8, 0))
plot(results$time, results$Mn_rd, type = "l", xlab = "time (s)", ylab = "Mn reduction rate (mol/s)", mgp = c(2, .8, 0))
# plot(calcite$time, calcite$SrCa_c * 1e3, type = "l", xlab = "time (s)", ylab = expression("Sr/Ca"[c]*" (mmol/mol)"), mgp = c(2, .8, 0))

ggplot(calcite) +
  geom_path(aes(x = SrCa_c * 1e3, y = MnCa_c * 1e3, color = time/3600)) +
  scale_color_distiller(palette = "RdBu") +
  geom_point(data = DB1, aes(x = SrCa, y = MnCa, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  theme(legend.title = element_text(size = 12, margin = margin(b = 15))) +
  labs(x = "Sr/Ca (mmol/mol)",
       y = "Mn/Ca (mmol/mol)",
       fill = "depth (cm)",
       color = "time (h)")

ggplot(calcite) +
  geom_path(aes(x = d13_c, y = d18_c, color = time)) +
  scale_color_distiller(palette = "RdBu") +
  geom_point(data = DB1, aes(x = d13, y = d18, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  labs(x = expression(delta^"13"*"C"[c]*" (\u2030, VPDB)"),
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"),
       fill = "depth (cm)",
       color = "time (s)") 

# Medium version ----
# no CO2 degassing, variable Mn reduction rate
source('model/box_model_ode.R')
time = seq(0, 1e4, 10)
# parms["F_in"] = 0 # 3e-5
# parms["k_degas"] = 0 # 2e-7
# parms["k_diss"] = 0
state["DIC_sw"] = 2.5e-3
state["Mn_sw"] = 0
parms["peak"] = 0.7
parms["r_rd"] = 2
parms["xmid"] = 0.8
parms["steep"] = 5
# parms["lamda"] = 3e-3
results = ode(y = state,
              times = time,
              parms = parms,
              func = soil_water_chemistry)
results = as.data.frame(results) |> filter(time > 100)
calcite = results |> filter(Jp > 0)
par(mfrow = c(3,2), mar = margin(3, 3, 1, 1))
plot(results$time, results$pH, type = "l", xlab = "time (s)", ylab = "pH", mgp = c(2, .8, 0))
plot(results$time, results$Alk, type = "l", xlab = "time (s)", ylab = "Alk (mol/L)", mgp = c(2, .8, 0))
plot(results$time, results$DIC_sw, type = "l", xlab = "time (s)", ylab = "DIC (mol/L)", mgp = c(2, .8, 0))
plot(calcite$time, 1e3*calcite$MnCa_c, type = "l", xlab = "time (s)", ylab = expression("Mn/Ca"[c]*" (mmol/mol)"), mgp = c(2, .8, 0))
plot(calcite$time, calcite$Jp, type = "l", xlab = "time (s)", ylab = "calcite (mol/s)", mgp = c(2, .8, 0))
plot(results$time, results$Mn_rd, type = "l", xlab = "time (s)", ylab = "Mn reduction rate (mol/s)", mgp = c(2, .8, 0))
# plot(calcite$time, 1e3*calcite$SrCa_c, type = "l", xlab = "time (s)", ylab = expression("Sr/Ca"[c]*" (mmol/mol)"), mgp = c(2, .8, 0))

ggplot(calcite) +
  geom_path(aes(x = SrCa_c * 1e3, y = MnCa_c * 1e3, color = time/3600)) +
  scale_color_distiller(palette = "RdBu") +
  geom_point(data = DB1, aes(x = SrCa, y = MnCa, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  theme(legend.title = element_text(size = 12, margin = margin(b = 15))) +
  labs(x = "Sr/Ca (mmol/mol)",
       y = "Mn/Ca (mmol/mol)",
       fill = "depth (cm)",
       color = "time (h)")

ggplot(calcite) +
  geom_path(aes(x = d13_c, y = d18_c, color = time)) +
  scale_color_distiller(palette = "RdBu") +
  geom_point(data = DB1, aes(x = d13, y = d18, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  labs(x = expression(delta^"13"*"C"[c]*" (\u2030, VPDB)"),
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"),
       fill = "depth (cm)",
       color = "time (s)") 

# No degassing version ----
source('model/box_model_ode.R')
time = seq(0, 1e4, 1)
# parms["F_in"] = 0 # 3e-5
parms["k_degas"] = 0 # 2e-7
# parms["k_diss"] = 1e-2
parms["kp"] = .1
state["DIC_sw"] = 2.6e-3
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
# plot(results$time, results$degas, type = "l", xlab = "time (s)", ylab = "Degass (mol/s)", mgp = c(2, .8, 0))
# plot(results$time, results$HCO3_rd, type = "l", xlab = "time (s)", ylab = "Mn_reduction (mol/s)", mgp = c(2, .8, 0))
# plot(results$time, results$CO2_diss, type = "l", xlab = "time (s)", ylab = "CO2 dissolution (mol/s)", mgp = c(2, .8, 0))
plot(calcite$time, 1e3*calcite$MnCa_c, type = "l", xlab = "time (s)", ylab = expression("Mn/Ca"[c]*" (mmol/mol)"), mgp = c(2, .8, 0))
plot(calcite$time, calcite$Jp, type = "l", xlab = "time (s)", ylab = "calcite (mol/s)", mgp = c(2, .8, 0))
# plot(results$time, results$Mn_rd, type = "l", xlab = "time (s)", ylab = "Mn reduction rate (mol/s)", mgp = c(2, .8, 0))
plot(calcite$time, 1e3*calcite$SrCa_c, type = "l", xlab = "time (s)", ylab = expression("Sr/Ca"[c]*" (mmol/mol)"), mgp = c(2, .8, 0))

p1 = ggplot(calcite) +
  geom_point(aes(x = SrCa_c * 1e3, y = MnCa_c * 1e3, color = time/3600), size = 1) +
  scale_color_distiller(palette = "RdBu") +
  geom_point(data = DB1, aes(x = SrCa, y = MnCa, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  theme(legend.title = element_text(size = 12, margin = margin(b = 15))) +
  labs(x = "Sr/Ca (mmol/mol)",
       y = "Mn/Ca (mmol/mol)",
       fill = "depth (cm)",
       color = "time (h)")

p2 = ggplot(calcite) +
  geom_point(aes(x = d13_c, y = d18_c, color = time/3600), size = 1) +
  scale_color_distiller(palette = "RdBu") +
  geom_point(data = DB1, aes(x = d13, y = d18, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  theme(legend.title = element_text(size = 12, margin = margin(b = 15))) +
  labs(x = expression(delta^"13"*"C"[c]*" (\u2030, VPDB)"),
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"),
       fill = "depth (cm)",
       color = "time (h)") 
p1 + p2 + plot_layout(guides = "collect") &
  theme(legend.position = "top")
