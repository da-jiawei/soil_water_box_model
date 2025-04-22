rm(list = ls())
library(tidyverse)
library(ggpubr)
library(readxl)
library(scales)
theme = theme(axis.text.x = element_text(margin = margin(t = 0.1, unit = "cm")),
              axis.text.y = element_text(margin = margin(r = 0.1, unit = "cm")),
              axis.ticks.length=unit(0.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              text = element_text(color = "black", size = 10),
              axis.title = element_text(size = 12), 
              axis.text = element_text(color = "black", size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())

# data grooming ---- 
DB = read_xlsx("data/DanceBayou_all_carb_data.xlsx")
DB1 = data.frame(DB$ID, DB$site, DB$depth, DB$d13C, DB$d18O, DB$`Sr, L1`, DB$`Mg, L1`, DB$`Mn, L1`)
names(DB1) = c("ID", "site", "depth", "d13", "d18", "SrCa_mass", "MgCa_mass", "MnCa_mass")
DB1 = DB1 %>%
  filter(site == "1") %>%
  mutate(MgCa = 1e3 * MgCa_mass * 40.078 / 24.305,
         SrCa = 1e3 * SrCa_mass * 40.078 / 87.62,
         MnCa = 1e3 * MnCa_mass * 40.078 / 54.94)
DB1$depth = as.numeric(DB1$depth)
source('box_model_Mn_redox.R')
cat("\014")

# fraction of remaining soil water at max respiration ----
vars = ctrl()
vars$k_degas = 0
vars$F_evap = 1e-3
# vars$F_in = 0
for (i in 1:5) {
  vars$xmax = 0.5 + 0.1 * i
  dat = SWTS_bm(vars) |>
    filter(Jp != 0)
  dat$xmax = vars$xmax
  if (i == 1) {
    sims = dat
  } else {
    sims = rbind(sims, dat)
  }
}

p1 = ggplot(sims) +
  geom_path(aes(x = SrCa_c * 1e3, y = MnCa_c * 1e3, color = xmax, group = xmax)) +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  geom_point(data = DB1, aes(x = SrCa, y = MnCa, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  # scale_y_continuous(limits = c(0, max(DB1$MnCa))) +
  theme_bw() + theme +
  guides(fill = "none") +
  labs(x = "Sr/Ca (mmol/mol)",
       y = "Mn/Ca (mmol/mol)",
       color = expression(italic(f)[sw]*" w/ max R"[s]))
p1

ggplot(sims) +
  geom_point(aes(x = SrCa_c * 1e3, y = MgCa_c * 1e3, color = xmax, group = xmax)) +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  geom_point(data = DB1, aes(x = SrCa, y = MgCa, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme

# mid-point of the sigmoidal curve ----
vars = ctrl()
vars$k_degas = 0
vars$F_evap = 1e-3
# vars$F_in = 0
for (i in 1:5) {
  vars$xmid = 0.9 + 0.02 * i
  dat = SWTS_bm(vars) |>
    filter(Jp != 0)
  dat$xmid = vars$xmid
  if (i == 1) {
    sims = dat
  } else {
    sims = rbind(sims, dat)
  }
}

p2 = ggplot(sims) +
  geom_path(aes(x = SrCa_c * 1e3, y = MnCa_c * 1e3, color = xmid, group = xmid)) +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  geom_point(data = DB1, aes(x = SrCa, y = MnCa, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  # scale_y_continuous(limits = c(0, max(DB1$MnCa))) +
  theme_bw() + theme +
  guides(fill = "none") +
  labs(x = "Sr/Ca (mmol/mol)",
       y = "Mn/Ca (mmol/mol)",
       color = expression("xmid"))
p2
# steepness of the sigmoidal curve ----
vars = ctrl()
vars$k_degas = 0
vars$F_evap = 1e-3
# vars$F_in = 0
for (i in 1:4) {
  vars$steep = 10 + 5*i
  dat = SWTS_bm(vars) |>
    filter(Jp != 0)
  dat$steep = vars$steep
  if (i == 1) {
    sims = dat
  } else {
    sims = rbind(sims, dat)
  }
}

p3 = ggplot(sims) +
  geom_path(aes(x = SrCa_c * 1e3, y = MnCa_c * 1e3, color = steep, group = steep)) +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  geom_point(data = DB1, aes(x = SrCa, y = MnCa, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  # scale_y_continuous(limits = c(0, max(DB1$MnCa))) +
  theme_bw() + theme +
  guides(fill = "none") +
  labs(x = "Sr/Ca (mmol/mol)",
       y = "Mn/Ca (mmol/mol)",
       color = expression("steep"))
p3

# respired CO2 ----
vars = ctrl()
vars$k_degas = 0
vars$F_evap = 1e-3
# vars$F_in = 0
for (i in 1:5) {
  vars$res_Co = i * 2e-6
  dat = SWTS_bm(vars) |>
    filter(Jp != 0)
  dat$res_Co = vars$res_Co
  if (i == 1) {
    sims = dat
  } else {
    sims = rbind(sims, dat)
  }
}

p4 = ggplot(sims) +
  geom_path(aes(x = SrCa_c * 1e3, y = MnCa_c * 1e3, color = log10(res_Co), group = res_Co)) +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  geom_point(data = DB1, aes(x = SrCa, y = MnCa, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  # scale_y_continuous(limits = c(0, max(DB1$MnCa))) +
  theme_bw() + theme +
  guides(fill = "none") +
  labs(x = "Sr/Ca (mmol/mol)",
       y = "Mn/Ca (mmol/mol)",
       color = expression("log"[10]*"(C"[res]*") (mol/s)"))
p4
ggplot(sims) +
  geom_path(aes(x = d13c, y = d18c, color = log10(res_Co), group = res_Co)) +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  geom_point(data = DB1, aes(x = d13, y = d18, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  # scale_x_continuous(limits = c(-11, -2)) +
  # scale_y_continuous(limits = c(-5, 0)) +
  guides(fill = "none") +
  labs(x = expression(delta^"13"*"C"[c]*" (\u2030, VPDB)"),
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"),
       fill = "depth (cm)",
       color = expression("log"[10]*"(C"[res]*") (mol/s)"))

ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, align = "hv")
ggsave("figures/sens_test_Mn_wo_degassing.jpg", width = 8.6, height = 6)


# turn off degassing
vars = ctrl()
# vars$k_degas = 0
# vars$F_evap = 1e-3
# vars$lamda = 0
vars$F_in = 0
dat = SWTS_bm(vars) |>
  filter(Jp != 0)

ggplot(dat) +
  geom_path(aes(x = time, y = Jp))

ggplot(dat) +
  geom_path(aes(x = SrCa_c * 1e3, y = MnCa_c * 1e3, color = fraction)) +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  geom_point(data = DB1, aes(x = SrCa, y = MnCa, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  # scale_y_continuous(limits = c(0, max(DB1$MnCa))) +
  theme_bw() + theme +
  guides(fill = "none") +
  labs(x = "Sr/Ca (mmol/mol)",
       y = "Mn/Ca (mmol/mol)")

ggplot(dat) +
  geom_path(aes(x = d13c, y = d18c, color = fraction)) +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  geom_point(data = DB1, aes(x = d13, y = d18, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  # scale_x_continuous(limits = c(-11, -2)) +
  # scale_y_continuous(limits = c(-5, 0)) +
  guides(fill = "none") +
  labs(x = expression(delta^"13"*"C"[c]*" (\u2030, VPDB)"),
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"),
       fill = "depth (cm)")
