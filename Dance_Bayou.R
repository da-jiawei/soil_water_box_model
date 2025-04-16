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
DB1 = data.frame(DB$ID, DB$site, DB$depth, DB$d13C, DB$d18O, DB$`Sr, L1`, DB$`Mg, L1`)
names(DB1) = c("ID", "site", "depth", "d13", "d18", "SrCa_mass", "MgCa_mass")
DB1 = DB1 %>%
  filter(site == "1") %>%
  mutate(MgCa = 1000 * MgCa_mass * 40.078 / 24.305,
         SrCa = 1000 * SrCa_mass * 40.078 / 87.62)
DB1$depth = as.numeric(DB1$depth)
cat("\014")
source('box_model_Mn_redox.R')

# run the box model ----
vars = ctrl()
# vars$res_Co = 1e-6
# vars$d13_r = -14
# vars$Tsoil = 20
# vars$k_degas = 0
# vars$F_evap = 5e-1

dat = SWTS_bm(vars) %>%
  filter(Jp != 0)
# dat = dat %>%
#   filter(d18c <= max(DB1$d18))

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


# respired CO2 input ----
vars = ctrl()
sims = list()
for (i in 1:5) {
  vars$res_Co = i * 2e-5
  dat = SWTS_bm(vars)
  dat$res_Co = vars$res_Co
  sims[[i]] = dat
}
selected_cols = lapply(sims, function(df) df[, c("res_Co", "fraction", "SrCa_c", "MgCa_c", "d13c", "d18c")])
dat = do.call(rbind, selected_cols) %>% drop_na()

p1 = ggplot(dat) +
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
p1 

# d13C of respired CO2 ----
vars = ctrl()
sims = list()
for (i in 1:5) {
  vars$d13_r = -24 + 2*i
  dat = SWTS_bm(vars)
  dat$d13_r = vars$d13_r
  sims[[i]] = dat
}
selected_cols = lapply(sims, function(df) df[, c("d13_r", "fraction", "SrCa_c", "MgCa_c", "d13c", "d18c")])
dat = do.call(rbind, selected_cols) %>% drop_na()

p2 = ggplot(dat) +
  geom_path(aes(x = d13c, y = d18c, color = d13_r, group = d13_r)) +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  geom_point(data = DB1, aes(x = d13, y = d18, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  scale_x_continuous(limits = c(-11, -2)) +
  scale_y_continuous(limits = c(-5, 0)) +
  guides(fill = "none") +
  labs(x = expression(delta^"13"*"C"[c]*" (\u2030, VPDB)"),
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"),
       fill = "depth (cm)",
       color = expression(delta^"13"*"C"[res]*"(\u2030)"))
p2

# degassing rate ----
vars = ctrl()
sims = list()
for (i in 1:5) {
  vars$k_degas = i * 2e-8
  dat = SWTS_bm(vars)
  dat$k_degas = vars$k_degas
  sims[[i]] = dat
}
selected_cols = lapply(sims, function(df) df[, c("k_degas", "fraction", "SrCa_c", "MgCa_c", "d13c", "d18c")])
dat = do.call(rbind, selected_cols) %>% drop_na()

ggplot(dat) +
  geom_path(aes(x = d13c, y = d18c, color = log10(k_degas), group = k_degas)) +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  geom_point(data = DB1, aes(x = d13, y = d18, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  scale_x_continuous(limits = c(-11, -2)) +
  scale_y_continuous(limits = c(-5, 0)) +
  labs(x = expression(delta^"13"*"C"[c]*" (\u2030, VPDB)"),
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"),
       fill = "depth (cm)",
       color = expression("log"[10]*"(k"[degassing]*") (mol/s)"))

# calcite precipitation rate ----
vars = ctrl()
sims = list()
for (i in 1:5) {
  vars$kp = 10 ^ - (1 + i)
  dat = SWTS_bm(vars)
  dat$kp = vars$kp
  sims[[i]] = dat
}
selected_cols = lapply(sims, function(df) df[, c("kp", "fraction", "SrCa_c", "MgCa_c", "d13c", "d18c")])
dat = do.call(rbind, selected_cols) %>% drop_na()

ggplot(dat) +
  geom_path(aes(x = d13c, y = d18c, color = log10(kp), group = kp)) +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  geom_point(data = DB1, aes(x = d13, y = d18, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  scale_x_continuous(limits = c(-11, -2)) +
  scale_y_continuous(limits = c(-5, 0)) +
  labs(x = expression(delta^"13"*"C"[c]*" (\u2030, VPDB)"),
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"),
       fill = "depth (cm)",
       color = expression("log"[10]*"(k"[p]*") (mol/s)"))

# evaporation rate ----
vars = ctrl()
sims = list()
for (i in 1:4) {
  vars$F_evap = i * 2e-3
  dat = SWTS_bm(vars)
  dat$F_evap = vars$F_evap
  sims[[i]] = dat
}
selected_cols = lapply(sims, function(df) df[, c("F_evap", "fraction", "SrCa_c", "MgCa_c", "d13c", "d18c")])
dat = do.call(rbind, selected_cols) %>% drop_na()

p3 = ggplot(dat) +
  geom_path(aes(x = d13c, y = d18c, color = log10(F_evap), group = F_evap)) +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  geom_point(data = DB1, aes(x = d13, y = d18, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  scale_x_continuous(limits = c(-11, -2)) +
  scale_y_continuous(limits = c(-5, 0)) +
  guides(fill = "none") +
  labs(x = expression(delta^"13"*"C"[c]*" (\u2030, VPDB)"),
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"),
       fill = "depth (cm)",
       color = expression("log"[10]*"(F"[evap]*") (mol/s)"))
p3

# rainfall influx ----
vars = ctrl()
sims = list()
for (i in 1:4) {
  vars$F_in = i * 2e-4
  dat = SWTS_bm(vars)
  dat$F_in = vars$F_in
  sims[[i]] = dat
}
selected_cols = lapply(sims, function(df) df[, c("F_in", "fraction", "SrCa_c", "MgCa_c", "d13c", "d18c")])
dat = do.call(rbind, selected_cols) %>% drop_na()

p4 = ggplot(dat) +
  geom_path(aes(x = d13c, y = d18c, color = log10(F_in), group = F_in)) +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  geom_point(data = DB1, aes(x = d13, y = d18, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  scale_x_continuous(limits = c(-11, -2)) +
  scale_y_continuous(limits = c(-5, 0)) +
  guides(fill = "none") +
  labs(x = expression(delta^"13"*"C"[c]*" (\u2030, VPDB)"),
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"),
       fill = "depth (cm)",
       color = expression("log"[10]*"(F"["in"]*") (mol/L/s)"))
p4 

# d13C of influx-DIC ----
vars = ctrl()
sims = list()
for (i in 1:4) {
  vars$d13_DIC_p = -5 + i
  dat = SWTS_bm(vars)
  dat$d13_DIC_p = vars$d13_DIC_p
  sims[[i]] = dat
}
selected_cols = lapply(sims, function(df) df[, c("d13_DIC_p", "fraction", "SrCa_c", "MgCa_c", "d13c", "d18c")])
dat = do.call(rbind, selected_cols) %>% drop_na()

p5 = ggplot(dat) +
  geom_path(aes(x = d13c, y = d18c, color = d13_DIC_p, group = d13_DIC_p)) +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  geom_point(data = DB1, aes(x = d13, y = d18, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  scale_x_continuous(limits = c(-11, -2)) +
  scale_y_continuous(limits = c(-5, 0)) +
  guides(fill = "none") +
  labs(x = expression(delta^"13"*"C"[c]*" (\u2030, VPDB)"),
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"),
       fill = "depth (cm)",
       color = expression(delta^"13"*"C"[DIC]*" (\u2030)"))
p5 

# temp ----
vars = ctrl()
sims = list()
for (i in 1:4) {
  vars$Tsoil = 12 + 2 * i
  dat = SWTS_bm(vars)
  dat$Tsoil = vars$Tsoil
  sims[[i]] = dat
}
selected_cols = lapply(sims, function(df) df[, c("Tsoil", "fraction", "SrCa_c", "MgCa_c", "d13c", "d18c")])
dat = do.call(rbind, selected_cols) %>% drop_na()

p6 = ggplot(dat) +
  geom_path(aes(x = d13c, y = d18c, color = Tsoil, group = Tsoil)) +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  geom_point(data = DB1, aes(x = d13, y = d18, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  scale_x_continuous(limits = c(-11, -2)) +
  scale_y_continuous(limits = c(-5, 0)) +
  guides(fill = "none") +
  labs(x = expression(delta^"13"*"C"[c]*" (\u2030, VPDB)"),
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"),
       fill = "depth (cm)",
       color = expression(paste("T (", degree, "C)")))
p6

ggarrange(p1, p2, p3, p4, p5, p6, nrow = 3, ncol = 2, align = "hv")
ggsave("figures/sens_test_isotopes.jpg", width = 8, height = 7.2)


## plot ----
p1 = ggplot(dat) +
  geom_point(aes(x = d13c, y = d18c, color = fraction), 
             shape = 21, size = 3) +
  scale_color_distiller(palette = "RdBu", direction = 1, 
                        breaks = seq(min(dat$fraction), max(dat$fraction), length.out = 4),
                        labels = label_number(accuracy = 0.01)) +
  geom_point(data = DB1, aes(x = d13, y = d18, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  # scale_x_continuous(limits = c(-11, -2)) +
  # scale_y_continuous(limits = c(-5, 0)) +
  labs(x = expression(delta^"13"*"C"[c]*" (\u2030, VPDB)"),
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"),
       fill = "depth (cm)")
p1
p2 = ggplot(dat) +
  geom_point(aes(x = 1e3 * MgCa_c, y = 1e3 * SrCa_c, color = fraction),
             shape = 21, size = 3) +
  scale_color_distiller(palette = "RdBu", direction = 1, 
                        breaks = seq(min(dat$fraction), max(dat$fraction), length.out = 4),
                        labels = label_number(accuracy = 0.01)) +
  geom_point(data = DB1, aes(x = MgCa, y = SrCa, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.3)) +
  labs(x = "Mg/Ca (mmol/mol)",
       y = "Sr/Ca (mmol/mol)",
       fill = "depth (cm)")
p2
ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv", common.legend = TRUE)
ggsave("figures/degassing_evaporation.jpg", width = 9, height = 4)
