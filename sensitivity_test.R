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
source('box_model.R')
cat("\014")

# calcite precipitation rate ----
vars = ctrl()
sims = list()
for (i in 1:5) {
  vars$kp = 2e-3 * i
  dat = SWTS_bm(vars)
  dat$kp = vars$kp
  sims[[i]] = dat
}
selected_cols = lapply(sims, function(df) df[, c("kp", "fraction", "SrCa_c", "MgCa_c", "d13c", "d18c")])
dat = do.call(rbind, selected_cols) %>% drop_na()
ggplot(dat, aes(x = 1e3 * SrCa_c, y = d18c, color = kp, group = kp)) +
  geom_point() +
  theme_bw() + theme +
  labs(x = "Sr/Ca (mmol/mol)",
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"))

ggplot(dat, aes(x = d13c, y = d18c, color = kp, group = kp)) +
  geom_point() +
  theme_bw() + theme +
  labs(x = expression(delta^"13"*"C"[c]*" (\u2030, VPDB)"),
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"))

ggplot(dat, aes(x = 1e3 * SrCa_c, y = MgCa_c, color = kp, group = kp)) +
  geom_point() +
  theme_bw() + theme +
  labs(x = "Sr/Ca (mmol/mol)",
       y = "Mg/Ca (mol/mol)")

# CO2 degassing rate ----
vars = ctrl()
sims = list()
for (i in 1:4) {
  vars$k_degas = 10^(-i-5)
  dat = SWTS_bm(vars)
  dat$kd = vars$k_degas
  sims[[i]] = dat
}
selected_cols = lapply(sims, function(df) df[, c("kd", "fraction", "SrCa_c", "MgCa_c", "d13c", "d18c")])
dat = do.call(rbind, selected_cols) %>% drop_na()
dat = dat %>% filter(SrCa_c < 0.1)
ggplot(dat, aes(x = 1e3 * SrCa_c, y = d18c, color = -log10(kd), group = kd)) +
  geom_point() +
  theme_bw() + theme +
  labs(x = "Sr/Ca (mmol/mol)",
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"))

ggplot(dat, aes(x = d13c, y = d18c, color = -log10(kd), group = kd)) +
  geom_path() +
  theme_bw() + theme +
  labs(x = expression(delta^"13"*"C"[c]*" (\u2030, VPDB)"),
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"))

ggplot(dat, aes(x = 1e3 * SrCa_c, y = MgCa_c, color = -log10(kd), group = kd)) +
  geom_point() +
  theme_bw() + theme +
  labs(x = "Sr/Ca (mmol/mol)",
       y = "Mg/Ca (mol/mol)")

# soil temperature ----
vars = ctrl()
sims = list()
for (i in 1:5) {
  vars$Tsoil = 5 * i
  dat = SWTS_bm(vars)
  dat$Tsoil = vars$Tsoil
  sims[[i]] = dat
}
selected_cols = lapply(sims, function(df) df[, c("Tsoil", "fraction", "SrCa_c", "MgCa_c", "d13c", "d18c")])
dat = do.call(rbind, selected_cols) %>% drop_na()
dat = dat %>% filter(SrCa_c < 1e-3)
p1 = ggplot(dat, aes(x = 1e3 * SrCa_c, y = d18c, color = Tsoil, group = Tsoil)) +
  geom_path() +
  theme_bw() + theme +
  labs(x = "Sr/Ca (mmol/mol)",
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"))

p2 = ggplot(dat, aes(x = d13c, y = d18c, color = Tsoil, group = Tsoil)) +
  geom_path() +
  theme_bw() + theme +
  labs(x = expression(delta^"13"*"C"[c]*" (\u2030, VPDB)"),
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"))

ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv", common.legend = TRUE)

ggplot(dat, aes(x = 1e3 * SrCa_c, y = MgCa_c, color = Tsoil, group = Tsoil)) +
  geom_point() +
  theme_bw() + theme +
  labs(x = "Sr/Ca (mmol/mol)",
       y = "Mg/Ca (mol/mol)")

