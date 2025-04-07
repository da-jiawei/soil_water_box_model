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

holocene = read_xlsx("data/Holocene_nodules_CLP.xlsx")
ggplot(holocene, aes(x = d13, y = d18, fill = T47)) +
  geom_point(shape = 21, size = 3) +
  theme_bw() + theme

# box model ----
source("box_model.R")
# temp
vars = ctrl()
vars$d18p = -10
vars$d18so = -10
# vars$F_in = 
sims = list()
for (i in 1:4) {
  vars$Tsoil = 10 + 5 * i
  dat = SWTS_bm(vars)
  dat$Tsoil = vars$Tsoil
  sims[[i]] = dat
}
selected_cols = lapply(sims, function(df) df[, c("Tsoil", "fraction", "SrCa_c", "MgCa_c", "d13c", "d18c")])
dat = do.call(rbind, selected_cols) %>% drop_na()

ggplot(dat) +
  geom_path(aes(x = d13c, y = d18c, color = Tsoil, group = Tsoil)) +
  scale_color_distiller(palette = "RdBu", direction = -1) +
  geom_point(data = holocene, aes(x = d13, y = d18, fill = T47), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  # scale_x_continuous(limits = c(-11, -2)) +
  # scale_y_continuous(limits = c(-5, 0)) +
  labs(x = expression(delta^"13"*"C"[c]*" (\u2030, VPDB)"),
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"),
       fill = expression(paste("T"[47]*"(", degree, "C)")),
       color = expression(paste("T (", degree, "C)")))
