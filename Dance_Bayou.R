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

# run the box model ----
source('box_model.R')
vars = ctrl()
# parameters that affect d13Cc
vars$F_in = 5e-3
vars$d13_DIC_p = -5
# parameters that affect d18Oc
vars$Tsoil = 18
dat = SWTS_bm(vars) %>%
  filter(Jp != 0)

## plot ----
# ggplot(dat, aes(x = time, y = Jp)) +
#   geom_line()

p1 = ggplot(dat) +
  geom_point(aes(x = d13c, y = d18c, color = fraction), 
             shape = 21, size = 3) +
  scale_color_distiller(palette = "RdBu", direction = 1, 
                        breaks = seq(min(dat$fraction), max(dat$fraction), length.out = 4),
                        labels = label_number(accuracy = 0.1)) +
  geom_point(data = DB1, aes(x = d13, y = d18, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  scale_x_continuous(limits = c(-11, -2)) +
  scale_y_continuous(limits = c(-5, 0)) +
  labs(x = expression(delta^"13"*"C"[c]*" (\u2030, VPDB)"),
       y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"),
       fill = "depth (cm)")
p1
p2 = ggplot(dat) +
  geom_point(aes(x = 1e3 * MgCa_c, y = 1e3 * SrCa_c, color = fraction),
             shape = 21, size = 3) +
  scale_color_distiller(palette = "RdBu", direction = 1, 
                        breaks = seq(min(dat$fraction), max(dat$fraction), length.out = 4),
                        labels = label_number(accuracy = 0.1)) +
  geom_point(data = DB1, aes(x = MgCa, y = SrCa, fill = depth), shape = 22, size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw() + theme +
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.3)) +
  labs(x = "Mg/Ca (mmol/mol)",
       y = "Sr/Ca (mmol/mol)",
       fill = "depth (cm)")
p2
# p3 = ggplot(dat) +
#   geom_point(aes(x = 1e3 * SrCa_c, y = d18c, color = fraction), 
#              shape = 21, size = 3) +
#   scale_color_distiller(palette = "RdBu", direction = 1, 
#                         breaks = seq(min(dat$fraction), max(dat$fraction), length.out = 4),
#                         labels = label_number(accuracy = 0.1)) +
#   geom_point(data = DB1, aes(x = SrCa, y = d18, fill = depth), shape = 22, size = 3) +
#   scale_fill_viridis_c(direction = -1) +
#   theme_bw() + theme +
#   labs(x = "Sr/Ca (mmol/mol)",
#        y = expression(delta^"18"*"O"[c]*" (\u2030, VPDB)"),
#        fill = "depth (cm)")
# p3
ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv", common.legend = TRUE)
ggsave("figures/degassing_evaporation.jpg", width = 9, height = 4)
