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

lime = read_xlsx("data/lime_mortar.xlsx")
lime = lime[1:56, c(2:4, 6, 7)]
names(lime) = c("ID", "d13", "d18", "SrCa", "MgCa")

# model ----
source('box_model.R')
vars = ctrl()
vars$d18so = -15
vars$d13_co2_initial = -8
dat = SWTS_bm(vars) %>%
  filter(Jp != 0) 

# plot ----
ggplot(dat, aes(x = d13c, y = d18c)) +
  geom_point(shape = 21, size = 2) +
  geom_point(data = lime, aes(x = d13, y = d18, color = ID))

ggplot(lime, aes(x = SrCa, y = MgCa, color = ID)) +
  geom_point()

