rm(list = ls())
library(tidyverse)
library(ggpubr)
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
cat("\014")

# root finding function for [H+]
# calcite_eq = function(H, ALK, DIC){
#   term1 = (2 * k1 * k2) / (H^2) + (k1 / H)
#   term2 = 1 + (k1 * k2) / (H^2) + (k1 / H)
#   return(ALK * term2 - DIC * term1)
# }

solve_H = function(H, Alk, pCO2) {
  H2CO3 = kH * pCO2 * 1e-6 
  HCO3 = k1 * H2CO3 / H
  CO3 = k2 * HCO3 / H
  OH = kw / H
  Alk_calc = HCO3 + 2 * CO3 + OH 
  return(Alk_calc - Alk)
}

# the first and second disassociation constants for carbonic acid
kw = 10^-14
kH = 10^-1.5 # Henry's law constant
k1 = 10^-6.3
k2 = 10^-10.3

# test
dat = data.frame(matrix(nrow = 5, ncol = 8))
names(dat) = c("alk","pCO2", "H2CO3", "HCO3", "CO3", "OH", "pH", "alkalinity")
for (i in 1:5) {
  dat$alk[i] = 10^-i
  dat$pCO2[i] = 10^4
  dat$H2CO3[i] = kH * dat$pCO2[i] * 1e-6
  H = uniroot(solve_H, c(1e-14, 1e-1), Alk = dat$alk[i], pCO2 = dat$pCO2[i])$root
  dat$pH[i] = -log10(H)
  dat$HCO3[i] = k1 * dat$H2CO3[i] / H
  dat$CO3[i] = k2 * dat$HCO3[i] / H
  dat$OH[i] = kw / H
  dat$alkalinity[i] = dat$HCO3[i] + 2 * dat$CO3[i] + dat$OH[i]
}
dat



ksp = 3.3e-9 # solubility product of calcite
kd_Sr = 0.057 # partition coefficient for Sr
kd_Mg = 0.031 # partition coefficient for Mg

kp = 1e-3 # rate constant for carbonate precipitation (made up)

### inputs ----
Vo = 10 # initial pore water volume
P = 0 # rainfall
E = 0 # evaporation
# Q = P - E # outflow into deeper soil layers
Q = 0
DCO2 = 100 # degassing_rate: ppmv/time

# rainfall chemistry
Ca_p = 1e-4 # unit - mol/volume
Mg_p = 1e-5
Sr_p = 1e-6
DIC_p = 2.3e-4

# weathering input
Ca_w = 0 # unit - mol/time
Mg_w = 0
Sr_w = 0
DIC_w = 0

# parameterization
dt = 1  # time
time = 60 * 60
num = time / dt

V = rep(0, num)
pCO2 = rep(0, num)
Ca_s = rep(0, num) 
Mg_s = rep(0, num)
Sr_s = rep(0, num)
MgCa_s = rep(0, num)
SrCa_s = rep(0, num)
DIC_s = rep(0, num)
CO3_s = rep(0, num)
HCO3_s = rep(0, num)
H2CO3_s = rep(0, num)
H_s = rep(0, num)
pH = rep(0, num)
Jp = rep(0, num)
SrCa_c = rep(0, num)
MgCa_c = rep(0, num)

# initial soil solution
V[1] = Vo
Ca_s[1] = 1e-2 # unit - mol/volume
Mg_s[1] = 1e-3
Sr_s[1] = 1e-4
SrCa_s[1] = Sr_s[1] / Ca_s[1]
MgCa_s[1] = Mg_s[1] / Ca_s[1]


# loop ----
for (i in 1:(num-1)) {
  Alk = 2 * (Ca_s[i] + Mg_s[i] + Sr_s[i])
  H2CO3_s[i] = pCO2[i] * 1e-6 * kH
  results = uniroot(alk_co2, c(1e-14, 1e-1), ALK = Alk, H2CO3 = H2CO3_s[i])
  H_s[i] = results$root
  pH[i] = -log10(H_s[i])
  HCO3_s[i] = k1 * H2CO3_s[i] / H_s[i]
  CO3_s[i] = k2 * HCO3_s[i] / H_s[i]
  DIC_s[i] = H2CO3_s[i] + HCO3_s[i] + CO3_s[i]
  omega = Ca_s[i] * CO3_s[i] / ksp # the saturation state of calcite
  if(omega > 1) {
    Jp[i] = kp * (omega - 1) # precipitation flux of calcite - unit: mol/time
  } else {
    Jp[i] = 0
  }
  dV = dt * (P - E - Q) # unit: L/time
  # DIC is not driven by precipitation, but pCO2
  # dDIC = (dt / V[i]) * (P * (DIC_p - DIC_s[i]) + E * DIC_s[i] + DIC_w - Jp[i])
  dCa = (dt / V[i]) * (P * (Ca_p - Ca_s[i]) + E * Ca_s[i] + Ca_w - Jp[i])
  dMg = (dt / V[i]) * (P * (Mg_p - Mg_s[i]) + E * Mg_s[i] + Mg_w - (Jp[i] * kd_Mg * Mg_s[i] / Ca_s[i]))
  dSr = (dt / V[i]) * (P * (Sr_p - Sr_s[i]) + E * Sr_s[i] + Sr_w - (Jp[i] * kd_Sr * Sr_s[i] / Ca_s[i]))
  if(Jp[i] == 0) {
    SrCa_c[i] = NA
    MgCa_c[i] = NA
  } else {
    SrCa_c[i] = kd_Mg * Mg_s[i] / Ca_s[i]
    MgCa_c[i] = kd_Sr * Sr_s[i] / Ca_s[i]
  }
  V[i+1] = V[i] + dV
  if(V[i+1] <= 0) {
    print("V is negative, stopping loop.")
    break
  } 
  pCO2[i+1] = pCO2[i] - DCO2
  if(pCO2[i+1] <= 400) {
    pCO2[i+1] = 400
  }
  
  # DIC_s[i+1] = DIC_s[i] + dDIC
  Ca_s[i+1] = Ca_s[i] + dCa
  Mg_s[i+1] = Mg_s[i] + dMg
  Sr_s[i+1] = Sr_s[i] + dSr
  SrCa_s[i+1] = Sr_s[i+1] / Ca_s[i+1]
  MgCa_s[i+1] = Mg_s[i+1] / Ca_s[i+1]
}
dat = data.frame(time = seq(dt, time, dt), V = V, Jp = Jp, pCO2 = pCO2,
                 DIC = DIC_s, H2CO3_s = H2CO3_s, HCO3_s = HCO3_s, CO3_s = CO3_s,
                 Ca = Ca_s, Mg = Mg_s, Sr = Sr_s, H = H_s, pH = pH,
                 MgCa_s = MgCa_s, SrCa_s = SrCa_s, MgCa_c = MgCa_c, SrCa_c = SrCa_c)
dat = dat[seq(1,i),]