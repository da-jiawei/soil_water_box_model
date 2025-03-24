library(tidyverse)
library(ggpubr)
library(readxl)
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

DB = read_xlsx("data/DanceBayou_all_carb_data.xlsx")
DB1 = data.frame(DB$ID, DB$site, DB$d18O, DB$`Sr, L1`, DB$`Mg, L1`)
names(DB1) = c("ID", "site", "d18", "SrCa_mass", "MgCa_mass")
DB1 = DB1 %>%
  mutate(MgCa = 1000 * MgCa_mass * 40.078 / 24.305,
         SrCa = 1000 * SrCa_mass * 40.078 / 87.62)

# root finding function for [H+]
calcite_eq = function(H, ALK, DIC) {
  term1 = (2 * k1 * k2) / (H^2) + (k1 / H)
  term2 = 1 + (k1 * k2) / (H^2) + (k1 / H)
  return(ALK * term2 - DIC * term1)
}

## constants ----
# the first and second disassociation constants for carbonic acid
k1 = 10^-6.3
k2 = 10^-10.3

ksp = 3.3e-9 # solubility product of calcite
kd_Sr = 0.057 # partition coefficient for Sr
kd_Mg = 0.031 # partition coefficient for Mg

kp = 1e-3 # rate constant for carbonate precipitation (made up) - mol/L/s

# isotope standards
R18smow = 0.0020052
R18vpdb = 0.0020672

# molarity of water
H2O_mol = 55.51 # mol/L

## inputs ----
# isotope model
d18p = -3
RH = 0.7
Tsoil = 20
w = 0.5
Tsoil.K = Tsoil + 273.15
alpha18_c_w = exp((1.61e4 / Tsoil.K - 24.6) / 1000) # Tremaine (2011)
alpha18_l_v = exp(-2.0667 * 10^-3 - 0.4156 / Tsoil.K + (1.137 * 10 ^ 3) / (Tsoil.K ^ 2)) # Majoube (1971)
alpha18_diff = 1.028489
alpha18_diff = alpha18_diff * w + (1-w)
R18p = (d18p / 1000 + 1) * R18smow
# assuming atmospheric vapor in equilibrium with rainfall under soil temperature
R18v = R18p / alpha18_l_v

P = 0 # rainfall - L/time
E = 0.1 # evaporation - L/time
# Q = P - E # outflow into deeper soil layers
Q = 0

# rainfall chemistry
Ca_p = 1e-4 # unit - mol/L
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

# setting the initial soil solution at the same time
V = rep(10, num) # unit - L
d18_s = rep(d18p, num)
d18_c = rep(0, num)
Ca_s = rep(1e-3, num) # unit - mol/L
Mg_s = rep(1e-4, num)
Sr_s = rep(2e-6, num)
MgCa_s = rep(Mg_s[1] / Ca_s[1], num)
SrCa_s = rep(Sr_s[1] / Ca_s[1], num)
DIC_s = rep(2.3e-3, num)
CO3_s = rep(0, num)
HCO3_s = rep(0, num)
H2CO3_s = rep(0, num)
H_s = rep(0, num)
pH = rep(0, num)
Jp = rep(0, num)
SrCa_c = rep(0, num)
MgCa_c = rep(0, num)

## loop ----
for (i in 1:(num-1)) {
  R18s = (d18_s[i]/1000 + 1) * R18smow
  R18e = (R18s - RH * R18v * alpha18_l_v) / ((1 - RH) * alpha18_l_v * alpha18_diff) # Craig-Gordon model, assuming zero turbulence
  d18e = (R18e / R18smow - 1) * 1000
  R18c = R18s * alpha18_c_w
  d18_c[i] = (R18c / R18vpdb - 1) * 1000
  Alk = 2 * (Ca_s[i] + Mg_s[i] + Sr_s[i]) # calculating alkalinity
  results = uniroot(calcite_eq, c(1e-14, 1e-1), ALK = Alk, DIC = DIC_s[i], tol = 1e-14)
  H_s[i] = results$root
  pH[i] = -log10(H_s[i])
  CO3_s[i] = DIC_s[i] / ((H_s[i]^2) / (k1 * k2) + H_s[i] / k2 + 1) # carbonate ion
  HCO3_s[i] = H_s[i] * CO3_s[i] / k2
  H2CO3_s[i] = H_s[i] * HCO3_s[i] / k1
  omega = Ca_s[i] * CO3_s[i] / ksp # the saturation state of calcite
  if(omega > 1) {
    Jp[i] = kp * (omega - 1) # precipitation flux of calcite - unit: mol/L
  } else {
    Jp[i] = 0
  }
  dV = dt * (P - E - Q) # unit: L/time
  dDIC = (dt / V[i]) * (P * (DIC_p - DIC_s[i]) + E * DIC_s[i] + DIC_w - Jp[i])
  dCa = (dt / V[i]) * (P * (Ca_p - Ca_s[i]) + E * Ca_s[i] + Ca_w - Jp[i])
  dMg = (dt / V[i]) * (P * (Mg_p - Mg_s[i]) + E * Mg_s[i] + Mg_w - (Jp[i] * kd_Mg * Mg_s[i] / Ca_s[i]))
  dSr = (dt / V[i]) * (P * (Sr_p - Sr_s[i]) + E * Sr_s[i] + Sr_w - (Jp[i] * kd_Sr * Sr_s[i] / Ca_s[i]))
  dd18_s = dt * ((d18p - d18_s[i]) * (P / V[i]) - (d18e - d18_s[i]) * (E / V[i]) - (d18_c[i] - d18_s[i]) * (Jp[i] / (V[i] * H2O_mol)))
  if(Jp[i] == 0) {
    SrCa_c[i] = NA
    MgCa_c[i] = NA
  } else {
    SrCa_c[i] = kd_Mg * SrCa_s[i]
    MgCa_c[i] = kd_Sr * MgCa_s[i]
  }
  V[i+1] = V[i] + dV 
  if(V[i+1] <= 0) {
    print("V is negative, stopping loop.")
    break
  } 
  DIC_s[i+1] = DIC_s[i] + dDIC
  Ca_s[i+1] = Ca_s[i] + dCa
  Mg_s[i+1] = Mg_s[i] + dMg
  Sr_s[i+1] = Sr_s[i] + dSr
  SrCa_s[i+1] = Sr_s[i+1] / Ca_s[i+1]
  MgCa_s[i+1] = Mg_s[i+1] / Ca_s[i+1]
  d18_s[i+1] = d18_s[i] + dd18_s
}

dat = data.frame(time = seq(dt, time, dt), V = V, Jp = Jp,
                 DIC = DIC_s, H2CO3_s = H2CO3_s, HCO3_s = HCO3_s, CO3_s = CO3_s,
                 Ca = Ca_s, Mg = Mg_s, Sr = Sr_s, H = H_s, pH = pH,
                 MgCa_s = MgCa_s, SrCa_s = SrCa_s, MgCa_c = MgCa_c, SrCa_c = SrCa_c,
                 d18s = d18_s, d18c = d18_c)
dat = dat[seq(1,i),] %>% drop_na()


## plot ----
ggplot(dat) +
  geom_point(aes(x = time, y = Jp))

ggplot(dat) +
  geom_line(aes(x = V/V[1], y = SrCa_s), color = "red") +
  annotate("text", x = 0.2, y = 0.12, label = "soil water", color = "red") +
  geom_line(aes(x = V/V[1], y = SrCa_c), color = "blue") +
  annotate("text", x = 0.3, y = 0.05, label = "pedogenic carbonate", color = "blue") +
  theme_bw() + theme +
  scale_x_reverse(limits = c(1,0)) +
  labs(x = "fraction of remaining soil water",
       y = "Sr/Ca (mol/mol)")

dat$fraction = dat$V / dat$V[1]
p1 = ggplot(dat) +
  geom_point(aes(x = d18c, y = 1e3 * SrCa_c, fill = fraction),
             shape = 21, size = 3) +
  scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(0, 1)) +
  geom_point(data = DB1, aes(x = d18, y = SrCa), shape = 22, size = 3, color = "darkseagreen") +
  theme_bw() + theme +
  labs(x = expression(delta^"18"*"O"),
       y = "Sr/Ca (mmol/mol)") +
  scale_y_log10()
p1
p2 = ggplot(dat) +
  geom_point(aes(x = 1e3 * MgCa_c, y = 1e3 * SrCa_c, fill = fraction),
             shape = 21, size = 3) +
  scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(0, 1)) +
  geom_point(data = DB1, aes(x = MgCa, y = SrCa), shape = 22, size = 3, color = "darkseagreen") +
  theme_bw() + theme +
  labs(x = "Mg/Ca (mmol/mol)",
       y = "Sr/Ca (mmol/mol)") +
  scale_x_continuous(limits = c(0, 25)) +
  scale_y_continuous(limits = c(0, 0.3))
p2
ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv", common.legend = TRUE)
ggsave("figures/evaporation_model_0.1.jpg", width = 6.2, height = 4)

# DIC species
par(mar = c(4, 4, 1, 4))
plot(0, 0, xlim = c(5.5, 8), ylim = c(0, 3), axes = FALSE,
     xlab = "", ylab = "")

yext = range(log10(dat$H2CO3_s))
tix = seq(floor(min(yext)), 
          ceiling(max(yext)), by = 1)
H2CO3.rs = cbind(dat$pH,
               2 + (log10(dat$H2CO3_s) - min(tix)) / diff(range(tix)))
lines(lowess(H2CO3.rs[, 1], H2CO3.rs[, 2]))
axis(2, 2 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("log"[10]*"[H"[2]*"CO"[3]*"]"), 2, line = 2.5, at = 2.5)

yext = range(log10(dat$HCO3_s))
tix = seq(floor(min(yext)), 
          ceiling(max(yext)), by = 1)
HCO3.rs = cbind(dat$pH,
                 1 + (log10(dat$HCO3_s) - min(tix)) / diff(range(tix)))
lines(lowess(HCO3.rs[, 1], HCO3.rs[, 2]))
axis(4, 1 + (tix - min(tix)) / diff(range(tix)), tix)

mtext(expression("log"[10]*"[HCO"[3]^"-"*"]"), 4, line = 2.5, at = 1.5)

yext = range(log10(dat$CO3_s))
tix = seq(floor(min(yext)), 
          ceiling(max(yext)), by = 1)
CO3.rs = cbind(dat$pH,
                0 + (log10(dat$CO3_s) - min(tix)) / diff(range(tix)))
lines(lowess(CO3.rs[, 1], CO3.rs[, 2]))
axis(2, 0 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression("log"[10]*"[HCO"[3]*" "^"2-"*"]"), 2, line = 2.5, at = 0.5)

axis(1)
mtext("pH", 1, line = 2)
dev.off()
