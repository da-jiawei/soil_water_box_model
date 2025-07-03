rm(list = ls())
#### Instructions ----
# CO2 + H2O = H2CO3      - kH
# H2CO3 = HCO3- + H+     - k1
# HCO3- = CO32- + H+     - k2

# input parameters ----
temp = 20
temp_kelvin = temp + 273.15
CO2_g = 400 # ppm

# mass action constants ----
kH = 10^-1.5 # mol/atm
k1 = 10^-6.3
k2 = 10^-10.3
kcc = 3.3e-9
kw = 1e-14

# open system ----
# DIC species in equilibrium with CO2
# kH = H2CO3 / CO2
# k1 = (H * HCO3) / H2CO3
# k2 = (H * CO3) / HCO3
# k1 * k2 = (H^2 * CO3) / H2CO3
# kw = H * OH
# charge balance
# H = HCO3 + 2*CO3 + OH
# H = HCO3
# HCO3 = sqrt(k1 * kH * CO2)

open_system = function(CO2_g){
  HCO3 = sqrt(k1 * kH * CO2_g * 1e-6)
  pH = -log10(HCO3)
  return(results = data.frame(pH = pH, HCO3 = HCO3))
}

open_system(1e2)

CO2_g = seq(1e3, 1e4, 1e3)
for (i in 1:length(CO2_g)) {
  result = open_system(CO2_g[i])
  result$CO2 = CO2_g[i]
  if (i == 1) {
    sims = result
  } else {
    sims = rbind(sims, result)
  }
}
plot(sims$CO2, sims$pH, type = "b")
plot(sims$CO2, sims$HCO3, type = "b")

# closed system ----
find_proton = function(proton, DIC, k1, k2, kw){
  # k1 = (proton * HCO3) / H2CO3
  # k2 = (proton * CO3) / HCO3
  # k1 * k2 = (proton^2 * CO3) / H2CO3
  # DIC = H2CO3 + HCO3- + CO32-
  # H2CO3 = DIC - (H2CO3 * k1) / proton - (H2CO3 * k1 * k2) / (proton^2)
  denom = 1 + k1/proton + (k1*k2)/(proton^2) 
  H2CO3 = DIC / denom
  # H+ = HCO3- + 2*CO32- + OH-
  neutral_charge = (H2CO3*k1)/proton + 2*(k1*k2*H2CO3)/(proton^2) + kw/proton - proton
  return(neutral_charge)
}

closed_system = function(CO2_g) {
  DIC = CO2_g * 1e-6 * kH
  proton = uniroot(find_proton, lower = 1e-14, upper = 1e-1, DIC = DIC, k1 = k1, k2 = k2, kw = kw)$root
  denom = 1 + k1/proton + (k1*k2)/(proton^2)
  H2CO3 = DIC / denom
  HCO3 = H2CO3 * k1 / proton
  CO3 = HCO3 * k2 / proton
  pH = -log10(proton)
  return(results = data.frame(pH = pH, H2CO3 = H2CO3, HCO3 = HCO3, CO3 = CO3))
}

closed_system(1e2)

CO2_g = seq(5e2, 5e4, 5e2)
for (i in 1:length(CO2_g)) {
  result = closed_system(CO2_g[i])
  result$CO2 = CO2_g[i]
  if (i == 1) {
    sims = result
  } else {
    sims = rbind(sims, result)
  }
}

plot(sims$CO2, sims$pH, type = "l")
plot(sims$CO2, sims$HCO3, type = "l")

# alkalinity ----
find_proton = function(proton, Alk, k1, k2, kw){
  # k1 = (proton * HCO3) / H2CO3
  # k2 = (proton * CO3) / HCO3
  # k1 * k2 = (proton^2 * CO3) / H2CO3
  # kw = proton * OH
  # Alk = HCO3 + 2*CO3
  denom = 1 + 2*k2/proton
  HCO3 = Alk / denom
  H2CO3 = (HCO3 * proton) / k1
  CO3 = (HCO3 * k2) / proton
  # H = HCO3 + 2CO3 + OH
  neutral_charge = HCO3 + 2*CO3 + (kw/proton) - proton
  return(neutral_charge)
}

Alk = seq(1e-3, 5e-2, 1e-3) # mmol/L
for (i in 1:length(Alk)) {
  proton = uniroot(find_proton, lower = 1e-14, upper = 1e-1, Alk = Alk[i], k1 = k1, k2 = k2, kw = kw)$root
  denom = 1 + 2*k2/proton
  HCO3 = Alk / denom
  H2CO3 = (HCO3 * proton) / k1
  CO3 = (HCO3 * k2) / proton
  pH = -log10(proton)
  results = data.frame(pH = pH, H2CO3 = H2CO3, HCO3 = HCO3, CO3 = CO3)
  if (i == 1) {
    sims = results
  } else {
    sims = rbind(sims, results)
  }
}


