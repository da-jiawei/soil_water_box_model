rm(list = ls())
#### Instructions ----
# CO2 + H2O = H2CO3      - kH
# H2CO3 = HCO3- + H+     - k1
# HCO3- = CO32- + H+     - k2
# Ca2+ + CO32- = CaCO3   - kcc
# Ca2+ + 2HCO3- = CaCO3 + H2O + CO2 

# input parameters ----
temp = 20
temp_kelvin = temp + 273.15
CO2_g = 400 # ppm

# constants and coefficients ----
kH = 3.4e-2 # mol/atm
k1 = 10^((-356.3094 - 0.06091964*temp_kelvin + 21834.37/temp_kelvin + 126.8339*log10(temp_kelvin) - 1684915/(temp_kelvin^2)))
k2 = 10^((-107.8871 - 0.032528493*temp_kelvin + 5151.79/temp_kelvin + 38.92561*log10(temp_kelvin) - 563713.9/(temp_kelvin^2)))
kcc = 3.3e-9
kw = 1e-14

# species abundance ----
find_proton = function(proton, H2CO3, k1, k2, kw, kcc){
  # k1 = (proton * HCO3) / H2CO3
  # k2 = (proton * CO3) / HCO3
  # k1 * k2 = (proton^2 * CO3) / H2CO3
  # kw = proton * OH
  # kcc = Ca * CO3
  # assuming calcite dissolution in pure water
  # 2*mCa2+ + mH+ = mHCO3- + 2*mCO32- + mOH-
  term1 = proton + 2 * (kcc * proton^2) / (k1 * k2 * H2CO3) 
  term2 = (H2CO3 * k1 / proton) + 2 * (k1 * k2 * H2CO3) / (proton^2) + (kw / proton)
  return(term1 - term2)
}

open_system = function(CO2_g){
  H2CO3 = kH * CO2_g * 1e-6
  proton = uniroot(find_proton, lower = 10^-14, upper = 10^-1, H2CO3 = H2CO3, k1 = k1, k2 = k2, kw = kw, kcc = kcc)$root
  pH = -log10(proton)
  HCO3 = H2CO3 * k1 / proton
  CO3 = HCO3 * k2 / proton
  Ca = kcc / CO3
  return(results = data.frame(pH = pH, H2CO3 = H2CO3, HCO3 = HCO3, CO3 = CO3, Ca = Ca))
}

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


