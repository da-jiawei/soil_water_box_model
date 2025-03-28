# input
Tair = 10
d13C_co2 = -8
pH = 7
DIC_p = 2.2e-4

# coefficients
k1 = 10^-6.3
k2 = 10^-10.3

# carbon isotope 
Tair_K = Tair + 273.15
alpha13_HCO3_CO3 = exp((-0.87 * 1e3 / Tair_K + 2.52) * 1e-3) # Mook et al. (1974)
alpha13_H2CO3_CO2 = exp((6e3 / Tair_K^2 - 0.9) * 1e-3) # Deines et al. (1974)
alpha13_HCO3_CO2 = exp((1.1e6 / Tair_K^2 - 4.5) * 1e-3) # Deines et al. (1974)

d13C_H2CO3 = (d13C_co2 + 1000) * alpha13_H2CO3_CO2 - 1000
d13C_HCO3 = (d13C_co2 + 1000) * alpha13_HCO3_CO2 - 1000
d13C_CO3 = (d13C_HCO3 + 1000) / alpha13_HCO3_CO3 - 1000

# calculating the abundance of DIC species
H_p = 10 ^ (-pH)
CO3_p = DIC_p / ((H_p^2) / (k1 * k2) + H_p / k2 + 1) # carbonate ion
HCO3_p = H_p * CO3_p / k2
H2CO3_p = H_p * HCO3_p / k1 
X_CO3 = CO3_p / DIC_p
X_HCO3 = HCO3_p / DIC_p
X_H2CO3 = H2CO3_p / DIC_p

# mass balance
d13C_DIC = d13C_H2CO3 * X_H2CO3 + d13C_HCO3 * X_HCO3 + d13C_CO3 * d13C_CO3
