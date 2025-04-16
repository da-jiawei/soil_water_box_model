rm(list = ls())
# Distribution coefficients as a function of calcite precipitation rate - Lorens (1981) 

omega = seq(1, 1.5, 0.01) # saturation state
cal_rate = rep(0, length(omega))
for (i in 1:length(omega)) {
  if(omega[i] <= 1.5) {
    cal_rate[i] = 0.28 * omega[i] # calcite precipitation rate - nmol/mg/min
  } else {
    cal_rate[i] = 7.5 * omega[i] - 15 # omega < 5.5
  }
}

k_Sr = 10 ^ (0.249 * log10(cal_rate) - 1.57)
k_Mn = 10 ^ (-0.266 * log10(cal_rate) + 1.35)
