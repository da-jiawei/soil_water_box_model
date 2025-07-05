rm(list = ls())
k1 = 10^-6.3
k2 = 10^-10.3

pH_finder = function(pH, Alk, DIC, k1, k2){
  # Alk = HCO3 + 2CO3 = HCO3 + 2*HCO3*k2/H
  # HCO3 = Alk / (1 + 2*k2/H)
  # DIC = H2CO3 + HCO3 + CO3 = H*HCO3/k1 + HCO3 + HCO3*k2/H
  # HCO3 = DIC / (H/k1 + 1 + k2/H)
  H = 10^(-pH)
  term1 = Alk / (1 + 2*k2/H)
  term2 = DIC / (H/k1 + 1 + k2/H)
  return(term1 - term2)
}

pH_finder2 = function(pH, Ca, Mg, Mn, Sr, DIC, k1, k2){
  Alk = 2 * (Ca + Mg + Mn + Sr)
  H = 10^(-pH)
  term1 = Alk / (1 + 2*k2/H)
  term2 = DIC / (H/k1 + 1 + k2/H)
  return(term1 - term2)
}

# sensitivity tests ----
# DIC
for (i in 1:1e2) {
  DIC = 1.1e-3 + 1e-4*i
  pH = uniroot(pH_finder, lower = 1, upper = 14, Alk = 2e-3, DIC = DIC, k1 = k1, k2 = k2)$root
  result = data.frame(DIC = DIC, pH = pH)
  if (i == 1) {
    sims = result
  } else {
    sims = rbind(sims, result)
  }
}
plot(sims$DIC * 1e3, sims$pH, type = "l",
     xlab = "DIC (mmol/L)", ylab = "pH")
text(8, 10, label = "Alk = 2 mmol/L", font = 2)

# Alkalinity
uniroot(pH_finder, lower = 1, upper = 14, Alk = 1e-3, DIC = 5e-3, k1 = k1, k2 = k2)$root
for (i in 1:80) {
  Alk = 1e-3 + 1e-4*i
  pH = uniroot(pH_finder, lower = 1, upper = 14, Alk = Alk, DIC = 5e-3, k1 = k1, k2 = k2)$root
  result = data.frame(Alk = Alk, pH = pH)
  if (i == 1) {
    sims = result
  } else {
    sims = rbind(sims, result)
  }
}
plot(sims$Alk * 1e3, sims$pH, type = "l",
     xlab = "Alkalinity (mmol/L)", ylab = "pH")
text(3, 10, label = "DIC = 5 mmol/L", font = 2)


# Calculate pH of various water types (unit: 1e-6 mol/L) ----
# rain ----
# Ca ~ [1, 10]
# Mg ~ [0.5, 5]
# Mn ~ [0, 0.1]
# Sr ~ [0, 0.1]
uniroot(pH_finder2, lower = 1, upper = 14,
             Ca = 5e-6, Mg = 2e-6, Mn = 1e-7, Sr = 1e-7, 
             DIC = 1.7e-5, k1 = k1, k2 = k2)$root
# river ----
# Ca ~ [50, 1500]
# Mg ~ [30, 800]
# Mn ~ [0.01, 1]
# Sr ~ [0.05, 5]
uniroot(pH_finder2, lower = 1, upper = 14,
        Ca = 1e-3, Mg = 2e-4, Mn = 1e-7, Sr = 5e-7, 
        DIC = 2.5e-3, k1 = k1, k2 = k2)$root

# Mn reduction ----
pH_finder3 = function(pH, Alk, DIC, k1, k2, amount){
  Alk_new = Alk + 2 * 2 * amount
  DIC_new = DIC + amount
  H = 10^(-pH)
  term1 = Alk_new / (1 + 2*k2/H)
  term2 = DIC_new / (H/k1 + 1 + k2/H)
  return(term1 - term2)
}
uniroot(pH_finder, c(1, 14), Alk = 4e-3, DIC = 5e-3, k1 = k1, k2 = k2)$root
uniroot(pH_finder3, c(1, 14), Alk = 4e-3, DIC = 5e-3, k1 = k1, k2 = k2, amount = 1e-3)$root
for (i in 1:100) {
  amount = 1e-4 + 1e-5 * i
  pH = uniroot(pH_finder3, c(1, 14), Alk = 4e-3, DIC = 5e-3, k1 = k1, k2 = k2, amount = amount)$root
  result = data.frame(amount = amount, pH = pH)
  if (i == 1) {
    sims = result
  } else {
    sims = rbind(sims, result)
  }
}
plot(sims$amount * 1e3, sims$pH, type = "l",
     xlab = expression("MnO"[2]*" (mmol/L)"), ylab = "pH")
text(.8, 8, label = "Alk = 4 mmol/L", font = 2)
text(.8, 7.5, label = "DIC = 5 mmol/L", font = 2)
