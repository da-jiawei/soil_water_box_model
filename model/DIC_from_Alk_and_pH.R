rm(list = ls())
k1 = 10^-6.3
k2 = 10^-10.3


DIC_finder = function(DIC, Alk, pH, k1, k2){
  # Alk = HCO3 + 2CO3 = HCO3 + 2*HCO3*k2/H
  # DIC = H2CO3 + HCO3 + CO3 = H*HCO3/k1 + HCO3 + HCO3*k2/H
  # HCO3 = DIC / (H/k1 + 1 + k2/H)
  H = 10^(-pH)
  HCO3 = Alk / (1 + 2*k2/H)
  term2 = DIC / (H/k1 + 1 + k2/H)
  return(HCO3 - term2)
}


Alk = 1.04e-2
pH = 6.67
uniroot(DIC_finder, c(1e-3, 1e-1), Alk = Alk, pH = pH, k1 = k1, k2 = k2)$root
