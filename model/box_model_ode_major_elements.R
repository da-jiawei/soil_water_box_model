pacman::p_load(deSolve)
# functions ----
# find proton
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
# find d13C_HCO3
find_d13HCO3 = function(X_HCO3, X_CO3, X_H2CO3, d13_DIC, d13_HCO3, alpha13_HCO3_H2CO3, alpha13_HCO3_CO2) {
  term1 = d13_HCO3 * X_HCO3
  term2 = ((d13_HCO3 + 1e3) / alpha13_HCO3_H2CO3 - 1e3) * X_H2CO3
  term3 = ((d13_HCO3 + 1e3) / alpha13_HCO3_CO2 - 1e3) * X_CO3
  return(term1 + term2 + term3 - d13_DIC)
}

# parameters and initial state ---- 
parms = c(
  # parameters associated with Mn reduction
  peak = 0.8, # volumetric water content with maximum respiration
  r_rd = 5, # respiration reduction rate
  steep = 10, # parameters determining the relative amount of heterotrophic CO2 respiration
  xmid = 0.93, # parameters determining the relative amount of heterotrophic CO2 respiration
  Mn_rd = 5e-11, # Mn reduction rate - mol/s
  
  # rate constant 
  kp = 3e-3, # rate constant for carbonate precipitation - mol/s
  k_degas = 2e-7, # CO2 degassing constant - mol/s
  k_diss = 1e-2, # CO2 dissolution constant
  lamda = 1e-3, # Mn reduction reaction rate constant
  
  # environmental parameters
  CO2_atm = 400, # atmospheric CO2 - ppmv
  d13_r = -20, # the d13C of the respired CO2 in equilibrium with DIC
  CO2_r_max = 5e-8, # the maximum amount of respired CO2 being added into the DIC - mol/s
  d18p = -4, # rainfall d18O
  RH = 0.7, # relative humidity
  Tsoil = 20, # soil temperature
  w = 0.8, # the relative contribution of diffusion to evaporation
  
  # mass action constants
  k1 = 10^-6.3,
  k2 = 10^-10.3,
  kh = 10^-1.5,
  ksp = 3.3e-9, # solubility product of calcite
  
  # isotope standards
  R18smow = 0.0020052,
  R18vpdb = 0.0020672,
  
  # flux and concentration
  F_in = 3e-5, # rainfall - L/s (based on MAP)
  F_evap = 5e-4, # evaporation - L/s
  F_out = 0, # outflow of soil waters (e.g., to deeper soils, rootwater uptake)
  # rainfall chemistry
  Ca_p = 1e-4, # unit - mol/L
  Mg_p = 1e-5,
  Sr_p = 1e-6,
  Mn_p = 1e-6,
  DIC_p = 2.2e-4,
  d13_DIC_p = -3,
  # weathering input
  Ca_w = 0, # unit - mol/s
  Mg_w = 0,
  Sr_w = 0,
  Mn_w = 0,
  DIC_w = 0,
  d13_DIC_w = -25,
  
  Vt = 11, # maximum soil water volume
  
  # distribution coefficients
  kd_Mg = 0.031,
  kd_Sr = 0.057,
  kd_Mn = 2
)

state = c(V = 10, # initial soil water volume - L
          d18_sw = -4,
          d13_DIC = -15,
          Ca_sw = 4e-3,
          Mg_sw = 4e-4,
          Na_sw = 6e-4,
          K_sw = 2e-4,
          Sr_sw = 4e-6,
          Mn_sw = 1e-5,
          DIC_sw = 1.48e-2)

# ode function ----
soil_water_chemistry = function(t, state, parms){
  with(as.list(c(state, parms)), {
    Tsoil.K = Tsoil + 273.15
    # oxygen isotope
    alpha18_c_w = exp((1.61e4 / Tsoil.K - 24.6) / 1e3) # Tremaine (2011)
    alpha18_l_v = exp(-2.0667 * 10^-3 - 0.4156 / Tsoil.K + (1.137 * 10 ^ 3) / (Tsoil.K ^ 2)) # Majoube (1971)
    alpha18_diff_p = 1.028489
    alpha18_diff = alpha18_diff_p * w + (1-w)
    R18p = (d18p / 1e3 + 1) * R18smow
    # assuming atmospheric vapor in equilibrium with rainfall under soil temperature
    R18v = R18p / alpha18_l_v
    
    # carbon isotope fractionation factor
    alpha13_HCO3_CO3 = exp((-0.87 * 1e3 / Tsoil.K + 2.52) * 1e-3) # Mook et al. (1974)
    alpha13_H2CO3_CO2 = exp((6e3 / Tsoil.K^2 - 0.9) * 1e-3) # Deines et al. (1974)
    alpha13_HCO3_CO2 = exp((1.1e6 / Tsoil.K^2 - 4.5) * 1e-3) # Deines et al. (1974)
    alpha13_HCO3_H2CO3 = alpha13_HCO3_CO2 / alpha13_H2CO3_CO2
    alpha13_cal_CO2 = (11.98 - 0.12 * Tsoil) * 1e-3 + 1 # Romanek et al. (1992)
    
    # time series parameters ----
    R18s = (d18_sw/1e3 + 1) * R18smow
    R18e = (R18s - RH * R18v * alpha18_l_v) / ((1 - RH) * alpha18_l_v * alpha18_diff) # Craig-Gordon model
    d18e = (R18e / R18smow - 1) * 1e3
    
    Alk = 2 * (Ca_sw + Mg_sw + Sr_sw + Mn_sw + K_sw + Na_sw) 
    pH = uniroot(pH_finder, c(1, 14), Alk = Alk, DIC = DIC_sw, k1 = k1, k2 = k2)$root
    H_sw = 10^(-pH)
    # DIC = H2CO3 + HCO3 + CO3 = H*HCO3/k1 + HCO3 + HCO3*k2/H
    # HCO3 = DIC / (H/k1 + 1 + k2/H)
    HCO3_sw = DIC_sw / (H_sw/k1 + 1 + k2/H_sw)
    CO3_sw = HCO3_sw * k2 / H_sw
    H2CO3_sw = H_sw * HCO3_sw / k1
    CO2_sw = 1e6 * H2CO3_sw / kh
    X_CO3 = CO3_sw / DIC_sw
    X_HCO3 = HCO3_sw / DIC_sw
    X_H2CO3 = H2CO3_sw / DIC_sw
    d13_HCO3 = uniroot(find_d13HCO3, X_HCO3 = X_HCO3, X_CO3 = X_CO3, X_H2CO3 = X_H2CO3, d13_DIC = d13_DIC, 
                          alpha13_HCO3_H2CO3 = alpha13_HCO3_H2CO3, alpha13_HCO3_CO2 = alpha13_HCO3_CO2, c(-20, 10))$root
    d13_CO3 = (d13_HCO3 + 1e3) / alpha13_HCO3_CO3 - 1e3
    d13_H2CO3 = (d13_HCO3 + 1e3) / alpha13_HCO3_H2CO3 - 1e3
    d13_co2 = (d13_H2CO3 + 1e3) / alpha13_H2CO3_CO2 - 1e3
    
    # calcite precipitation 
    omega = Ca_sw * CO3_sw / ksp 
    if(omega > 1) {
      Jp = kp * (omega - 1) # precipitation flux of calcite - mol/s
      # cal_rate = ifelse(omega <= 2.1, 0.28 * omega, (7.5 * omega - 15)) # nmol/mg seed crystal/min (Lorens, 1981)
      # Jp = cal_rate * 1e-9 / 60
      # distribution coefficients - Lorens 1981
      # kd_Sr = exp(0.249 * log(cal_rate) - 1.57)
      # kd_Mn = exp(-0.266 * log(cal_rate) + 1.35)
    } else {
      Jp = 0
    }
    
    SrCa_sw = Sr_sw / Ca_sw
    MgCa_sw = Mg_sw / Ca_sw
    MnCa_sw = Mn_sw / Ca_sw
    
    if(Jp == 0) {
      SrCa_c = 0
      MgCa_c = 0
      MnCa_c = 0
      d13_c = 0
      d18_c = 0
    } else {
      SrCa_c = kd_Sr * SrCa_sw
      MgCa_c = kd_Mg * MgCa_sw
      MnCa_c = kd_Mn * MnCa_sw
      d13_c = (d13_co2 + 1e3) * alpha13_cal_CO2 - 1e3
      R18c = R18s * alpha18_c_w
      d18_c = (R18c / R18vpdb - 1) * 1e3
    }
    
    # Mn reduction
    # C input from soil respiration as a function of remaining soil water
    theta = V / Vt
    b = (r_rd/peak) - r_rd
    c =  CO2_r_max / (peak^r_rd * (1 - peak)^b)
    # total respired CO2 being added into the DIC - mol/s
    CO2_r = c * theta^r_rd * (1 - theta)^b 
    # relative amount of anaerobic CO2 
    CO2_ana = CO2_r * (1 / (1 + exp(-steep * (theta - xmid))))
    # relative amount of Mn reduction
    # 4H + CH2O + 2MnO2 = 2Mn + CO2 + 3H2O
    # Mn_rd = 2 * lamda * CO2_ana
    HCO3_rd = 0.5 * Mn_rd
    
    degas = k_degas * (CO2_sw - CO2_atm)
    CO2_diss = k_diss * CO2_r
    
    # ODE
    dV = F_in - F_evap - F_out # unit: L/s
    dDIC = 1/V * (F_in * (DIC_p - DIC_sw) + F_evap * DIC_sw + DIC_w - Jp - degas + CO2_diss + HCO3_rd)
    dCa = 1/V * (F_in * (Ca_p - Ca_sw) + F_evap * Ca_sw + Ca_w - Jp)
    dK = 1/V * (F_in * (0 - K_sw) + F_evap * K_sw + 0)
    dNa = 1/V * (F_in * (0 - Na_sw) + F_evap * Na_sw + 0)
    dMg = 1/V * (F_in * (Mg_p - Mg_sw) + F_evap * Mg_sw + Mg_w - (Jp * MgCa_c))
    dSr = 1/V * (F_in * (Sr_p - Sr_sw) + F_evap * Sr_sw + Sr_w - (Jp * SrCa_c))
    dMn = 1/V * (F_in * (Mn_p - Mn_sw) + F_evap * Mn_sw + Mn_w - (Jp * MnCa_c) + Mn_rd) 
    dd18_sw = 1/V * ((d18p - d18_sw) * F_in - (d18e - d18_sw) * F_evap)
    dd13_DIC = (1 / (V * DIC_sw)) * ((d13_DIC_p - d13_DIC) * F_in - (d13_co2 - d13_DIC) * degas + 
                                       (d13_r - d13_DIC) * (CO2_diss + HCO3_rd) - 
                                       (d13_c - d13_DIC) * Jp + (d13_DIC_w - d13_DIC) * DIC_w)
    
    list(c(dV, dd18_sw, dd13_DIC, dCa, dK, dNa, dMg, dSr, dMn, dDIC),
         CO2_r = CO2_r, Mn_rd = Mn_rd, Jp = Jp, pH = pH, Alk = Alk, omega = omega,
         degas = degas, HCO3_rd = HCO3_rd, CO2_diss = CO2_diss, 
         SrCa_c = SrCa_c, MgCa_c = MgCa_c, MnCa_c = MnCa_c,
         d13_c = d13_c, d18_c = d18_c)
  })
}

# parms["Vt"] = 11
# time = seq(0, 1e4, 1)
# results = ode(y = state,
#               times = time,
#               parms = parms,
#               func = soil_water_chemistry)
# results = as.data.frame(results)
# plot(results$time, results$CO2_r, type = "l")
# plot(results$time, results$d13_c, type = "l")
# plot(results$time, results$d18_c, type = "l")
# plot(results$time, results$Mn_rd, type = "l")
# plot(results$time, results$MnCa_c, type = "l")
