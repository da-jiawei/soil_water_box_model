# root finding function for [H+] and d13C_HCO3
calcite_eq = function(H, ALK, DIC, k1, k2) {
  term1 = (2 * k1 * k2) / (H^2) + (k1 / H)
  term2 = 1 + (k1 * k2) / (H^2) + (k1 / H)
  return(ALK * term2 - DIC * term1)
}

find_d13HCO3 = function(X_HCO3, X_CO3, X_H2CO3, d13_DIC, d13_HCO3, alpha13_HCO3_H2CO3, alpha13_HCO3_CO2) {
  term1 = d13_HCO3 * X_HCO3
  term2 = ((d13_HCO3 + 1000) / alpha13_HCO3_H2CO3 - 1000) * X_H2CO3
  term3 = ((d13_HCO3 + 1000) / alpha13_HCO3_CO2 - 1000) * X_CO3
  return(term1 + term2 + term3 - d13_DIC)
}

# input parameters
ctrl = function(){
  vars = list(
    kp = 1e-3, # rate constant for carbonate precipitation (made up) - mol/L/s
    k_degas = 2e-7, # CO2 degassing constant - mol/s
    d13_co2_initial = -15, # the d13C of the initial soil CO2 in equilibrium with DIC
    d18p = -4, # rainfall d18O
    RH = 0.7, # relative humidity
    Tsoil = 15, # soil temperature
    w = 0.9, # the relative contribution of diffusion to evaporation
    # flux and concentration
    F_in = 0, # rainfall - L/s
    F_evap = 1e-2, # evaporation - L/s
    F_out = 0, # outflow of soil waters (e.g., to deeper soils, rootwater uptake)
    # rainfall chemistry
    Ca_p = 1e-4, # unit - mol/L
    Mg_p = 1e-5,
    Sr_p = 1e-6,
    DIC_p = 2.3e-4,
    d13_DIC_p = 0,
    # weathering input
    Ca_w = 0, # unit - mol/s
    Mg_w = 0,
    Sr_w = 0,
    DIC_w = 0,
    d13_DIC_w = -10,
    # parameterization
    dt = 1,  # unit - s
    time = 60 * 60,
    Vo = 10, # initial soil water volume - L
    d18so = -4, # initial soil water d18O
    Ca_so = 1e-3,
    Mg_so = 1e-4,
    Sr_so = 2e-6,
    DIC_so = 2.5e-3
  )
}

# box model for trace elements and stable isotopes in soil waters
SWTS_bm = function(vars) {
  ## Unpack variables
  list2env(vars, environment())
  ## constants ----
  # the first and second disassociation constants for carbonic acid
  k1 = 10^-6.3
  k2 = 10^-10.3
  kh = 10^-1.5
  
  ksp = 3.3e-9 # solubility product of calcite
  kd_Sr = 0.057 # partition coefficient for Sr
  kd_Mg = 0.031 # partition coefficient for Mg
  
  # isotope standards
  R18smow = 0.0020052
  R18vpdb = 0.0020672
  
  CO2_atm = 400
  
  # molarity of water
  H2O_mol = 55.51 # mol/L
  
  ## stable isotope ----
  Tsoil.K = Tsoil + 273.15
  # oxygen isotope
  alpha18_c_w = exp((1.61e4 / Tsoil.K - 24.6) / 1000) # Tremaine (2011)
  alpha18_l_v = exp(-2.0667 * 10^-3 - 0.4156 / Tsoil.K + (1.137 * 10 ^ 3) / (Tsoil.K ^ 2)) # Majoube (1971)
  alpha18_diff_p = 1.028489
  alpha18_diff = alpha18_diff_p * w + (1-w)
  R18p = (d18p / 1000 + 1) * R18smow
  # assuming atmospheric vapor in equilibrium with rainfall under soil temperature
  R18v = R18p / alpha18_l_v
  
  # carbon isotope 
  alpha13_HCO3_CO3 = exp((-0.87 * 1e3 / Tsoil.K + 2.52) * 1e-3) # Mook et al. (1974)
  alpha13_H2CO3_CO2 = exp((6e3 / Tsoil.K^2 - 0.9) * 1e-3) # Deines et al. (1974)
  alpha13_HCO3_CO2 = exp((1.1e6 / Tsoil.K^2 - 4.5) * 1e-3) # Deines et al. (1974)
  alpha13_HCO3_H2CO3 = alpha13_HCO3_CO2 / alpha13_H2CO3_CO2
  alpha13_cal_CO2 = (11.98 - 0.12 * Tsoil) * 1e-3 + 1 # Romanek et al. (1992)
  
  
  # parameterization ----
  num = time / dt
  V = rep(Vo, num) # unit - L
  d18_s = rep(d18so, num)
  d18_c = rep(0, num)
  d13_DIC = rep(0, num)
  d13_co2 = rep(0, num) # degassing CO2 in equilibrium with H2CO3
  d13_H2CO3 = rep(0, num)
  d13_HCO3 = rep(0, num)
  d13_CO3 = rep(0, num)
  d13_c = rep(0, num)
  Ca_s = rep(Ca_so, num) # unit - mol/L
  Mg_s = rep(Mg_so, num)
  Sr_s = rep(Sr_so, num)
  DIC_s = rep(DIC_so, num)
  CO3_s = rep(0, num)
  HCO3_s = rep(0, num)
  H2CO3_s = rep(0, num)
  CO2_s = rep(0, num)
  H_s = rep(0, num)
  pH = rep(0, num)
  MgCa_s = rep(0, num)
  SrCa_s = rep(0, num)
  Jp = rep(0, num)
  SrCa_c = rep(0, num)
  MgCa_c = rep(0, num)
  
  ## loop ----
  tryCatch(
    {
      for (i in 1:(num-1)) {
        # calculating water isotope
        R18s = (d18_s[i]/1000 + 1) * R18smow
        R18e = (R18s - RH * R18v * alpha18_l_v) / ((1 - RH) * alpha18_l_v * alpha18_diff) # Craig-Gordon model
        d18e = (R18e / R18smow - 1) * 1000
        
        # solve [H+] and DIC species using alkalinity and DIC
        Alk = 2 * (Ca_s[i] + Mg_s[i] + Sr_s[i]) 
        H_s[i] = uniroot(calcite_eq, c(1e-14, 1e-1), ALK = Alk, DIC = DIC_s[i], k1 = k1, k2 = k2, tol = 1e-14)$root
        pH[i] = -log10(H_s[i])
        CO3_s[i] = DIC_s[i] / ((H_s[i]^2) / (k1 * k2) + H_s[i] / k2 + 1) # carbonate ion
        HCO3_s[i] = H_s[i] * CO3_s[i] / k2
        H2CO3_s[i] = H_s[i] * HCO3_s[i] / k1
        CO2_s[i] = 1e6 * H2CO3_s[i] / kh
        X_CO3 = CO3_s[i] / DIC_s[i]
        X_HCO3 = HCO3_s[i] / DIC_s[i]
        X_H2CO3 = H2CO3_s[i] / DIC_s[i]
        
        if(i == 1) {
          # solve DIC-d13C using mass balance model
          d13_HCO3[i] = alpha13_HCO3_CO2 * (d13_co2_initial + 1000) - 1000
          d13_CO3[i] =  (d13_HCO3[i] + 1000) / alpha13_HCO3_CO3 - 1000
          d13_H2CO3[i] = alpha13_H2CO3_CO2 * (d13_co2_initial + 1000) - 1000
          d13_DIC[i] = d13_CO3[i] * X_CO3 + d13_HCO3[i] * X_HCO3 + d13_H2CO3[i] * X_H2CO3
        } else {
          d13_HCO3[i] = uniroot(find_d13HCO3, X_HCO3 = X_HCO3, X_CO3 = X_CO3, X_H2CO3 = X_H2CO3, d13_DIC = d13_DIC[i], 
                                alpha13_HCO3_H2CO3 = alpha13_HCO3_H2CO3, alpha13_HCO3_CO2 = alpha13_HCO3_CO2, c(-20, 10))$root
          d13_CO3[i] = (d13_HCO3[i] + 1000) / alpha13_HCO3_CO3 - 1000
          d13_H2CO3[i] = (d13_HCO3[i] + 1000) / alpha13_HCO3_H2CO3 - 1000
        }
        d13_co2[i] = (d13_H2CO3[i] + 1000) / alpha13_H2CO3_CO2 - 1000
        
        omega = Ca_s[i] * CO3_s[i] / ksp # the saturation state of calcite
        if(omega > 1) {
          Jp[i] = kp * (omega - 1) # precipitation flux of calcite - unit: mol/L/s
        } else {
          Jp[i] = 0
        }
        # box model
        dV = dt * (F_in - F_evap - F_out) # unit: L/s
        degas = k_degas * (CO2_s[i] - CO2_atm)
        dDIC = (dt / V[i]) * (F_in * (DIC_p - DIC_s[i]) + F_evap * DIC_s[i] + DIC_w - Jp[i] - degas)
        dCa = (dt / V[i]) * (F_in * (Ca_p - Ca_s[i]) + F_evap * Ca_s[i] + Ca_w - Jp[i])
        dMg = (dt / V[i]) * (F_in * (Mg_p - Mg_s[i]) + F_evap * Mg_s[i] + Mg_w - (Jp[i] * kd_Mg * Mg_s[i] / Ca_s[i]))
        dSr = (dt / V[i]) * (F_in * (Sr_p - Sr_s[i]) + F_evap * Sr_s[i] + Sr_w - (Jp[i] * kd_Sr * Sr_s[i] / Ca_s[i]))
        dd18_s = (dt / V[i]) * ((d18p - d18_s[i]) * F_in - (d18e - d18_s[i]) * F_evap - d18_c[i] * (3 * Jp[i] / H2O_mol))
        dd13_DIC = (dt / (V[i] * DIC_s[i])) * ((d13_DIC_p - d13_DIC[i]) * F_in - (d13_co2[i] * degas) - (d13_c[i] * Jp[i]) + (d13_DIC[i] * DIC_s[i] * F_evap) - (d13_DIC[i] * V[i] * dDIC / dt))
        
        if(Jp[i] == 0) {
          SrCa_c[i] = NA
          MgCa_c[i] = NA
          d13_c[i] = NA
          d18_c[i] = NA
        } else {
          SrCa_s[i] = Sr_s[i] / Ca_s[i]
          MgCa_s[i] = Mg_s[i] / Ca_s[i]
          SrCa_c[i] = kd_Mg * SrCa_s[i]
          MgCa_c[i] = kd_Sr * MgCa_s[i]
          d13_c[i] = (d13_co2[i] + 1000) * alpha13_cal_CO2 - 1000
          R18c = R18s * alpha18_c_w
          d18_c[i] = (R18c / R18vpdb - 1) * 1000
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
        d13_DIC[i+1] = d13_DIC[i] + dd13_DIC
      }
      }, error = function(e) {
      message("Error encountered: ", e$message)
      message("Returning results up to the last successful iteration.")
    })
  results = data.frame(time = seq(dt, time, dt), V = V, fraction = V/V[1], DIC = DIC_s, CO2_s = CO2_s,
                   pH = pH, MgCa_s = MgCa_s, SrCa_s = SrCa_s, MgCa_c = MgCa_c, SrCa_c = SrCa_c,
                   Jp = Jp, d18s = d18_s, d18c = d18_c, d13_DIC = d13_DIC, d13_co2 = d13_co2, d13c = d13_c)[1:(i-1), ]
}

