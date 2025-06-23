rm(list = ls())
pacman::p_load(deSolve)
# with(as.list(c(state, parms)), {})
# state = c(A, B) corresponds to list(c(dA, dB))
# t is time step, irrelevant to time

# initial values of parameters
state = c(V = 1) # unit: L
parms = c(F_in = 0.03, F_out = 0.01, k_evap = 0.01) # unit: L/s
time = seq(0, 1e3, 10) # unit: s

soil_water_volume = function(t, state, parms){
  with(as.list(c(state, parms)), {
    dV = F_in - F_out - k_evap * V
    list(c(dV))
  })
}

result = ode(y = state,
             times = time,
             parms = parms,
             func = soil_water_volume)
plot(time, result[,2], type = "l")
