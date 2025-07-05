# bell-shaped distribution
theta = seq(0, 0.9, 0.01)
C_atm = 4e2
peak = 0.7
a = 2
CO2_max = 2e4

b =  (a/peak) - a
peak = a / (a + b)
c = CO2_max / (peak^a * (1 - peak)^b)
CO2 = c * theta^a * (1 - theta)^b + C_atm
plot(theta, CO2, type = "l",
     xlab = expression(theta),
     ylab = expression("[CO"[2]*"]"[res]*" (ppmv)"))


# logistic function for sigmoidal curve
lamda = 1e-11 # maximum value
steep = 5 # steepness of the curve
x_mid = .8 # midpoint
Mn_rd = CO2 * lamda / (1 + exp(-steep * (theta - x_mid)))
plot(theta, Mn_rd, type = "l",
     xlab = expression(theta),
     ylab = expression(italic(f)[ana]))
