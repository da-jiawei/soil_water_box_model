# bell-shaped distribution
theta = seq(0, 1, 0.01)
C_atm = 4e2
peak = 0.8
a = 5
CO2_max = 2e4

b =  (a/peak) - a
peak = a / (a + b)
c = CO2_max / (peak^a * (1 - peak)^b)
CO2 = c * theta^a * (1 - theta)^b + C_atm
plot(theta, CO2, type = "l",
     xlab = expression(theta),
     ylab = expression("[CO"[2]*"]"[res]*" (ppmv)"))


# logistic function for sigmoidal curve
L = 1 # maximum value
k = 50 # steepness of the curve
xo = .9 # midpoint
y = L / (1 + exp(-k * (x - xo)))
plot(theta, y, type = "l",
     xlab = expression(theta),
     ylab = expression(italic(f)[ana]))
