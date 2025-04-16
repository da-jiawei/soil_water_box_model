# bell-shaped distribution
x = seq(0, 1, 0.01)
x.max = 0.8 # volumetric soil water content with maximum respiration
n = 1 / x.max - 1
a = 1 / (x.max * (1 - x.max) ^ n)
y = a * x * (1 - x) ^ n
plot(x, y, type = "l")


# logistic function for sigmoidal curve
x = seq(0, 1, 0.01)
L = 1 # maximum value
k = 20 # steepness of the curve
xo = 0.8 # midpoint
y = L / (1 + exp(-k * (x - xo)))
plot(x, y, type = "l")
