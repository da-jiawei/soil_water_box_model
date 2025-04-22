# bell-shaped distribution
x = seq(0, 1, 0.01)
x.max = 0.8 # volumetric soil water content with maximum respiration
n = 1 / x.max - 1
a = 1 / (x.max * (1 - x.max) ^ n)
y = a * x * (1 - x) ^ n
plot(x, y, type = "l",
     xlab = expression(italic(f)[sw]),
     ylab = expression("[CO"[2]*"]"[res]))


# logistic function for sigmoidal curve
x = seq(0, 1, 0.01)
L = 1 # maximum value
k = 50 # steepness of the curve
xo = 0.9 # midpoint
y = L / (1 + exp(-k * (x - xo)))
plot(x, y, type = "l",
     xlab = expression(italic(f)[sw]),
     ylab = expression("[CO"[2]*"]"[hetero]*" / [CO"[2]*"]"[total]))
