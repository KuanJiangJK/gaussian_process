plot(density(abs(rnorm(1e5, 0, 1))))
lines(density(rnorm(1e5, 0, 1)), col = "red")

plot(density(abs(rt(1e5, 1))))
