library("bio3d")
library("rsvd")

X <- matrix(c(2, 8, 3, 5, 1, 3, 9, 2, 1, 3, 4, 2, 0, 0, 0),
              nrow = 5, byrow = T)

Y <- matrix(c(1, 5, 7, 2, 2, 2, 1, 4, 2, 8, 3, 1, 5, 4, 3),
            nrow = 5, byrow = T)

X <- matrix(c(3.74099994, -11.49300003,   6.4289999 , 
  7.13500023, -12.20800018,   8.04599953,
  8.02799988, -14.58100033,   5.24499989,
  7.1869998 , -11.95499992,   2.68099999,
  9.98799992,  -9.79500008,   3.98900008),
  nrow = 5, byrow = T, dimnames = list(NULL, c("x", "y", "z")) )

Y <- matrix(c(  1.83371258, -12.48775291,   8.70929527,
  5.34590292, -13.30946255,   7.76916122,
  8.41334724, -12.24099731,   5.81806755,
  8.36140156, -11.38534451,   2.10474873,
 12.12463951, -10.60843945,   1.98874044),
  nrow = 5, byrow = T, dimnames = list(NULL, c("x", "y", "z")) )
n <- nrow(X)
X
Y
# X <- c(2, 8, 3)
# Y <- c(1, 5, 7)
#Yr <- c(1.97385508,  7.89542034,  2.96078263)

X_Y <- X - Y

# compute the difference and get new 5 x 3
# square the difference
# sum row wise -> 5 x 1
# divide by 1 / n
# take the sqroot
# sum for every element
sum(sqrt((1 / n) * rowSums(X_Y ^ 2)))

# in the single row case. we have 1 x 3
# we should just sum the square diff
# 1 / n becomes one and can be removed
# since n = 1, there is only one element
# thus no final sum is required
sqrt(sum(X_Y[1,] ^ 2))
rmsd(a = X[1,], b = Y[1,]) # shuld be the same as above

acum <- 0
for (i in 1:nrow(X)) {
  new <- rmsd(a = X[i,], b = Y[i,])
  print(new)
  acum <- acum + new
}
acum / sqrt(5)

# X <- matrix(X)
# Y <- matrix(Y, ncol = length(Y))

R <- X %*% Y
result <- rsvd(R)

V <- result$u
Wt <- result$v
S <- result$d

n <- ncol(X)
sqrt((1 / n) * sum(X ^ 2 + t(Y ^ 2)) - 2 * sum(S))
