set.seed(1)

# Per group sample size
n <- 50

x <- cbind(
  intercept = 1,
  exposure = rep(c(0,1), n)
  )

beta <- c(1.2, 2)

eta <- exp(x%*%beta)

lib.size <- sample(x = 200:400,
                   size = n*2,
                   replace = TRUE)

y <- rpois(n = n*2, lambda = lib.size * eta)


fit <- glm(y~x-1, family = poisson)
fit_1 <- glm(y~x-1, family = poisson, offset = lib.size)
fit_2 <- glm(y~x-1, family = poisson, offset = log(lib.size))
fit_3 <- glm(y~x-1 + offset(lib.size), family = poisson)
fit_4 <- glm(y~x-1 + offset(log(lib.size)), family = poisson)
