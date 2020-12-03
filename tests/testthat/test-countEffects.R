test_that("two-group intercept only Poisson", {
  skip()
  fit <- countEffects(y = "dv",
                      x = "treat",
                      k = "k",
                      data = d,
                      method = "poisson")
  fit
  fit@results@est
})


test_that("two-group single manifest covariate Poisson", {
  skip()
  fit <- countEffects(y = "dv",
                      x = "treat",
                      k = "k",
                      z = "z12",
                      data = d,
                      method = "poisson")
  fit
  fit@results@est
})


test_that("two-group three manifest covariate Poisson", {
  skip()
  fit <- countEffects(y = "dv",
                      x = "treat",
                      k = "k",
                      z = c("z12", "z11", "z21"),
                      data = d,
                      method = "poisson")
  fit
  fit@results@est
})



test_that("two-group one latent covariate Poisson", {
  skip()
  fit <- countEffects(y = "dv",
                      x = "treat",
                      k = "k",
                      z = c("eta1"),
                      measurement = list(eta1 = c("z41", "z42", "z43")),
                      data = d,
                      method = "negbin")
  fit
  fit@results@est
})



test_that("two-group two latent covariate Poisson", {
  skip()
  fit <- countEffects(y = "dv",
                      x = "treat",
                      k = "k",
                      z = c("eta1", "eta2"),
                      measurement = list(eta1 = c("z41", "z42", "z43"),
                                         eta2 = c("z21", "z22")),
                      data = d,
                      method = "poisson")
  fit
  fit@results@est
})

