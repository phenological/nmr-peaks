test_that("signal shapes lorentzian", {
  xx <- 0:100
  meanx <- 50
  maxy <- 123
  fwhm <- 6
  y <- florentzian(xx, meanx, maxy, fwhm)
  expect_equal(max(y), maxy, tolerance=1e-3)
  expect_equal(xx[[which(y == max(y))]], meanx)
  expect_equal(y[[meanx - fwhm/2 + 1]], maxy / 2, tolerance=1e-2)
  expect_equal(y[[meanx + fwhm/2 + 1]], maxy / 2, tolerance=1e-2)
  
  xx <- 0:200
  meanx <- 80
  maxy <- 1000
  fwhm <- 16
  y <- florentzian(xx, meanx, maxy, fwhm)
  expect_equal(max(y), maxy, tolerance=1e-3)
  expect_equal(xx[[which(y == max(y))]], meanx)
  expect_equal(y[[meanx - fwhm/2 + 1]], maxy / 2, tolerance=1e-2)
  expect_equal(y[[meanx + fwhm/2 + 1]], maxy / 2, tolerance=1e-2)
})

test_that("signal shapes gaussian", {
  xx <- 100:200
  meanx <- 150
  maxy <- 123
  fwhm <- 4
  y <- fgaussian(xx, meanx, maxy, fwhm)
  expect_equal(max(y), maxy, tolerance=1e-3)
  expect_equal(xx[[which(y == max(y))]], meanx)
  print(meanx - fwhm/2 + 1)
  expect_equal(y[[meanx - fwhm/2 + 1 - 100]], maxy / 2, tolerance=1e-2)
  expect_equal(y[[meanx + fwhm/2 + 1 - 100]], maxy / 2, tolerance=1e-2)
  
  xx <- 100:300
  meanx <- 180
  maxy <- 1000
  fwhm <- 12
  y <- fgaussian(xx, meanx, maxy, fwhm)
  expect_equal(max(y), maxy, tolerance=1e-3)
  expect_equal(xx[[which(y == max(y))]], meanx)
  expect_equal(y[[meanx - fwhm/2 + 1 - 100]], maxy / 2, tolerance=1e-2)
  expect_equal(y[[meanx + fwhm/2 + 1 - 100]], maxy / 2, tolerance=1e-2)
})


test_that("signal shapes pseudoVoigt", {
  xx <- 0:100
  meanx <- 50
  maxy <- 123
  fwhm <- 10
  y <- fpseudoVoigt(xx, 0, meanx, maxy, fwhm)
  expect_equal(max(y), maxy, tolerance=1e-3)
  expect_equal(xx[[which(y == max(y))]], meanx)
  expect_equal(y[[meanx - fwhm/2 + 1]], maxy / 2, tolerance=1e-2)
  expect_equal(y[[meanx + fwhm/2 + 1]], maxy / 2, tolerance=1e-2)
  
  xx <- 0:200
  meanx <- 80
  maxy <- 1000
  fwhm <- 20
  y <- fpseudoVoigt(xx, 1, meanx, maxy, fwhm)
  expect_equal(max(y), maxy, tolerance=1e-3)
  expect_equal(xx[[which(y == max(y))]], meanx)
  expect_equal(y[[meanx - fwhm/2 + 1]], maxy / 2, tolerance=1e-2)
  expect_equal(y[[meanx + fwhm/2 + 1]], maxy / 2, tolerance=1e-2)
})
