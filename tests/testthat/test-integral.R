test_that("integral works", {
  mu <- .85
  fwhm <- .0015
  h <- 2
  sucrose <- new("NMRSignal1D"
                 ,id="Sucrobation"
                 ,chemicalShift=5.416
                 ,shape=list(name="pseudoVoigt"
                             ,params = list(mu=mu,fwhm=fwhm,base=0)
                 )
                 ,peaks = list(new("NMRPeak1D",x=5.4113, y=h)
                               ,new("NMRPeak1D",x=5.4207, y=h)
                               )
                 )
  ppm <- seq(5.2,5.6,length.out=500)
  handmade <- fpseudoVoigt(ppm,5.4113,h
                           ,fwhm,mu) + fpseudoVoigt(ppm,5.4207,h,fwhm,mu)
  expect_equal(integral(sucrose,"fwhm")
               ,2 * mu * fwhm * h * pi/2 + 2 * (1-mu) * fwhm * h * 1.064467
  )
  expect_equal(round(integral(sucrose,"rect",offset=10),5),0.00876)
  expect_equal(integral(list(x=ppm,y=handmade),"rect")
               , sum(handmade) * (ppm[2]-ppm[1])
               )
})