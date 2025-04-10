test_that("signalDomain works", {
  sucrose <- new("NMRSignal1D"
                 ,id="Sucro"
                 ,chemicalShift=5.416
                 ,shape=list(name="pseudoVoigt"
                             ,params = list(mu=0.85,fwhm=0.0015,base=0)
                 )
                 ,peaks = list(new("NMRPeak1D",x=5.4113, y=1)
                               ,new("NMRPeak1D",x=5.4207, y=1)
                 )
  )
  expect_equal(signalDomain(sucrose),c(5.4113,5.4207) + 3*.0015*c(-1,1))
  quad <- new("NMRSignal1D"
              ,id="quadruplet"
              ,chemicalShift=3
              ,shape=list(name="pseudoVoigt"
                          ,params = list(mu=0.85,fwhm=0.002,base=0)
              )
              ,peaks = list(new("NMRPeak1D",x=2.85, y=.25)
                            ,new("NMRPeak1D",x=2.95, y=1)
                            ,new("NMRPeak1D",x=3.05, y=1)
                            ,new("NMRPeak1D",x=3.15, y=.25)
              )
  )
  expect_equal(signalDomain(quad,n=5),c(2.85,3.15) + 5*.002*c(-1,1))
})

test_that("signalToY works", {
  sucrose <- new("NMRSignal1D"
                 ,id="Sucrobation"
                 ,chemicalShift=5.416
                 ,shape=list(name="pseudoVoigt"
                             ,params = list(mu=0.85,fwhm=0.0015,base=0)
                 )
                 ,peaks = list(new("NMRPeak1D",x=5.4113, y=1)
                               ,new("NMRPeak1D",x=5.4207, y=1)
                 )
  )
  ppm <- seq(5.2,5.6,length.out=500)
  handmade <- fpseudoVoigt(xx=ppm,mean=5.4113,max=1,fwhm=0.0015,mu=0.85
                           ) + fpseudoVoigt(xx=ppm,mean=5.4207,max=1,fwhm=0.0015,mu=0.85)
  expect_equal(signalToY(sucrose,ppm),handmade)
})

test_that("scaleSignal works", {
  sucrose <- new("NMRSignal1D"
                 ,id="Sucrobation"
                 ,chemicalShift=5.416
                 ,shape=list(name="pseudoVoigt"
                             ,params = list(mu=0.85,fwhm=0.0015,base=0)
                 )
                 ,peaks = list(new("NMRPeak1D",x=5.4113, y=0.5)
                               ,new("NMRPeak1D",x=5.4207, y=0.5)
                 )
  )
  expect_equal(scaleSignal(sucrose,2)@peaks[[1]]@y, scaleSignal(sucrose,2)@peaks[[2]]@y)
  expect_equal(scaleSignal(sucrose,2)@peaks[[1]]@y, 2)
  expect_equal(scaleSignal(sucrose,by=2)@peaks[[1]]@y, scaleSignal(sucrose,by=2)@peaks[[2]]@y)
  expect_equal(scaleSignal(sucrose,by=2)@peaks[[1]]@y, 1)
})

test_that("shiftSignal works", {
  sucrose <- new("NMRSignal1D"
                 ,id="Sucrobation"
                 ,chemicalShift=5.416
                 ,shape=list(name="pseudoVoigt"
                             ,params = list(mu=0.85,fwhm=0.0015,base=0)
                 )
                 ,peaks = list(new("NMRPeak1D",x=5.4113, y=1)
                               ,new("NMRPeak1D",x=5.4207, y=1)
                 )
  )
  expect_equal(shiftSignal(sucrose,5.5)@peaks[[1]]@x, 5.5 -0.0047)
  expect_equal(shiftSignal(sucrose,5.5)@peaks[[2]]@x, 5.5 + 0.0047)
  expect_equal(shiftSignal(sucrose,5.5)@chemicalShift, 5.5)
  expect_equal(shiftSignal(sucrose,by=0.5)@peaks[[1]]@x, 5.4113 + 0.5)
  expect_equal(shiftSignal(sucrose,by=0.5)@peaks[[2]]@x, 5.4207 + 0.5)
  expect_equal(shiftSignal(sucrose,by=0.5)@chemicalShift, 5.416 + 0.5)
})
