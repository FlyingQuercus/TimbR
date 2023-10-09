
# Function for estimating free-free modal frequencies of a uniform prismatic specimen of wood with rectangular cross-section, properly taking into account G
# By Dan Ridley-Ellis (d.ridleyellis@napier.ac.uk)
# Last updated 13/09/2023
#
# Args:
# dims = c(#1,#2,#3) for rectangular cross-section 3 dimensions: length, width and thickness (mm), specified in any order
# mass = mass (g) or
# density = density (kg/m^3)
# E = Modulus of elasticity (kN/mm^2) {along the length, x direction}
# G = Shear modulus (kN/mm^2) {xy and xz planes} (use G = "ISO" if you want the value based on E, Poisson's ratio and assumption of isotropic material)
# Poisson.ratio = Poisson's ratio (0.4 is default, used for the longitudinal mode correction factor, it is Poisson.ratio.xy = Poisson.ratio.xz) (also used in the form factor if using https://doi.org/10.1007/s11012-022-01618-1)
# mode.count = #n or c(#n1,#n2,#n3) the number of frequencies to calculate for each mode type: on number for all or specified in the order flexural, torsional and longitudinal (12 is default)
# torsion.fast = TRUE or FALSE {if TRUE the much faster approximate method is used for the torsional calculation}
# form.factor.minor = 1 / Timoshenko shear coefficient minor axis flexural (Defaults to 1.2. Use NULL if you want to use https://doi.org/10.1007/s11012-022-01618-1)
# form.factor.major = 1 / Timoshenko shear coefficient minor axis flexural (Defaults to 1.2. Use NULL if you want to use https://doi.org/10.1007/s11012-022-01618-1)
# eta.squared = Material parameter for torsion (0.98 is default) 
#               see https://doi.org/10.1007/s11029-008-9007-z
#               "transtropic and orthotropic fibrous composites whose longitudinal axis coincides with the reinforcement direction, this parameter - as shown by calculation - changes in the narrow range 0.969 to 0.999"
#
# Returns:
# A list of lists of frequencies (Hz)
#    Flexural minor
#    Flexural major
#    Torsional
#    Longitudinal
#    These as a single named list in frequency order
# The values on which the calculation was based

#
# Notes:
# The flexural frequencies are estimates using Loic Brancheriau's method in https://doi.org/10.1007/s00226-014-0672-x
# Set the form.factor.minor and form.factor.major to NULL to use shear coefficients calculated according to https://doi.org/10.1007/s11012-022-01618-1 based on Poisson's ratio and aspect ratio
# The torsional frequencies use Ekel'chik's method in https://doi.org/10.1007/s11029-008-9007-z (the accurate method is the default)
#     Applied assuming that G in the radial-longitudinal plane = G in the tangential-longitudinal plane (Gxy = Gxz, where x is longitudinal)
#     If asking for many torsional frequencies and/or the specimen is not slender the search needed to solve the accurate method may take a very long time, or just stop immediately with an error
#     ...in which case you can use torsion.fast = TRUE to speed things up, but be aware a slow accurate calculation is a sign that the fast method will not be very accurate
# The longitudinal frequencies use Loic Brancheriau's method for correction in ISBN-13:978-1616682170 pp 205-224


IEfreqs = function (dims = NULL, mass = NULL, density = NULL, E = NULL, G = NULL, Poisson.ratio = 0.4, mode.count = 12, torsion.fast = FALSE, form.factor.minor = 1.2, form.factor.major = 1.2, eta.squared = 0.98, maxfreq = 22200) {
  
  # Make sure the required values are given 

  if (is.null(E)) 
    stop("You must specify E")
  if (is.null(dims)) 
    stop("You must specify dims")
  if (length(dims) != 3)
    stop("You need dims of length 3")
  if (length(mode.count) == 1)
  mode.count = c(mode.count, mode.count, mode.count)  
  if (length(mode.count) != 3)
    stop("You need mode count of length 3 (flexural, torsional, longitudinal")
  if (min(mode.count) < 2)
    stop("You need to ask for at least 2 frequencies per mode")

  if (!is.null(mass) && !is.null(density))   
    stop("Specify either 'mass' or 'density', but not both")
  if (is.null(mass) && is.null(density)) 
    stop("You must specify either 'mass' or 'density'")

  # Flexural Bernoulli solutions (solutions of cos(K) = 1 / cosh(K) )
  # For n of 5 and above the approximation is used: K = pi * (1 + 2 * n) / 2
  K = seq(from = 1, to = mode.count[1])
  K = pi * (1 + 2 * K) / 2
  # These are more accurate solutions for the first four
  K[1] = 4.730040745
  K[2] = 7.853204624
  K[3] = 10.99560784
  K[4] = 14.13716549
  
  # Preliminary calculations
  
  L = max(dims) # the length
  h = min(dims) # the thickness
  b = median(dims) # the width

  L = L / 1000 # metres
  h = h / 1000 # metres
  b = b / 1000 # metres
  
  # https://doi.org/10.1007/s11012-022-01618-1
  if(is.null(form.factor.major)){
  s = b/h
  form.factor.major = 1/( (20*(1+Poisson.ratio)^2*s^2*(1+2*s^2)^2) / (Poisson.ratio*(5+39*Poisson.ratio)+(24+Poisson.ratio*(57+79*Poisson.ratio))*s^2+4*(4+3*Poisson.ratio)*(6+7*Poisson.ratio)*s^4+12*(1+Poisson.ratio)*(8+7*Poisson.ratio)*s^6) )
  cat("form.factor.major is calculated to be", form.factor.major ,"\n")
  }
  if(is.null(form.factor.minor)){
  s = h/b
  form.factor.minor = 1/( (20*(1+Poisson.ratio)^2*s^2*(1+2*s^2)^2) / (Poisson.ratio*(5+39*Poisson.ratio)+(24+Poisson.ratio*(57+79*Poisson.ratio))*s^2+4*(4+3*Poisson.ratio)*(6+7*Poisson.ratio)*s^4+12*(1+Poisson.ratio)*(8+7*Poisson.ratio)*s^6) )
  cat("form.factor.minor is calculated to be", form.factor.minor ,"\n")
  }

  if(L/b < 10 ) cat("Warning: the specimen aspect ratio is poor for these equations \n")
  if(L/b < 10 && L/h >= 10 ) cat("Warning: the specimen aspect ratio is poor for major axis flexural \n")
    
  if (is.null(mass)) {
    mass = density * 1000 * (L) * (b) * (h)     
  }
    
  if (is.null(density)) {
    density = (mass / 1000) / ((L) * (b) * (h))
  }
    
  if (is.null(G) && !is.null(E)) {
    EG.ratio = 16
    cat("You did not specify G, so using E/G =", EG.ratio, "(an assumption for wood) \n")
    G = E / EG.ratio
    cat("This means G is", G, "kN/mm^2 \n")
  }

  if (G == "ISO") {
    cat("You asked for G to be calculated as if an isotropic material \n")
    G = E / (2 * (1 + Poisson.ratio))
    cat("This means G is", G, "kN/mm^2 \n")
  }
  
  # convert E and G to compatible units (N / m^2)
  E = E * 10^9
  G = G * 10^9      
  
  Area = b * h #m^2
  #cat("Cross section area", Area * 1000^2, "\n")
  I.minor = (1/12) * b * h^3 #m^4
  #cat("I minor axis", I.minor * 1000^4, "\n")
  I.major = (1/12) * h * b^3 #m^4
  #cat("I major axis", I.major * 1000^4, "\n")
  i.minor = (I.minor / Area)^0.5
  #cat("rad gyration minor axis", i.minor * 1000, "\n")  
  i.major = (I.major / Area)^0.5
  #cat("rad gyration major axis", i.major * 1000, "\n") 
  EoverG = E/G
  
  # Flexural

  Cminor = NULL
  Cmajor = NULL
  flex.minor.Hz = NULL
  flex.major.Hz = NULL

  for (n in 1:mode.count[1]){
  Cminor[n] = 0.5+0.5*(K[n]^2)*((i.minor/L)^2)*(1+(form.factor.minor)*EoverG)+(0.25+0.5*(K[n]^2)*((i.minor/L)^2)*(1+(form.factor.minor)*EoverG)+(0.5*(K[n]^2)*((i.minor/L)^2)*(1-(form.factor.minor)*EoverG))^2)^0.5
  Cmajor[n] = 0.5+0.5*(K[n]^2)*((i.major/L)^2)*(1+(form.factor.major)*EoverG)+(0.25+0.5*(K[n]^2)*((i.major/L)^2)*(1+(form.factor.major)*EoverG)+(0.5*(K[n]^2)*((i.major/L)^2)*(1-(form.factor.major)*EoverG))^2)^0.5
  flex.minor.Hz[n] = (E/(4*(L^4)*density*Cminor[n]*(pi^2)/((K[n]^4)*(i.minor^2))))^0.5
  flex.major.Hz[n] = (E/(4*(L^4)*density*Cmajor[n]*(pi^2)/((K[n]^4)*(i.major^2))))^0.5
  }

  # Torsional
  
  a. = b/2 # half the width
  b. = h/2 # half the thickness
  s = a. / b. # aspect ratio
  g1 = 1 # 1 when Gxy = Gxz
  d = s*g1
  Jp = (4*s*b.^4)*(1+s^2)/3 # polar second moment of cross-sectional area
  Jg = (4*s*b.^4)*(1+d^2)/3 # modified moment
  
  
  # (calculate the equivalent of (1+A)/B in ASTM E1876-1  
  kg.sum = 0
  for(iterate.i in 1:10) {
    iterate.n = iterate.i*2-1 
    kg.sum = kg.sum + (1/iterate.n^5)*tanh(0.5*pi*iterate.n*d)
  }
  kg = (4/(1+d^2))*(1-(192/(pi^5*d))*kg.sum) 

  # calculate mu squared 
  mu.sum = 0
  for(iterate.i in 1:10) {
    iterate.n = iterate.i*2-1 
    mu.sum = mu.sum + (1/iterate.n^6)*((6*tanh(0.5*pi*iterate.n*d))/(pi*iterate.n*d)+(tanh(0.5*pi*iterate.n*d))^2-3)
  }
  mu.squared = ((b.^2)/(1+d^2))*((s^2)/3 + mu.sum*(768)/(pi^6*g1^2))
  
  C = kg*G*Jg
  kp = kg

  torsional.Hz = NULL
    for (tor.i in 1:mode.count[2]) {
    pc = ((tor.i * pi)/L)*(kp*G/density)^0.5 # circular frequency
    mu = (mu.squared)^0.5
    h.tor = tor.i*pi*mu/L
    Aquad = 1
    Bquad = - ( (E/(G*kp)) + eta.squared*(Jg/(Jp*kp)) + (1-kg)/(kp*h.tor^2)  )
    Cquad = + (Jg/Jp) * (eta.squared*E/(G*kp^2) + (kg*(1-kg))/((h.tor^2)*(kp^2))) 
    N.solution1 = (-Bquad + sqrt(Bquad^2 - 4*Aquad*Cquad)) / (2*Aquad)
    N.solution2 = (-Bquad - sqrt(Bquad^2 - 4*Aquad*Cquad)) / (2*Aquad)
    N.solution = min(N.solution1, N.solution2)
    torsional.Hz[tor.i] = pc / (2*pi) # convert to Hz
    torsional.Hz[tor.i] = N.solution^0.5 * torsional.Hz[tor.i]
    }
  
  if(torsion.fast){
    cat("Using the approximate method for torsional calculation - this is not so good for high E/G, non-slender specimens and high mode numbers. \n")
    cat("Change this with torsion.fast = FALSE \n")  
  }else{
    cat("Using the more accurate method for the torsional calculation, this might take a while \n")
    cat("Change this with torsion.fast = TRUE \n")
    
    solve.this.even = function(f){
    p = f*2*pi  
    a1 = 0.5*(L^2)*(((density*p^2)/G)*(Jp/(Jg*eta.squared)+G/E)-((kg*(1-kg))/(eta.squared*mu.squared))*(G/E))
    a2 = -(L^4)*((((1-kg)*(density*Jp*p^2)  )/(eta.squared*mu.squared*E*Jg)) - (((density*p^2)^2)/(eta.squared*E*G))*(Jp/Jg))
    alpha = (-a1+(a1^2-a2)^0.5)^0.5
    beta =  (a1+(a1^2-a2)^0.5)^0.5
    b1 = (1-(density*mu.squared*p^2)/((G*(1-kg))))^-1
    b2 = density*p^2*L^2*Jp / (G*Jg)
    a3 = (b1*L^2)*(C - (1+E/G)*(((p^2)*density*mu.squared*Jg)/(1-kg)))
    a4 = E*mu.squared*Jg*b1*Jg / ((1-kg)*Jp)
    A1 = alpha* (a3-a4*alpha^2)*(b2- beta^2)
    A2 = beta * (a3+a4* beta^2)*(b2+alpha^2)
    a12 = A1/A2
    freqeq1 = tan(0.5*beta)+a12*tanh(0.5*alpha)
    return(freqeq1)
    }

    solve.this.odd = function(f){
    p = f*2*pi  
    a1 = 0.5*(L^2)*(((density*p^2)/G)*(Jp/(Jg*eta.squared)+G/E)-((kg*(1-kg))/(eta.squared*mu.squared))*(G/E))
    a2 = -(L^4)*((((1-kg)*(density*Jp*p^2)  )/(eta.squared*mu.squared*E*Jg)) - (((density*p^2)^2)/(eta.squared*E*G))*(Jp/Jg))
    alpha = (-a1+(a1^2-a2)^0.5)^0.5
    beta =  (a1+(a1^2-a2)^0.5)^0.5
    b1 = (1-(density*mu.squared*p^2)/((G*(1-kg))))^-1
    b2 = density*p^2*L^2*Jp / (G*Jg)
    a3 = (b1*L^2)*(C - (1+E/G)*(((p^2)*density*mu.squared*Jg)/(1-kg)))
    a4 = E*mu.squared*Jg*b1*Jg / ((1-kg)*Jp)
    A1 = alpha* (a3-a4*alpha^2)*(b2- beta^2)
    A2 = beta * (a3+a4* beta^2)*(b2+alpha^2)
    a12 = A1/A2
    freqeq2 = tan(0.5*beta)-(1/a12)*tanh(0.5*alpha)
    return(freqeq2)
    }    
    
    gap = 0.4 * (torsional.Hz[2] - torsional.Hz[1]) # the +/- range to search
    
  for(do.mode in 1:mode.count[2]){
    
  lower.range = torsional.Hz[do.mode] - gap
  upper.range = torsional.Hz[do.mode] + gap

  if(do.mode %% 2 == 0){
  solved = uniroot(solve.this.even, lower = lower.range, upper = upper.range)
  torsional.Hz[do.mode] = solved$root
  }else{
  
  # this is a little tricky because of the discontinuous part so we need to figure out the range 
  
  try.counts = c(10, 10^3, 10^6, 10^9, 10^12, 10^15)
  for (try.i in 1:length(try.counts)){
  
  if(try.i == 3) cat("The calculation cannot go at the fastest rate \n")  
  if(try.i == 4) cat("The calculation needs a lot of time for the",do.mode,"th torsional mode \n")
  if(try.i == 5) cat("Go and do something else... this will take a long time \n")
  if(try.i == 6) cat("This is not going to finish any time soon \n")
    
  ftry = seq(from = lower.range, to = upper.range, length.out = try.counts[try.i])
  ftry.result = solve.this.odd(ftry)
  ftry.result.nan = na.omit(max(ftry.result))
  new.upper.range = max(ftry.result.nan)
  if(new.upper.range > 0) break
  }
  
  solved = uniroot(solve.this.odd, lower = lower.range, upper = ftry[which(ftry.result == new.upper.range)])
  
  
  torsional.Hz[do.mode] = solved$root
  }  
    
  }  
  }


  # Longitudinal
  
  correction.T = 1 + ((pi^2)/(Area*L^2))*(Poisson.ratio^2 * I.major + Poisson.ratio^2 * I.minor)
  
  longitudinal.Hz = seq(from = 1, to = mode.count[3])
  longitudinal.Hz = longitudinal.Hz * (E / ((2*L)^2 * density * correction.T))^0.5

  result = NULL
  result$flex.minor.Hz = flex.minor.Hz
  result$flex.major.Hz = flex.major.Hz
  result$torsional.Hz = torsional.Hz
  result$longitudinal.Hz = longitudinal.Hz
  
  all.freq = c(flex.minor.Hz, flex.major.Hz, torsional.Hz, longitudinal.Hz)
  names(all.freq) = c(rep("Flexural minor", length.out = mode.count[1]),
                      rep("Flexural major", length.out = mode.count[1]),
                      rep("Torsional", length.out = mode.count[2]),
                      rep("Longitudinal", length.out = mode.count[3]))
  all.freq = all.freq[order(all.freq)]
  result$all.Hz = all.freq
  
  result$values = c(L*1000, b*1000, h*1000, density, E/10^9, G/10^9, Poisson.ratio, form.factor.minor, form.factor.major, eta.squared)
  names(result$values) = c("Length (mm)", "Width (mm)", "Thickness (mm)", "Density (kg/m^3)", "E (kN/mm^2)", "G (kN/mm^2)", "Poisson's ratio", "Form factor minor", "Form factor major", "eta squared")

  return(result)
  
} # End of function

