
# Function for finding peaks and bandwidth Q using an fft spectrum (in dB)
# By Dan Ridley-Ellis
# Last updated 24/09/2017
#
# Requires r packages "tuneR" and "seewave"
#
# Args:
# spectrum = fft spectrum in dB, with frequency in first column and dB in second column
# minf = lower end of band in which to search the spectrum (kHz) 
# maxf = upper end of band in which to search the spectrum (kHz) 
# minQ = lowest Q value of interest (peaks with lower Q are ignored) 
# maxQ = highest Q value of interest (peaks with higher Q are ignored) 
# mindB = required amplitude for peak to be recorded 
# peakdist = peaks are only recorded if there is no taller peak within this distance (kHz)
# mindBrelative = required amplitude above average across peakdist for peak to be recorded  
# level = dB difference on which to calculate bandwidth Q (standard is -3  dB)
# minpoints = lower limit on absolute peak width (in spectral points) - use to avoid spurious peaks
# peakscount = maximum number of peaks to be reported (lower frequency peaks first)
# peakfitpoints = number of points involved in the peak top fitting spline (odd and at least 5)
# plot = TRUE (default) or FALSE : Whether graphs are plotted or not
# 
# Returns
# A matrix of peak frequencies and damping (as bandwidth Q, and decay constant, lambda) 
# Results are reported according to two methods (see below)
# For each peak found:
#   f.simple.Hz
#   f.spline.Hz
#   peak.amplitude.simple.dB
#   Q.simple.peak
#   Q.simple.geom
#   Q.spline.peak
#   Q.spline.geom
#   lambda.simple.peak
#   lambda.simple.geom
#   lambda.spline.peak
#   lambda.spline.geom
# If no peaks are found, the function returns one row of zeros  
#   
# Notes:
# There are two methods in use:
# A) peaks correspond to fft spectrum points, and Q is calculated using linear interpolation
# B) peaks and Q are calculated from a spline fitted through spectral points
# Q is also calculated on the basis of the peak frequency and the geometric mean frequency
# Q and lamba are related (lamba = pi*f/Q)
# Units for lambda are /s
# 
# Noisy spectra may result in many peaks recorded, so use the parameters above to target desired peaks
# Note that the criteria for minf, maxf, minQ and maxQ are applied on method A calculation
# This version uses nested loops. Faster calculation may be achieved using apply
# If an exponential window has previously been applied on the wave, the lambda and Q will need to be corrected

peaksQ = function(spectrum,
                minf,
                maxf,
                minQ = 1,
                maxQ = 10000,
                mindB,
                peakdist,
                mindBrelative,
                level = -3,
                minpoints = 3,
                peakscount,
                peakfitpoints = 9,
                plot = TRUE) {
  
  # compute the range of frequency and amplitude values (kHz in column 1, amplitude in column 2)
  range = matrix(c(apply(spectrum, 2, FUN = min ), apply(spectrum, 2, FUN=max)), nrow = 2, ncol = 2, byrow=TRUE)
  if(missing(minf) == TRUE) minf = (range[1, 1])
  if(missing(maxf) == TRUE) maxf = (range[2, 1])
  if(missing(mindB) == TRUE) mindB = (range[1, 2])
  
  # error messages if data is not correct
  if (range[2, 2] == 1) stop("data must be in dB")
  if (minQ > maxQ) stop("minQ must be less than maxQ")
  if (minf > maxf) stop("minf must be less than maxf")
  if(missing(mindBrelative) == FALSE & missing(peakdist) == TRUE) stop("mindBrelative requires peakdist")
  
  # begin by making a copy of the spectrum and setting up the results container
  spec = spectrum
 
  results=matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow=1,ncol=11) # see end of script for column names
  
  # if minimum peak distance is active determine the number of points this represents
  if(missing(peakdist) != TRUE) fstep = ceiling(peakdist / (spec[2, 1] - spec[1, 1])) else fstep = -1
  
  # the top level loop, examines each point in the spectrum (with two points padding at each end)
  for(i in 1:(length(spectrum[, 1]) - 4)) {
    f.point = i + 2
    
    # make a note of the amplitude at the point being examined
    amp.ref = spec[f.point, 2]
    
    # set state
    go=TRUE
    
    fmin = -1 # assigns a temporary value of -1 
    f.point.relative = 1
    
    # some quick checks to increase speed
    # don't bother if the frequency is too low or too high
    if(spec[f.point, 1] < minf) {go = FALSE}
    if(spec[f.point, 1] > maxf) {go = FALSE}
    # don't bother if the peak is not high enough
    if((amp.ref < mindB) == 1) {go = FALSE}
    # don't bother if the neighbouring points are not both lower
    if(spec[(f.point - 1), 2] > amp.ref) {go = FALSE}
    if(spec[(f.point + 1), 2] > amp.ref) {go = FALSE}
    # don't bother if enough peaks are already found
    if(missing(peakscount) == FALSE) if(results[1] != 0) if(peakscount <= length(results[, 1])) go = FALSE
    
    # The second level loop, which runs only if the go state is TRUE
    # Searches left from this point for a crossing point for bandwidth Q
    while (go == TRUE) {
      
      # stop if the spectrum goes too far left
      if((f.point - f.point.relative) == 1) {go = FALSE}
      # stop if the spectrum goes up beyond the reference amplitude (can't be a peak)
      if(spec[(f.point - f.point.relative), 2] > amp.ref) {go = FALSE}
      
      # stop if the spectrum goes below the reference level (crossing point found)
      if(spec[(f.point - f.point.relative), 2] < (amp.ref + level)) {
        fmin = spec[(f.point - f.point.relative), 1]
        
        #calculate the lower frequency of band (linear interpolation)
        fmin.f.left = spec[(f.point - f.point.relative), 1]
        fmin.f.right = spec[(f.point - f.point.relative + 1), 1]
        fmin.a.left = spec[(f.point - f.point.relative), 2] - amp.ref
        fmin.a.right = spec[(f.point - f.point.relative + 1), 2] - amp.ref
        fmin = fmin.f.left + (fmin.f.right - fmin.f.left) * ((fmin.a.left - level)/(fmin.a.left - fmin.a.right))
        points.in.peak = f.point.relative
        Qfpoint.left = f.point.relative
        
        go = FALSE 
      }
      
      f.point.relative = f.point.relative + 1
      # end of searching left
    }
    
    # only bother to search right if the search left found a crossing point
    if(fmin > 0) {
      # search right
      # set state
      go = TRUE
      fmax = -1 # assigns a temporary value of -1 
      f.point.relative = 1
      
      # Another second level loop, which runs only if the go state is TRUE
      while (go == TRUE){
        
        # stop if the spectrum goes too far right
        if((f.point + f.point.relative) == (length(spectrum[, 1]))) {go = FALSE}
        
        # stop if the spectrum goes up beyond the reference amplitude
        if(spec[(f.point + f.point.relative), 2] > amp.ref) {go = FALSE}
        
        # stop if the spectrum goes below the reference level - crossing point found
        if(spec[(f.point + f.point.relative), 2] < (amp.ref + level)) {
          
          # a peak has been found. Work out the peak and the Q value
          Qfpoint.right = f.point.relative
          # calculate the width of the peak in spectrum points
          points.in.peak = points.in.peak + f.point.relative
          
          #calculate the upper frequency of band (linear interpolation)
          fmax.f.left = spec[(f.point + f.point.relative - 1), 1]
          fmax.f.right = spec[(f.point + f.point.relative), 1]
          fmax.a.left = spec[(f.point + f.point.relative - 1), 2] - amp.ref
          fmax.a.right = spec[(f.point + f.point.relative), 2] - amp.ref
          fmax = fmax.f.left + (fmax.f.right - fmax.f.left) * ((fmax.a.left - level)/(fmax.a.left - fmax.a.right))
          
          # peak frequency is the point being examined
          fcen = spec[f.point, 1]
          
          # compute results (method A)
          Q = (fcen / (fmax - fmin))
          fcentre = sqrt(fmax * fmin) # geometric mean
          
          f.simple.Hz = fcen * 1000
          Q.simple.peak = Q
          Q.simple.geom = (fcentre / (fmax - fmin))
          lambda.simple.peak = pi * fcen * 1000 / Q.simple.peak
          lambda.simple.geom = pi * fcen * 1000 / Q.simple.geom
          
          # Proceed further only if this peak passes criteria for recording
          record = 1
          if(Q < minQ) record = 0
          if(Q > maxQ) record = 0
          if(points.in.peak < minpoints) record = 0
          
          # if min peakdist is active, and the peak is not the tallest, don't record
          if(fstep > 0 & go == TRUE) {
            f.max.from = max((f.point - fstep), min(which((spectrum[, 1]) >= minf)))
            f.max.to = min((f.point + fstep), min(which((spectrum[, 1]) >= maxf)))
            if(amp.ref < max(spec[f.max.from:f.max.to, 2])) record = 0
            if(missing(mindBrelative) == FALSE) if(amp.ref < (mean(spec[f.max.from:f.max.to, 2])) + mindBrelative) {record = 0}
          }
          
          # round the peak amplitude
          peak.amplitude.simple.dB = round(amp.ref, 2)
          
          if(record == 1){
            
            
            Peakpoints = spec[(f.point - Qfpoint.left - 2):(f.point + Qfpoint.right + 2), 1:2]
            Qlowf = spec[(f.point - Qfpoint.left), 1]
            Qhighf = spec[(f.point + Qfpoint.right), 1]
            
            # Display this peak in detail
            if (plot == TRUE) plot(Peakpoints, xlab = "Frequency (kHz)", ylab = "Amplitude (dB)")
            if (plot == TRUE) lines(Peakpoints, col = "gray50")
            bandwidth = matrix(c(fmin, fmax, (amp.ref + level),(amp.ref + level)), nrow = 2, ncol = 2, byrow = FALSE)
            if (plot == TRUE) lines(bandwidth, col = "blue")
            if (plot == TRUE) points(f.simple.Hz / 1000, amp.ref, col = "blue", pch = 22)
            if (plot == TRUE) abline(v = f.simple.Hz / 1000, col="blue")
            
            # interpolate with a spline to find a better value of the peak frequency
            
            Peakspline = spline(Peakpoints, n = 10000) # number of points in the spline could be modified for speed vs precision
            if (plot == TRUE) lines(Peakspline, col = "green")
            
            # compute better results
            
            fcentre = Peakspline$x[which(Peakspline$y == max(Peakspline$y))] # a better indicator of peak frequency
            amp.ref = max(Peakspline$y) # update the peak amplitude
            
            Peakspline$y=(Peakspline$y - amp.ref-level)^2
            
            fmin = Peakspline$x[which(Peakspline$y == min((Peakspline$y[which(Peakspline$x < fcentre & Peakspline$x > Qlowf)])))]
            
            fmax = Peakspline$x[which(Peakspline$y == min((Peakspline$y[which(Peakspline$x > fcentre & Peakspline$x < Qhighf)])))]
            
            bandwidth = matrix(c(fmin, fmax, (amp.ref + level),(amp.ref + level)),nrow = 2,ncol = 2, byrow = FALSE)
            if (plot == TRUE) lines(bandwidth, col="red")
            if (plot == TRUE) points(fcentre, amp.ref, col = "red")
            if (plot == TRUE) abline(v = fcentre,col = "red")
            if (plot == TRUE) if(fstep > 0) lines(spec[f.max.from:f.max.to,1], spec[f.max.from:f.max.to, 2],col = "black")
            
            f.spline.Hz = fcentre * 1000
            Q.spline.peak = fcentre / (fmax-fmin)
            Q.spline.geom = sqrt(fmax * fmin) / (fmax-fmin)
            lambda.spline.peak = pi * fcentre * 1000 / Q.spline.peak
            lambda.spline.geom = pi * fcentre * 1000 / Q.spline.geom
            
            # add to the results matrix
            if(results[1] == 0) {
              results = matrix(c(f.simple.Hz,
                               f.spline.Hz,
                               peak.amplitude.simple.dB,
                               Q.simple.peak,
                               Q.simple.geom,
                               Q.spline.peak,
                               Q.spline.geom,
                               lambda.simple.peak,
                               lambda.simple.geom,
                               lambda.spline.peak,
                               lambda.spline.geom),
                               nrow = 1, ncol = 11)}
              else {results = rbind(results,c(f.simple.Hz,
                                            f.spline.Hz,
                                            peak.amplitude.simple.dB,
                                            Q.simple.peak,
                                            Q.simple.geom,
                                            Q.spline.peak,
                                            Q.spline.geom,
                                            lambda.simple.peak,
                                            lambda.simple.geom,
                                            lambda.spline.peak,
                                            lambda.spline.geom))}
            
          }
          
          go = FALSE
          
          # end of what to do if peak is found
        }
        
        f.point.relative = f.point.relative + 1
      }
      
    }
    # end of searching through points on the spectrum  
  }
  
  # plot the spectrum
  if (plot == TRUE) plot(spec,type = "l",col = "cornflowerblue", xlab = "Frequency (kHz)",ylab = "Amplitude (dB)", xlim = c(minf, maxf))
  if (plot == TRUE) if(mindB != (range[1,2])) abline(h = mindB, col = "grey")
  if (plot == TRUE) if(minf != (range[1,1])) abline(v = minf, col = "grey")
  if (plot == TRUE) if(maxf != (range[2,1])) abline(v = maxf, col = "grey")
  
  if (plot == TRUE) points(results[, 2] / 1000, results[, 3], col = "red")
  if (plot == TRUE) abline(v = results[, 2] / 1000, col = "red")
  
  colnames(results)=c(
    "f.simple.Hz",
    "f.spline.Hz",
    "peak.amplitude.simple.dB",
    "Q.simple.peak",
    "Q.simple.geom",
    "Q.spline.peak",
    "Q.spline.geom",
    "lambda.simple.peak",
    "lambda.simple.geom",
    "lambda.spline.peak",
    "lambda.spline.geom")
  
    return(results)
  
  # end of function
}