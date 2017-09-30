
# Function to calculate log decrement damping 
# Original by Melanie Libeau
# With further work by Dan Ridley-Ellis
# Last updated 25/09/2017
#
# Requires r packages "tuneR" and "seewave"
#
# Args:
# wave = the wave to be analysed (must contain sampling rate)
# Pmax = percentage of max of the peak the high value (linear model)
# Pmin = percentage of max of the peak for the low value (linear model)
# span = parameter for smoothing the wave envelope (width of sliding window in seconds)
# plot = TRUE (default) or FALSE : Whether graphs are plotted or not
#
# Returns:
# A matrix of results 
#   lambda.logdecrement
#   r.squared
#
# Notes:
# The input wave must be TuneR Wave class, mono, and have sample rate information included
# If there are several component waves, the wave should first be filtered for the frequency of interest
# If an exponential window has previously been applied on the wave, the lambda will need to be corrected
# If there is any zero offset on the wave, this should have first been removed

logdecrement = function(wave = NULL,
                        Pmax = 90,
                        Pmin = 10,
                        span = 1/40000,
                        plot = TRUE) {

  # Make sure that the correct values are specified
  
  if (is.null(wave)) 
    stop("you must specify 'wave'")

  # gather information from the wave file
  samplef = wave@samp.rate
  samplebit = wave@bit
  amplitude = as.matrix(wave@left)
  
  # make the corresponding time vector
  t = seq(from = 0, length.out = length(amplitude), by = (1/samplef))
  

  # the number of sampling points over which to make the envelope
  envelope.band.points.half = round(samplef * span / 2, 0) 
  
  # make an envelope of the amplitude
  n = length(amplitude) # the number of sampling points (time domain)
  
  amplitude.abs = abs(amplitude) # max and min amplitude count equally (if there is no zero offset)
  
  # find the maximum amplitude points within the span
  amplitude.env = sapply(seq(1:n), function(x) if(amplitude.abs[x] == max(amplitude.abs[(max(x - envelope.band.points.half, 0)):(min(x + envelope.band.points.half, n))])) {amplitude.abs[x]} else {NA})
  
  # make a matrix with time and amplitude
  wave.env = cbind(t, amplitude.env)
  # remove all the points that are not maximums
  wave.env = na.omit(wave.env)
  
  M = which.max(wave.env[,2])
  # extract the part of the envelope of interest
  wave.env2 = wave.env[M:length(wave.env[,2]),]
  
  wave.env2[,2] = sapply(wave.env2[,2], function(x) if(x < Pmax / 100 * max(wave.env[,2]) & x > Pmin / 100 * max(wave.env[,2])) {x} else {NA})
  
  wave.env2 = na.omit(wave.env2)
  
  if(plot == "TRUE"){
  # plot the results
  plot(t, amplitude, type="l", ylab = "Amplitude", xlab = "Time (s)", col = "gray80")
  lines(wave.env, col="red")
  lines(wave.env2, col="blue")
  points(wave.env[M,1], wave.env[M,2], col="red")
  abline(h = Pmax / 100 * max(wave.env[,2]), col = "gray80")
  abline(h = Pmin / 100 * max(wave.env[,2]), col = "gray80")
  }
  
  
  wave.env2[,2] = log(wave.env2[,2])
    
  
  # Make a linear model for decrement calculation
  wav.env2.lm = lm(wave.env2[,2] ~ wave.env2[,1])
  
  lambda.logdecrement = wav.env2.lm$coefficients[2] # gives the slope value of the linear model
  
  r.squared = summary(wav.env2.lm)$r.squared # gives the correlation coefficient

  
  if(plot == "TRUE"){  
  plot(wave.env2, ylab = "Natural log of amplitude", xlab = "Time (s)", col = "gray80")
  abline(wav.env2.lm, col="red") # show the line from the linear model 
  }

  
  results = matrix(c(lambda.logdecrement,
                   r.squared),
                   nrow = 1, ncol = 2)
  
  colnames(results)=c("lambda.logdecrement",
                      "r.squared")

  return(results)

 }

  