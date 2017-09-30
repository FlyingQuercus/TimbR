
# Function to apply an exponential window on a mono wave
# By Dan Ridley-Ellis
# Last updated 24/09/2017
#
# Requires r package "tuneR" and built for use with package "seewave"
#
# Args:
# wave = the wave to be modified (must contain sampling rate)
# peakt = the time at which the window is maximum (default is the start of the wave)
# decayin = the decay constant for the initial part of the hit
# decayout = the decay constant for the vibration response
# plot = TRUE (default) or FALSE : Whether graphs are plotted or not
#
# Returns:
# A Wave object with the window applied
#
# Notes:
# The input wave must be TuneR Wave class, mono, and have sample rate information included

expwin = function(wave = NULL, 
                  peakt = NULL,
                  decayin = NULL,
                  decayout = NULL,
                  plot = TRUE) {
  
  # Make sure that the correct values are specified
  
  if (is.null(wave)) 
    stop("'you must specify 'wave'")
  if (is.null(decayout)) 
    stop("You must specify 'decayout'")
 
  # if no peakt is specified, set it to the start of the wave
  if (is.null(peakt)) {
  peakt = 0
  decayin = 0
  } else {
  if (is.null(decayin)) 
  stop("You must specify 'decayin'")
  }
  
  # gather information from the wave file
  samplef = wave@samp.rate
  samplebit = wave@bit
  amplitude = as.matrix(wave@left)
  peaktpoints = round(peakt * samplef, 0)
  
  
  # the in part of the wave
  amplitudein = head(amplitude, n = peaktpoints) 
  t = seq(from = peakt, length.out = length(amplitudein), by = (-1/samplef))
  amplitudein = amplitudein * exp(-decayin * t)
  
  # the out part of the wave
  amplitudeout = tail(amplitude, n = length(amplitude) - peaktpoints) 
  t = seq(from = 0, length.out = length(amplitudeout), by = (1/samplef))
  amplitudeout = amplitudeout * exp(-decayout * t)
  
  # put the in and out parts back together
  amplitude2 = rbind(as.matrix(amplitudein), as.matrix(amplitudeout))
  
  # make the wave object
  wav = Wave(left = amplitude2, samp.rate=samplef, bit=samplebit)
  
  if(plot == "TRUE"){
  # plot the results
  t = seq(from = 0, length.out = length(amplitude), by = (1/samplef))
  plot(t, amplitude, type="l", ylab = "Amplitude", xlab = "Time (s)", col = "gray80")
  lines(t, amplitude2, col="red")
  abline(v = peakt, col = "blue")
  }
  
  return(wav)
  
}

  
