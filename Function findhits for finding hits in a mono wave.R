
# Function to detect hit start and end times in a mono wave 
# By Dan Ridley-Ellis
# Last updated 25/09/2017
#
# Requires r packages "tuneR" and built for use with package "seewave"
#
# Args:
# wave = the wave to be analysed (must contain sampling rate)
# threshold = the amplitude % considered to be a hit
# silence = the amplitude % below which is considered a hit can start and end
# span = parameter for smoothing the wave envelope (width of sliding window in seconds)
# plot = TRUE (default) or FALSE : Whether graphs are plotted or not (if no hits are found, a graph is drawn for diagnostics)
#
# Returns:
# A matrix with start and end times, a row for each hit
#   h.start = time at which a hit starts (s)
#   h.end = time at which a hit ends (s)
#   h.duration = length of the hit (s)
#   if no hits are found, the function returns NULL
#
# Notes:
# The input wave must be TuneR Wave class, mono, and have sample rate information included
# If there is any zero offset on the wave, this should have first been removed

findhits = function(wave = NULL,
                    threshold = 20,
                    silence = 5,
                    span = 0.0005,
                    plot = TRUE) {

  # Make sure that the correct values are specified
  
  if (is.null(wave)) 
    stop("you must specify 'wave'")

  # gather information from the wave file
  samplef = wave@samp.rate
  samplebit = wave@bit
  amplitude = as.matrix(wave@left)
  amplitude.scale = max(abs(amplitude))
  amplitude = 100 * amplitude / amplitude.scale
  
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
  
  wave.env.fitted = wave.env
  
  # remove any points that lie between silence and the threshold, since a hit is not allowed to start or end here
  wave.env = wave.env[which(wave.env[,2] >= threshold | wave.env[,2] <= silence),]
  
  # a hit is when the envelope goes from below the theshold to above threshold
  wave.hits = apply(as.matrix(1:(length(wave.env[, 1]) - 1)), 1, function(x) ifelse(wave.env[x,2] <= threshold & wave.env[x+1,2] > threshold, yes = 1, no = 0))  
  
  h.start = wave.env[which(wave.hits == 1), 1]
  
  # only proceed if some hits are found
  if(length(h.start) > 0) {

  h.end = tail(h.start, n = length(h.start) - 1)
  h.end = c(h.end, max(wave.env[, 1]))
  h.duration = h.end - h.start

  if(plot == "TRUE"){
  # plot the results
  plot(t, amplitude, type="l", ylab = "Amplitude (%)", xlab = "Time (s)", col = "gray80", ylim=c(-100, 100))
  lines(wave.env.fitted, col="black")
  lines(wave.env, col="blue")
  abline(h = silence, col = "gray80")
  abline(h = threshold, col = "gray80")
  abline(v = h.start, col = "gray50")
  points(h.start, seq(from = threshold, to = threshold, length.out = length(h.start)) , col="green")
  points(h.end, seq(from = silence, to = silence, length.out = length(h.start)) , col="red")
    }
  
  colnames(wave.env) = c("time (s)",
                         "amplitude envelope (%)")
  
  results = cbind(h.end, h.duration)
  results = cbind(h.start, results)
  colnames(results) = c("h.start (s)",
                        "h.end (s)",
                        "h.duration (s)")
  
  return(results)

 }  
 
  # if no hits are found
  print("no hits were found")
  
  # plot the results for diagnostic
  plot(t, amplitude, type="l", ylab = "Amplitude (%)", xlab = "Time (s)", col = "gray80", ylim=c(-100, 100))
  lines(wave.env.fitted, col="black")
  lines(wave.env, col="blue")
  abline(h = silence, col = "gray80")
  abline(h = threshold, col = "gray80")
  
  return(NULL)
 }

  