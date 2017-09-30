
# Function to simulate an impact excitation sound recording
# By Dan Ridley-Ellis
# Last updated 23/09/2017
#
# Requires r package "tuneR" and built for use with package "seewave"
#
# Args:
# samplef = sampling frequency in kHz
# e.g.
#     22.05 (maximum frequency is 10 kHz)
#     44.1 (audio CD, maximum frequency is 22 kHz, default) 
#     48.0 
#     88.2
#     96.0 (maximum frequency is 43.6 kHz)
#
# samplebit = sampling bit depth
# e.g.
#     8 (low quality)
#     16 (audio CD, default) 
#     24 
#     32 
#
# duration = total recording duration in seconds
# or
# decayed = % : the duration is set as the point at which the vibration has decayed to this much (max is 100)
# start = time at which the vibration begins
# f = frequency of the vibration in Hz
# a = relative amplitude (of the vibration's envelope)
# Q = bandwidth Q factor
# or
# lambda = decay constant in /s
#
# Returns:
# A Wave object
#
# Notes:
# There must be an equal number of start times, frequencies, relative amplitudes, and damping values (Q or lambda)
# You cannot specify both duration and decayed
# You cannot specify both Q and lamba (lamba = pi*f/Q)
 

simhit = function(samplef = 44.1, 
                  samplebit = 16,
                  duration = NULL,
                  decayed = NULL,
                  start = NULL,
                  f = NULL,
                  a = NULL,
                  Q = NULL,
                  lambda= NULL) {
  
  # Make sure that the correct values are specified
  
  if (!is.null(duration) && !is.null(decayed)) 
    stop("'duration' and 'decayed' cannot be used together")
  if (!is.null(Q) && !is.null(lambda)) 
    stop("'Q' and 'lambda' cannot be used together")
  if (is.null(duration) && is.null(decayed)) 
    stop("You must specify either 'duration' or 'decayed'")
  if (is.null(Q) && is.null(lambda)) 
    stop("You must specify either 'Q' or 'lambda'")
  if (is.null(start)) 
    stop("You must specify 'start'")
  if (is.null(f)) 
    stop("You must specify 'f'")
  if (is.null(a)) 
    stop("You must specify 'a'")
  
  if (is.null(Q)) {
    count = length(lambda)}
    else {count = length(Q)}
    
  if (max(count, length(start), length(f), length(a)) != min(count, length(start), length(f), length(a))) 
    stop("There must be an equal number of start, f, a and damping (Q or labda) values ")
  
  
  # if Q is specified, convert this to lambda
  if (is.null(lambda)) {
    lambda = pi * f / Q 
  }
  
  # if decayed is specified, calculate the duration
  if (is.null(duration)) {
    duration = max(start + log(100 / decayed) / lambda)
  }
    
  # create a time vector on which to calculate the waves
  samplef = samplef * 1000
  t = seq(0, duration, (1/samplef))
  
  # create an initial amplitude vector
  amplitude = 0 * t
  
  # make a function to calculate amplitude at time t
  AmplitudeCalc = function(start, a , lambda, f, t) {
    
    started = t-start
    started = sapply(started, function(x) if(x < 0) {0} else {1})
    amplitude = sum(started * a * exp(-lambda*(t - start)) * sin(2 * pi * f * (t - start)))  
  
    return(amplitude)
  }
  
  # calculate amplitude across all values of t
  amplitude = sapply(t, function(x) AmplitudeCalc(start = start, a = a, lambda = lambda, f = f, t = x))
  
  # make the wave object
  wav = Wave(left = amplitude, samp.rate=samplef, bit=samplebit)
  
  return(wav)
  
}

  
