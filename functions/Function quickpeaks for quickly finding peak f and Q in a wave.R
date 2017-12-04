
# Function giving a simple way to get peak frequencies and damping (from bandwidth as lambda and Q) 
# from a wave containing one or more hits
# By Dan Ridley-Ellis
# Last updated 27/11/2017
#
# Requires r packages "tuneR" and "seewave"
#
# Args:
# ##### The basic input
# wave = the wave to be analysed (must contain sampling rate)
# ID = an ID lable to attach to the results in the output
# ##### Related to finding the hits
# threshold = the amplitude % considered to be a hit
# silence = the amplitude % below which is considered a hit can start and end
# longesthit = the max allowed hit duration
# span = parameter for smoothing the wave envelope (width of sliding window in seconds)
# maxhitcount = the maximum number of hits to analyse (counting from the start)
# plothits = TRUE (default) or FALSE : Whether graphs are plotted or not (if no hits are found, a graph is drawn for diagnostics)
# ##### Related to the portion of the wave analysed as a hit, and the exponential window applied
# expwinends = percentage value (of maximum amplitude) that the window reduces to (at least) 
# ##### Related to the finding of peaks
# minf = lower end of band in which to search the spectrum (Hz) 
# maxf = upper end of band in which to search the spectrum (Hz) 
# minQ = lowest Q value of interest (peaks with lower Q are ignored) 
# maxQ = highest Q value of interest (peaks with higher Q are ignored) 
# mindB = required amplitude for peak to be recorded 
# peakdist = peaks are only recorded if there is no taller peak within this distance (kHz)
# mindBrelative = required amplitude above average across peakdist for peak to be recorded  
# minpoints = lower limit on absolute peak width (in spectral points) - use to avoid spurious peaks
# peakscount = maximum number of peaks to be reported (lower frequency peaks first)
# peakfitpoints = number of points involved in the peak top fitting spline (odd and at least 5)
# meanspec = TRUE or FALSE (default): Whether calculation is on mean spectrum or not
# plotpeaks = TRUE (default) or FALSE : Whether graphs are plotted or not
# plotsummary = TRUE (default) or FALSE : Whether summary graphs are plotted or not
# plotsummaryspectra = TRUE or FALSE (default): Whether a (slow) summary graph of spectra is plotted or not
#
# Returns:
# A data frame of results for each hit and peak found
# ID (if not, specified, the system clock time is used)
# Wave sample rate (Hz)
# The hit (as a number, counting from the beginning). Analysis of mean spectrum is given hit number -1 
# The start time of this hit
# The duration of this hit analysed (not counting any zero padding)
# The peak frequency (Hz) (uses spline interpolation)
# The peak amplitude (dB)
# The peak bandwith Q (uses spline interpolation and geometric mean)
# The equivalent damping decay constant lambda (/s) (calculated from Q)
# The version of this script used
# 
# Notes:
# The input wave must be TuneR Wave class, mono, and have sample rate information included
# Zero offset is removed as part of the processing
# An exponential window is applied by default, and damping corrected accordingly
# When no hits, or no peaks, are found zeros are returned

quickpeaks = function(wave = NULL, 
                      ID = NULL,
                      threshold = 20,
                      silence = 5,
                      longesthit = NULL,
                      maxhitcount = NULL,
                      span = 0.0005,
                      plothits = TRUE,
                      expwinends = 1,
                      minf = 20,
                      maxf = NULL,
                      minQ = 1,
                      maxQ = 100000,
                      mindB = NULL,
                      peakdist = NULL,
                      mindBrelative = NULL,  
                      minpoints = 3,
                      peakscount = NULL,
                      peakfitpoints = 9,
                      meanspec = FALSE,
                      plotpeaks = TRUE,
                      plotsummary = TRUE,
                      plotsummaryspectra = FALSE
                      ) {

  # script version information (embedded in results)
  scriptversion = "version_27Nov2017"

  # Get some basic information
  wavesamplerate = wave@samp.rate
  wavesamplebit = wave@bit
  if (is.null(maxf)) maxf = wavesamplerate / 2
  # because peaksQ takes frequencies in kHz
  maxf = maxf / 1000
  minf = minf / 1000 
  if (is.null(ID)) ID = Sys.time()
  
  # prepare the results table
  results = NULL
  spectra = list(1)
  results.header = c("ID",
                     "Wave sample rate (Hz)",
                     "Hit number",
                     "Start time (s)",
                     "Duration (s)",
                     "Peak frequency (Hz)",
                     "Peak amplitude (dB)",
                     "Bandwidth Q",
                     "lambda (/s)",
                     "R script")
  
    
  # validating the input 
  if (is.null(wave)) stop("you must specify 'wave'")
  if (minQ > maxQ) stop("minQ must be less than maxQ")
  if (minf > maxf) stop("minf must be less than maxf")
  if(is.null(mindBrelative) == FALSE & is.null(peakdist) == TRUE) stop("mindBrelative requires peakdist")
  
  
  # Phase 0: preparing the wave and printing some basic information
  
  # Removing any zero offset
  wave =  rmoffset(wave, fun=mean, output="Wave")

  # report some information to the user
  cat("Analysing Wave with sample rate",wavesamplerate,"Hz and bit depth",wavesamplebit,"\n")
  cat("Using R script",scriptversion,"\n")
  cat("The analysis may take some time...\n")
  
  # Phase 1: finding the hits
  cat("Phase 1: Looking for hits...\n")
  
  hits.list = findhits(wave, threshold = threshold, silence = silence, longesthit = longesthit, span = span, plot = plothits) 
  
  # If hits are found, proceed (otherwise just return zeros for results)
  if(is.null(hits.list)) {
    cat("Since no hits were found, the results will be a row of zeros \n")  
    
    results = data.frame(ID,
                   wavesamplerate,
                   0,
                   0,
                   0,
                   0,
                   0,
                   0,
                   0,
                   scriptversion)


  } else {
    number.of.hits = length(hits.list[,1])
    cat("Found",number.of.hits,"hits in this wave\n")
    
    # Phase 2: Searching for peaks
    cat("Phase 2: Looking for peaks...\n")
    
    if(is.null(maxhitcount) == FALSE) {
      hits.analyse = min(length(hits.list[,1]), maxhitcount) } else {
      hits.analyse = length(hits.list[,1])
      }
    
    # If mean spectrum is asked for
    # Work out what padding is required to make the spectra the same length
    if(meanspec == TRUE){
    maxhitlengthpoints = max(hits.list[,3])*wavesamplerate
    padding = ceiling(log(maxhitlengthpoints,2))
    if (padding > 16) {
      cat("WARNING! Spectral length is too high for mean spectrum calculation. Try reducing longesthit\n")
      padding = 9
      meanspec = FALSE
      }
    } else {padding = 9}
    
    
    
    # analyse each hit
    for(i in 1:hits.analyse){
    
      cat("Analysing hit",i,"of",number.of.hits, "found, analysing up to hit", hits.analyse ,"of wave ID:",ID,"\n")  
      # cut out the hit
      hit.wave = cutw(wave,from = max(hits.list[i,1],0),to = hits.list[i,2], output="Wave")
      
      expdamp = log(100 / expwinends) / hits.list[i,3]
      hit.wave = expwin(hit.wave, decayout = expdamp, plot= plothits, decayin = 10000, peakt = span)
      
      # find peaks
      spec = specFFT(hit.wave, plot = plotpeaks, padding = padding)
      peak.results = peaksQ(spec,
                            mindB = mindB,
                            minf = minf,
                            maxf = maxf,
                            minQ = minQ,
                            maxQ = maxQ,
                            peakscount = peakscount,
                            plot = plotpeaks)
      cat("Analysed", length(peak.results[,1]), "peaks\n")
      
      # save the spectrum
      spectra[[i]] = spec

      # collect results 
      calc.f.spline = peak.results[,2]
      calc.peak.amp = peak.results[,3]
      calc.lambda.spline = peak.results[,11] - expdamp
      calc.Q.spline =  pi * calc.f.spline / calc.lambda.spline 
      
      hit.number.col = rep(i, times = length(peak.results[,1]))
      ID.col =  rep(ID, times = length(peak.results[,1]))
      samp.col = rep(wavesamplerate, times = length(peak.results[,1]))
      start.col = rep(hits.list[i,1], times = length(peak.results[,1]))
      duration.col = rep(hits.list[i,3], times = length(peak.results[,1]))
      ver.col = rep(scriptversion, times = length(peak.results[,1]))
      
      results.hit = data.frame(ID.col,
                               samp.col,
                               hit.number.col,
                               start.col,
                               duration.col,
                               calc.f.spline,
                               calc.peak.amp,
                               calc.Q.spline,
                               calc.lambda.spline,
                               ver.col)

      if(is.null(results) == TRUE){
      results = results.hit  
      } else {
      results = rbind(results, results.hit)  
      }
      
      } # end of loop for hits
  }

  
  # If mean spectrum is asked for
  # calculate a mean spectrum and analyse it
  if(meanspec == TRUE){
    cat ("Analysing mean spectrum - shows in results table as hit -1\n")
    meanspectrum = spectra[[1]]
    
    for (k in 1:hits.analyse){
    speck = (spectra[[k]])
    if(k==1) {specamp = speck[,2]} else {specamp = cbind(specamp,speck[,2])}
    }
    
    meanspectrum[,2] = rowMeans(specamp)
    
    # find peaks
    peak.results = peaksQ(meanspectrum,
                          mindB = mindB,
                          minf = minf,
                          maxf = maxf,
                          minQ = minQ,
                          maxQ = maxQ,
                          peakscount = peakscount,
                          plot = plotpeaks)
    cat("Analysed", length(peak.results[,1]), "peaks\n")
    
    # collect results 
    calc.f.spline = peak.results[,2]
    calc.peak.amp = peak.results[,3]
    calc.lambda.spline = peak.results[,11] - expdamp
    calc.Q.spline =  pi * calc.f.spline / calc.lambda.spline 
    
    hit.number.col = rep(-1, times = length(peak.results[,1]))
    ID.col =  rep(ID, times = length(peak.results[,1]))
    samp.col = rep(wavesamplerate, times = length(peak.results[,1]))
    start.col = rep(hits.list[i,1], times = length(peak.results[,1]))
    duration.col = rep(hits.list[i,3], times = length(peak.results[,1]))
    ver.col = rep(scriptversion, times = length(peak.results[,1]))
    
    results.hit = data.frame(ID.col,
                             samp.col,
                             hit.number.col,
                             start.col,
                             duration.col,
                             calc.f.spline,
                             calc.peak.amp,
                             calc.Q.spline,
                             calc.lambda.spline,
                             ver.col)
    
    results = rbind(results, results.hit)  
    
  }
  

  names(results) = results.header

  if((plotsummary == TRUE) & (is.null(hits.list) == FALSE)) {
  # plot peak frequencies against time
  plot(results[,4], results[,6], xlim=c(min(results[,4]), max(results[,4]+results[,5])), ylim=c(minf*1000,maxf*1000), xlab="Time (s)", ylab="Frequency (Hz)", type="p", main = ID)  

  if(plotsummaryspectra == TRUE){
  for(i in 1:hits.analyse){
  spec = spectra[[i]]  
  x = spec[,2]/200 + hits.list[i,1]
  y = spec[,1]*1000
  lines(x, y, col = "grey")
  }}

  plot(results[,6], results[,9], xlab="Frequency (Hz)", ylab="Lambda (/s)", type="p", log = "y", main = ID)      
  }

  if(plotsummaryspectra == TRUE){
  for(i in 1:hits.analyse){
    spec = spectra[[i]]  
    x = spec[,1]*1000
    y = spec[,2]
    if(i == 1){
      plot(x, y, col = "grey", type = "l", xlab = "Frequency (Hz)", ylab = "Amplitude (dB)", xlim=c(minf*1000,maxf*1000))  
    } else {
      lines(x, y, col = "grey")        
    }
    if (meanspec == TRUE) lines(x = meanspectrum[,1]*1000, y = meanspectrum[,2], col = "blue")
    
    points(results[,6], results[,7], col = "red")
    
  }  }
  
  
    
  return(results)
  
  # end of function
  }