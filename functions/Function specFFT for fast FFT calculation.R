
# Function to calculate FFT spectrum of an impact excitation recording
# By Dan Ridley-Ellis
# Last updated 27/11/2017
#
# Args:
# wave = the wave to be analysed (must contain sampling rate)
# padding = the wave is zero padded to a multiple of this power of 2 (forced to be an integer, at least 1)
# plot = TRUE (default) or FALSE : Whether graphs are plotted or not
#
# Returns:
# An FFT frequency spectrum
#
# Notes:
# Based on the spec function of package "seewave", simplified and corrected for frequency accuracy
# The input wave must be TuneR Wave class, mono, and have sample rate information included
# The FFT is calculated with a rectangular window over the full duration, and dB as max 0
# The wave is zero padded to the next multiple of 2^9 for faster computation (can be configured)

      specFFT = function (wave, padding = 9, plot = "TRUE") {
       
        input = inputw(wave = wave)
        wave = input$w
        f = input$f
        rm(input)
        padding = max(round(padding,0),1)
        
        zeropadding = 2^padding - length(wave)%%(2^padding)
        
        wavepadded = c(wave, rep(0,zeropadding))
        
        y = Mod(fft(wavepadded))
        
        y = head(y,n = length(y)/2)
        
        y = 20 * log10(y/max(y))
        
        x = seq(0, (f/2), length.out = length(y) + 1)/1000 
        x = head(x, n = length(x)-1)
        
        if (plot == TRUE) {
          plot(x, y, type="l", xlab = "Frequency (kHz)", ylab = "Amplitude (dB)")
        }

        spec = cbind(x, y)
        return(spec)
        
      }

      
      
      
      
      
      