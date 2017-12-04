# Functions
This directory contains R functions required for running the scripts.  Look at the comments on each file for details on how to use the functions

## quickpeaks
This may be the only function that you need to use directly.  It is a simple and quick way to obtain peak frequencies and damping (from bandwidth method) from a sound wave containing one or more hits.  The basic input is the wave and an ID label to attach to the results in the output.  There are optional parameters to control the finding of hits and the finding of peaks.  This function calls on some of the functions below, so you will need findhits, expwin, specFFT and peaksQ also loaded.

## simhit
This function creates simulated impulse excitation hits that match theory (sine waves with exponential damping).

## findhits
This function exampines a sound recording and makes a list of times when the hits start and end.  This function cannot tell if there are additional subsequent hits below the trigger threshold, but this can be controlled using the longesthit parameter.

## expwin
This function applies an exponential window (both in and out) to a hit recording, to reduce problems caused by background noise, or the recording ending before the vibration has decayed sufficiently. Any subsequent damping calculations need to be corrected accordingly

## specFFT
This function calculates a frequency spectrum using Fast Fourier Transform (FFT)

## logdecrement
This function calculates damping (decrement, decay constant) on the time domain (only suitable when there is one mode of vibration present. 

## peaksQ
This function calculates peak frequencies and damping (bandwidth Q) on an FFT frequency spectrum (in dB). It can cope with multiple modes present, as long as the peak frequencies are distinct 
