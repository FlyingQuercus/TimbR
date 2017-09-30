# Functions
This directory contains R functions required for running the scripts.  Look at the comments on each file for details on how to use the functions

## simhit
This function creates simulated impulse excitation hits that match theory (sine waves with exponential damping).
The function is quite slow and could be made more efficient by only calculating the non-zero amplitudes.

## findhits
This function exampines a sound recording and makes a list of times when the hits start and end.
This function cannot tell if there are additional subsequent hits below the trigger threshold.

## expwin
This function applies an exponential window (both in and out) to a hit recording, to reduce problems caused by background noise, or the recording ending before the vibration has decayed sufficiently. Any subsequent damping calculations need to be corrected accordingly

## logdecrement
This function calculates damping (decrement, decay constant) on the time domain (only suitable when there is one mode of vibration present. 

## peaksQ
This function calculates peak frequencies and damping (bandwidth Q) on an FFT frequency spectrum (in dB). It can cope with multiple modes present, as long as the peak frequencies are distinct 
