
# load required libraries
library(seewave)
library(tuneR)

# first select any wav file in the directory to scan through 
wavname = file.choose()

# make a list of all the wav files in this directory
wavList = list.files(path = dirname(wavname),pattern="*.wav",full.names = TRUE)

allresults = -1

print(wavList)
cat("Found",length(wavList),"wav files \n")

for(i in 1:length(wavList)){
  
  cat("Working on file",i,"of",length(wavList),"  ",basename(wavList[i])," \n")  
  wav.recorded=readWave(filename=wavList[i])
  IDuse = basename(wavList[i])
 
  result = quickpeaks(wav.recorded,
                       
                      plotpeaks = FALSE,
                      plothits = FALSE,
                      mindB=-5,
                      peakscount = 5,
                      ID=IDuse
                      )

  if(i == 1) {
  allresults = result  
  } else {
  allresults = rbind(allresults, result)  
  }
  
}
  
write.csv(allresults, file = "allresults.csv", row.names=FALSE)
  

