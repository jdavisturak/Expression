## Analyze NfKB of a certain region

library('Biostrings')

### FUCNTION TO SCORE THE ENTIRE STRING, THEN (optionally) RE-ORDER BY SCORE .
scorePWM = function(mySeq,pwm,starts = 1:(length(mySeq)-ncol(pwm)+1),REORDER=TRUE){
  scores  = sapply(starts,function(start) PWMscoreStartingAt(log(pwm),mySeq,start))
  strings = sapply(starts,function(start) as.character(subseq(mySeq,start=start,width=ncol(pwm))))

  dat = data.frame(location=starts,score=scores, sequence=strings)

  if(REORDER)
    return(dat[order(dat$score,decreasing=TRUE),])

  dat
}

# Read in the PWM matrix
as.matrix(read.table('NfKB.pwm.txt',row.names=1))->pwm

# Read in DNA sequence
mySeq = read.DNAStringSet('IL10.fa')[[1]] # Reads in a fasta file: [[1]] keeps ONLY the first record
mySeq = read.DNAStringSet('IL10_5000.fa')[[1]] # Reads in a fasta file: [[1]] keeps ONLY the first record



### CODE TO SCORE A CERTAIN SECTION (e.g. starting at base pair 42):
#scorePWM(mySeq, pwm,42:50)

### CODE TO SCORE A CERTAIN SECTION (e.g. starting at base pair 42):
scorePWM(mySeq, pwm,REORDER=TRUE) -> topScores
head(topScores,20)


