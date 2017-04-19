rm(list=ls())
library("Biostrings")
library("annotate")
setwd("~/Projects/testProject/")

(sequences.seq <- readBStringSet("scaffold_prva_AT.fa","fasta"))

(sequences.ordered.seq <- sequences.seq[order(width(sequences.seq), decreasing=TRUE)] )

#(sub.sequence.seq <- sequences.seq[width(sequences.seq) > 15000])
#(sub.sequence.seq <- sub.sequence.seq[order(width(sub.sequence.seq), decreasing=TRUE)])

#writeXStringSet(sub.sequence.seq[order((width(sub.sequence.seq)), decreasing=TRUE)] , "sub.sequence.seq", format="fasta")

#sum(width(sequences.seq))
#sum(width(sub.sequence.seq))

#subseq(x, start=NA, end=NA, width=NA)
#subseq(sub.sequence.seq[order((width(sub.sequence.seq)), decreasing=TRUE)][1],1,38)

#(sub.sequence.seq[order((width(sub.sequence.seq)), decreasing=TRUE)])[1]


#TRY TO BLASTN A SEQUENCE
#(X <- sub.sequence.seq[order((width(sub.sequence.seq)), decreasing=TRUE)][1])
#(x <- toString(sub.sequence.seq[order((width(sub.sequence.seq)), decreasing=TRUE)][7]))

#system2('/home/stefano/Programs/ncbi-blast-2.6.0+/bin/blastn'
#                   ,c('-db',"nt"
#                   ,'-outfmt',"0"
#                   ,'-perc_identity',".90"
#                   ,'-remote'
#                    )
#                    ,input=x
#                    ,stdout=TRUE )


##compute the GC content in a sliding window (as a fraction) for a sequence
letterFrequency(sequences.ordered.seq[1], letters="NACGT", OR=0)
alphabetFrequency(DNAString(sequences.ordered.seq[[1]]))
hasOnlyBaseLetters(DNAString(sequences.ordered.seq[[1]]))
uniqueLetters(DNAString(sequences.ordered.seq[[1]]))
sum(alphabetFrequency(sequences.ordered.seq[[1]])) == length(sequences.ordered.seq[[1]])

maskedSeq <- maskMotif(sequences.ordered.seq[[1]], "N")


window <- 15000
AT <- rowSums(letterFrequencyInSlidingView(DNAString(maskedSeq), window, c("A", "T")))/window
plot(AT, type = 'l')
lines(lowess(x = 1:length(AT), y= AT, f = 0.10), col = 12, lwd = 2)
