setwd("/home/nikolay/Documents/Medium/DeepLearningNeanderthalIntrogression/Genes")

################################ PREPARE OLIGONUCLEOTIDE COMPOSITIONS IN GENIC AND INTERGENIC REGIONS FOR HMM ##############################

library("Biostrings")
geneseq <- readDNAStringSet("hg19_gene_clean.fa")
notgeneseq <- readDNAStringSet("hg19_notgene_clean.fa")

make_emission <- function(order, geneseq, notgeneseq)
{
  oligonucs <- colnames(oligonucleotideFrequency(geneseq, order))
  oligonucs_gene_freqs <- vector(); oligonucs_notgene_freqs <- vector()
  for(i in 1:length(oligonucs))
  {
    oligonucs_gene_freqs<-append(oligonucs_gene_freqs,
                                 sum(oligonucleotideFrequency(geneseq,order)[,oligonucs[i]])/sum(width(geneseq)))
    oligonucs_notgene_freqs<-append(oligonucs_notgene_freqs,
                                    sum(oligonucleotideFrequency(notgeneseq,order)[,oligonucs[i]])/sum(width(notgeneseq)))
  }
  emission <- t(data.frame(Gene = oligonucs_gene_freqs, NotGene = oligonucs_notgene_freqs))
  colnames(emission) <- oligonucs
  return(emission)
}

########################################################## PREPARE HMM ####################################################################
run_viterbi <- function(myseq, order, verbose = FALSE)
{
  # Define emission and transition matrices
  states <- c("Gene", "NotGene")
  emission <- as.matrix(read.delim(paste0("emission_", order, "mers.txt"), header = TRUE, row.names = 1, check.names = FALSE, sep = "\t"))
  #emission <- make_emission(order, geneseq, notgeneseq)
  if(verbose == TRUE)
  {
    print(paste0("Emission matrix of Hidden Markov Model (HMM) of ",order,"-order for gene recognition from nucleotide sequence"))
    print(emission)
  }
  transition <- matrix(c(c(0.59, 0.41), c(0.4, 0.6)), 2, 2, byrow = TRUE); rownames(transition) <- states; colnames(transition) <- states
  if(verbose == TRUE)
  {
    print("Transition matrix of Hidden Markov Model (HMM) for gene recognition from nucleotide sequence:")
    print(transition)
  }

  # Compute probabilities of genic and intergenic regions
  probs <- matrix(NA, nrow = length(myseq), ncol = length(states)); colnames(probs) <- states
  probs[1, "Gene"] <- 0.5 * emission["Gene", paste(myseq[1:order], collapse = "")]
  probs[1, "NotGene"] <- 0.5 * emission["NotGene", paste(myseq[1:order], collapse = "")]
  probs[1, ] <- probs[1, ] / sum(probs[1, ])
  for (i in 2:(length(myseq)-order+1))
  {
    probs[i, "Gene"] <- emission["Gene", paste(myseq[(i-1):(i+order-2)], collapse="")] * max(probs[(i-1), ] * transition[, "Gene"])
    probs[i, "NotGene"] <- emission["NotGene", paste(myseq[(i-1):(i+order-2)], collapse="")] * max(probs[(i-1), ] * transition[, "NotGene"])
    probs[i, ] <- probs[i, ] / sum(probs[i, ])
  }
  probs[is.na(probs[,"Gene"]),] <- matrix(rep(probs[is.na(probs[,"Gene"])==FALSE,][dim(probs[is.na(probs[,"Gene"])==FALSE,])[1],], 
                                              sum(!complete.cases(probs))), byrow=TRUE, nrow = sum(!complete.cases(probs)))
  if(verbose == TRUE)
  {
    print("Probabilities of Hidden Markov Model (HMM) for gene recognition from nucleotide sequence:")
    print(head(probs))
  }

  # Infer states from computed probabilities by tracing them back and applying correction via transition matrix
  inferred_state <- rep(NA, length(myseq))
  inferred_state[length(myseq)] <- names(which(probs[length(myseq),] == max(probs[length(myseq),])))
  for(i in (length(myseq)-1):1)
  {
    inferred_state[i] <- names(which(probs[i, ] * transition[inferred_state[i+1], ] == max(probs[i, ] * transition[inferred_state[i+1], ])))
  }
  inferred_state_num <- ifelse(inferred_state == "Gene", 1, 0)
  if(verbose == TRUE)
  {
    print("Inferred states computed by Hidden Markov Model (HMM): ")
    print(inferred_state_num)
  }
  return(list(inferred_states = inferred_state_num, score = probs[,"Gene"]))
}






########################################################### RUN HMM FOR SEQ1 ###############################################################
# Read test sequence with known ground truth for Gene (1) vs. NotGene (0) states
library("stringr")
myseq1 <- toupper(unlist(str_split(readLines("gene_notgene_seq1.fa"), pattern = "")))
myseq1
mytruth1 <- c(rep(1, 60), rep(0, 60), rep(1, 60), rep(0, 60))
mytruth1

order <- 1
inferred_states <- run_viterbi(myseq1, order)$inferred_states
par(mfrow = c(2, 2))
plot(inferred_states, cex = 0.5, type = "l", xlab = "Position on sequence", ylab = "State", 
     main = paste0("Hidden Markov Model of k = ", order," order: 0 - Not Gene, 1 - Gene"))
polygon(c(0, 0, 60, 60), c(0, 1, 1, 0), col=adjustcolor("green", alpha.f = 0.1), border = NA)
polygon(c(60, 60, 120, 120), c(0, 1, 1, 0), col=adjustcolor("red", alpha.f = 0.1), border = NA)
polygon(c(120, 120, 180, 180), c(0, 1, 1, 0), col=adjustcolor("green", alpha.f = 0.1), border = NA)
polygon(c(180, 180, 240, 240), c(0, 1, 1, 0), col=adjustcolor("red", alpha.f = 0.1), border = NA)
lines(inferred_states, lwd = 3)

table(inferred_states, mytruth1)
sum(diag(table(inferred_states, mytruth1))) / sum(table(inferred_states, mytruth1))


library("ROCit")
par(mfrow = c(2, 2))
color_vector<-c("blue", "red", "green", "orange")
plot(1, lwd=1, xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)", 
     type='n', xlim = c(0, 1), ylim = c(0, 1), main = "Sequence 1")
lines(c(0,1),c(0,1),lwd=2,lty=2,col="darkgreen")
for(order in 1:4)
{
  inferred_states <- run_viterbi(myseq1, order)$score
  myroc<-rocit(inferred_states, mytruth1, method="binormal")
  lines(myroc$TPR~myroc$FPR, col = color_vector[order])
}


########################################################### RUN HMM FOR SEQ2 ###############################################################
myseq2 <- toupper(unlist(str_split(readLines("gene_notgene_seq2.fa"), pattern = "")))
myseq2
mytruth2 <- c(rep(0, 120), rep(1, 120))
mytruth2

order <- 1
inferred_states <- run_viterbi(myseq2, order)$inferred_states
plot(inferred_states, cex = 0.5, type = "l", xlab = "Position on sequence", ylab = "State", 
     main = paste0("Hidden Markov Model of k = ", order," order: 0 - Not Gene, 1 - Gene"))
polygon(c(0, 0, 120, 120), c(0, 1, 1, 0), col=adjustcolor("red", alpha.f = 0.1), border = NA)
polygon(c(120, 120, 240, 240), c(0, 1, 1, 0), col=adjustcolor("green", alpha.f = 0.1), border = NA)
lines(inferred_states, lwd = 3)

table(inferred_states, mytruth2)
sum(diag(table(inferred_states, mytruth2))) / sum(table(inferred_states, mytruth2))

library("ROCit")
color_vector<-c("blue", "red", "green", "orange")
plot(1, lwd=1, xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)", 
     type='n', xlim = c(0, 1), ylim = c(0, 1), main = "Sequence 2")
lines(c(0,1),c(0,1),lwd=2,lty=2,col="darkgreen")
for(order in 1:4)
{
  inferred_states <- run_viterbi(myseq2, order)$score
  myroc<-rocit(inferred_states, mytruth2, method="binormal")
  lines(myroc$TPR~myroc$FPR, col = color_vector[order])
}




########################################################### RUN HMM FOR SEQ3 ###############################################################
myseq3 <- toupper(unlist(str_split(readLines("gene_notgene_seq3.fa"), pattern = "")))
myseq3
mytruth3 <- c(rep(0, 300), rep(1, 300), rep(0, 600))
mytruth3

order <- 1
inferred_states <- run_viterbi(myseq3, order)$inferred_states
plot(inferred_states, cex = 0.5, type = "l", xlab = "Position on sequence", ylab = "State", 
     main = paste0("Hidden Markov Model of k = ", order," order: 0 - Not Gene, 1 - Gene"))
polygon(c(0, 0, 300, 300), c(0, 1, 1, 0), col=adjustcolor("red", alpha.f = 0.1), border = NA)
polygon(c(300, 300, 600, 600), c(0, 1, 1, 0), col=adjustcolor("green", alpha.f = 0.1), border = NA)
polygon(c(600, 600, 1200, 1200), c(0, 1, 1, 0), col=adjustcolor("red", alpha.f = 0.1), border = NA)
lines(inferred_states, lwd = 3)

table(inferred_states, mytruth3)
sum(diag(table(inferred_states, mytruth3))) / sum(table(inferred_states, mytruth3))


library("ROCit")
color_vector<-c("blue", "red", "green", "orange")
plot(1, lwd=1, xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)", 
     type='n', xlim = c(0, 1), ylim = c(0, 1), main = "Sequence 3")
lines(c(0,1),c(0,1),lwd=2,lty=2,col="darkgreen")
for(order in 1:4)
{
  inferred_states <- run_viterbi(myseq3, order)$score
  myroc<-rocit(inferred_states, mytruth3, method="binormal")
  lines(myroc$TPR~myroc$FPR, col = color_vector[order])
}



########################################################### RUN HMM FOR SEQ4 ###############################################################
myseq4 <- toupper(unlist(str_split(readLines("gene_notgene_seq4.fa"), pattern = "")))
myseq4
mytruth4 <- c(rep(1, 120), rep(0, 60), rep(1, 120))
mytruth4

order <- 1
inferred_states <- run_viterbi(myseq4, order)$inferred_states
plot(inferred_states, cex = 0.5, type = "l", xlab = "Position on sequence", ylab = "State", 
     main = paste0("Hidden Markov Model of k = ", order," order: 0 - Not Gene, 1 - Gene"))
polygon(c(0, 0, 120, 120), c(0, 1, 1, 0), col=adjustcolor("green", alpha.f = 0.1), border = NA)
polygon(c(120, 120, 180, 180), c(0, 1, 1, 0), col=adjustcolor("red", alpha.f = 0.1), border = NA)
polygon(c(180, 180, 300, 300), c(0, 1, 1, 0), col=adjustcolor("green", alpha.f = 0.1), border = NA)
lines(inferred_states, lwd = 3)

table(inferred_states, mytruth4)
sum(diag(table(inferred_states, mytruth4))) / sum(table(inferred_states, mytruth4))

library("ROCit")
color_vector<-c("blue", "red", "green", "orange")
plot(1, lwd=1, xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)", 
     type='n', xlim = c(0, 1), ylim = c(0, 1), main = "Sequence 4")
lines(c(0,1),c(0,1),lwd=2,lty=2,col="darkgreen")
for(order in 1:4)
{
  inferred_states <- run_viterbi(myseq4, order)$score
  myroc<-rocit(inferred_states, mytruth4, method="binormal")
  lines(myroc$TPR~myroc$FPR, col = color_vector[order])
}
