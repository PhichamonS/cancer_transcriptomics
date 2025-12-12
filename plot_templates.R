####################################################################################
################## Function for complexity plot in R ###############################
####################################################################################


# You can use the function in this script by:
# (1) Copying the function to your R script
# OR
# (2) Using source("path_where_this_file_is/plot_templates.R", local = knitr::knit_global()) in your R script


# This plot shows the cumulative reads of the genes ordered by their expression, for each sample
# To use this function, please be aware that:
#     - data: should be a numeric data.frame() object, with genes as rows and samples as columns
#     - samples_to_plot: Which samples you desire to plot (names or indexes). If you want to use all samples, we suggest samples_to_plot=colnames(yourdata)
#     - reads_in_percentage: boolean, if you wish to plot the cumulative reads by absolute value (FALSE, default) or by percentage of the total number of reads per sample (TRUE)
# Output: list with the following entries:
#     - plot: ggplot object, that you can use to plot and change parameters
#     - data_to_plot: the data used to generate the plot


# Examples:
# plotcomplexity_1 <- complexityplot_byexprgenes(data=yourdata,samples_to_plot=c("TCGA-DM-A288-01A","TCGA-CM-6164-01A"), reads_in_percentage=FALSE)
# plotcomplexity_1$plot #gives you the ggplot object 
# plotcomplexity_1$data_to_plot # gives you the data.frame() object used to generate the plot

# plotcomplexity_2 <- complexityplot_byexprgenes(data=yourdata,samples_to_plot=c(1:12), reads_in_percentage=TRUE) #this option uses the first 12 samples only

library(ggplot2)
library(tximport)
library(data.table)
library(biomaRt)
library(dplyr)

complexityplot_byexprgenes <- function(data, samples_to_plot, reads_in_percentage=TRUE){
  
  final <- c()
  #create new dataset that contains cumulative reads
  for (sample in samples_to_plot){
    temporary<- data[order(data[,sample],decreasing = TRUE),]
    temporary <- cumsum(temporary[,sample])
    final <- cbind(final, temporary)
  }
  final <- as.data.frame(final)
  colnames(final) <- samples_to_plot
  final <- as.data.table(final)
  
  
  #this is the plot
  finalmelt <- melt(final)
  finalmelt <- cbind(finalmelt,rep(x=c(1:nrow(final)),length(unique(finalmelt$variable))))
  colnames(finalmelt) <- c("Sample","CumulativeReads","GeneIndex")
  
  if(!reads_in_percentage){
    (plt1 <- ggplot(data = finalmelt, aes(x=GeneIndex, y=log10(CumulativeReads), group=Sample, color=Sample)) + geom_line() +
       theme_bw() + theme(legend.position="none")   + xlab("Gene Index") + ylab("Number of Reads (log10)") + 
       theme(axis.line = element_line(colour = "black"),
             panel.grid.minor = element_blank(),
             panel.border = element_blank(),
             panel.background = element_blank()))
    data_to_plot <- finalmelt
  } else {
    orderfinal <- subset(finalmelt, GeneIndex == max(finalmelt$GeneIndex))
    #as percentage
    trial <- merge(orderfinal, finalmelt, by = "Sample", all=TRUE)
    trial <- trial[,-3]
    trial$percentage <- (trial$CumulativeReads.y / trial$CumulativeReads.x)
    trial$percentage <- trial$percentage * 100
    trial <- trial[,-c(2,3)]
    colnames(trial) <- c("Sample","GeneIndex","PercentageCumulativeReads")
    (plt1 <- ggplot(data = trial, aes(x=GeneIndex, y=PercentageCumulativeReads, group=Sample, color=Sample)) + geom_line() +
        theme_bw() +theme(legend.position="none") + xlab("Gene Index") + ylab("Number of Reads (%)")+
        theme(axis.line = element_line(colour = "black"),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank()))
    
    data_to_plot <- trial
  }
  
  return(list(plot=plt1,
              data_to_plot=data_to_plot))
}

# Load TCGA data
rawcounts <- as.data.frame(fread('C:/Users/Utilizador/Downloads/TCGA_KIDNEY_Gene_ReadCounts.txt'))
row.names(rawcounts) <- rawcounts[,1]
rawcounts <- rawcounts[,-1]
# Remove missing values
rawcounts <- rawcounts[!rowSums(is.na(rawcounts)) == ncol(rawcounts),]
rawcounts <- rawcounts[, colSums(is.na(rawcounts)) < ncol(rawcounts)]


plot = complexityplot_byexprgenes(data, colnames(data))
