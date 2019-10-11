library(dplyr)
library(readr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

input_file = args[1]
output = args[2]


data  = read_tsv(input_file, col_names= FALSE)


#Compute Media
data %>% filter(X4 != 0) %>% summarize( Median_depth = median(X4) ) -> median_estimation

#Compute distribution
data %>% filter(X4 != 0) %>% ggplot() + geom_density(aes(X4)) + geom_vline(xintercept = median_estimation$Median_depth, linetype="dotted", 
                color = "red", size=1) + theme_bw() -> distribution


#Compute genome length
data %>% mutate(Range= X3 - X2) %>% summarise(sum(Range)) -> total
total = total[[1]]
#Compute Genome coverage
data %>% mutate(Range= X3 - X2) %>% mutate(present = ifelse(X4 != 0, 1, 0)) %>%
 group_by(present) %>% summarise(coverage=sum(Range)/total) %>% filter(present=="1") -> Covered
Covered = Covered$coverage

#Two column output, Median Depth and Genome coverage
median_estimation %>% mutate(Coverage = Covered) -> result

 
 
write.table(result, file=output, quote=FALSE, sep=',', row.names = FALSE)
ggsave(filename=paste(output,"_distribution.pdf",sep = ""), plot=distribution)

