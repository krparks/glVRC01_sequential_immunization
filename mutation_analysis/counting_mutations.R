#load packages necessary for analysis
library(Biostrings)
library(seqinr)
library(tidyverse)
library(data.table)

###change this every time
##need to write what you want the name of your file name to be split into. I separate these
##items in my file names by spaces
col_names_for_compiled_sequences<-c("immunogen", "adjuvant", "allele", "chain")

#read in the amino acid alignments
filenames<-list.files(pattern = ".fasta", full.names=TRUE)
seq_list<-lapply(filenames, read.alignment, format = 'fasta')

#this turns your list of alignments into a list of matrices
align_mat<-lapply(seq_list, as.matrix)
names(align_mat)<-c(filenames)

#This function makes calculates the number of mutations in all of the sequences in regards 
#to the reference sequence used. This is either the KI gene or whatever the mouse 
#light chain derivative is.
#the first function makes all of the sequences one of the matrixes into a list, excluding 
#the reference sequence. The second part of the function compares the 'query' sequences
#to the reference sequence. You should always make the reference sequence the first sequence
#in your alignment.
mismatch_matrix<-function(alignment_matrix) { 
  query<-as.list(seq(2, nrow(alignment_matrix)))
  sapply(query, function(seq_list) {
    sum(alignment_matrix[1,] != 
          alignment_matrix[seq_list, 1:ncol(alignment_matrix)])
  })}


#this applies your function to your list of matrices
mut_summ<-lapply(align_mat, FUN = mismatch_matrix)

#this function gets the names of your sequences from each of your matrices 
#this excludes the reference sequence
names<-function(alignment_matrix){
  names<-row.names(alignment_matrix)
  names<-names[-1]
}

#apply the function to your matrices
names_for_dataframe<-lapply(align_mat, names)

#combine the names and the number of mutations and create a dataframe with them
mut_dataframe<-map2(names_for_dataframe, mut_summ, data.frame)

#save the mutation summaries as csv files
mapply(write.csv, mut_dataframe, file=paste0(filenames, '.csv'), row.names = FALSE)

#read all the files in
csv_files<-list.files(pattern = ".csv")
csv_list<-lapply(csv_files, read.csv)
#name the files
names(csv_list)<-str_replace(csv_files, pattern = ".fasta.csv", replacement = "")

#give the dataframes column names
colnames<-c("sequence", "mutations")
csv_list<-lapply(csv_list, setNames, colnames)

filenames2 <- sub(".fasta.csv", "", csv_files)

#this puts the fasta name into one column called in the data frames
csv_list2<-Map(cbind, csv_list, chain = filenames2)

#this separates the chain and group into separate columns
separate_columns<-function (list) {separate(list, chain, 
                      into = col_names_for_compiled_sequences, sep = " ")}

csv_list3<-lapply(csv_list2, FUN = separate_columns)

#bind all the dataframes in the list together to form one big dataframe
csv_list4<-bind_rows(csv_list3)

#create a new column that says if the sequence from prime or boost
csv_list5<-csv_list4 %>% mutate(immunization = csv_list4$immunogen)
csv_list5$immunization<-ifelse(csv_list5$immunization == "426c_Core_Fer", "prime", csv_list5$immunization)
csv_list5$immunization<-ifelse(csv_list5$immunization == "426c_Core_C4b", "prime", csv_list5$immunization)
csv_list5$immunization<-ifelse(csv_list5$immunization == "426c_Core_C4b_HxB2_Core_C4b", "boost", csv_list5$immunization)

#write a csv file with the results
write.csv(csv_list5, "mutation_summary.csv")

csv_list5 %>% group_by(immunization, chain, adjuvant, immunogen) %>% count

#summary table
mutation_summary<- csv_list5 %>% group_by(chain, immunization) %>% summarise(
  N    = length(mutations),
  mean = mean(mutations),
  sd   = sd(mutations),
  se   = sd / sqrt(N)
)



#make colors for the graphs
fillcolors <- c("red", "blue")

#label everything correctly
x_labels <- c("prime","prime + boost")
csv_list5$chain <- factor(csv_list5$chain,
                          labels = c("VH", "VL"))
mutation_summary$chain <- factor(mutation_summary$chain,
                          labels = c("VH", "VL"))
csv_list5$immunization<-factor(csv_list5$immunization, levels=c("prime","boost"))
mutation_summary$immunization<-factor(mutation_summary$immunization, levels=c("prime","boost"))

#graph the mutation summary
mutation_graph<-ggplot(data = mutation_summary)+
  geom_errorbar(aes(x = immunization, ymin = mean -sd, ymax = mean + sd ), width = .4)+
  geom_point(data = csv_list5, aes(y = mutations, x = immunization, fill = adjuvant), 
  position=position_jitter(h=0.1, w=0.1),
             shape = 21, alpha = 0.5, size = 6)+ 
  scale_fill_manual(values=fillcolors)+
  facet_wrap(~chain)+
  theme_bw()+
  ylab("# of Amino Acid Mutations")+
  ylim(-1,14)+
  #scale_y_continuous(breaks = c(0,2,4,6,8,12))+
  theme(axis.title.x= element_blank(), 
        axis.title.y=element_text(size =35, face = "bold"), axis.text.x = element_text( size = 30, face = "bold"), 
        axis.text.y=element_text(size =30, face = "bold"), legend.text = element_text(size = 30, face = "bold"), 
        legend.title = element_text(size = 35, face = "bold"),
        strip.text = element_text(size = 30, face = "bold"))+
  scale_x_discrete(labels= x_labels)
  
#save the graph as a pdf
ggsave("mutation_graph.tiff", device = "tiff", width = 14, height = 10)


#graph the mutation summary - do not color by adjuvant
mutation_graph_do_not_color_by_adjuvant<-ggplot(data = mutation_summary)+
  geom_errorbar(aes(x = immunization, ymin = mean -sd, ymax = mean + sd ), width = .4)+
  geom_point(data = csv_list5, aes(y = mutations, x = immunization, fill = immunization), 
             position=position_jitter(h=0.1, w=0.1),
             shape = 21, alpha = 0.5, size = 4)+ 
  scale_fill_manual(values=fillcolors)+
  facet_wrap(~chain)+
  theme_bw()+
  xlab("Immunization Group")+
  ylab("# of Amino Acid Mutations")+
  ylim(-1,14)+
  #scale_y_continuous(breaks = c(0,2,4,6,8,12))+
  theme(axis.title.x= element_text(size=20), 
        axis.title.y=element_text(size =20), axis.text.x = element_text( size = 15), 
        axis.text.y=element_text(size =15), legend.text = element_text(size = 15), 
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 20))+
  scale_x_discrete(labels= x_labels)

ggsave("mutation_graph_no_adjuvant_color.tiff",mutation_graph_do_not_color_by_adjuvant, device = "tiff", width = 10, height = 6)


##Heavy chain comparison using wilcoxon test (non-parametric, unpaired test)
prime_heavy<- csv_list5 %>% filter(chain == "VH", immunization == "prime")
boost_heavy<- csv_list5 %>% filter(chain == "VH", immunization == "boost")

t.test(prime_heavy$mutations, boost_heavy$mutations, paired = FALSE)

##Light chain comparison using wilcoxon test (non-parametric, unpaired test)
prime_light<- csv_list5 %>% filter(chain == "VL", immunization == "prime")
boost_light<- csv_list5 %>% filter(chain == "VL", immunization == "boost")

t.test(prime_light$mutations, boost_light$mutations, paired = FALSE)
