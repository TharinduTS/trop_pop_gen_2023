## trop_pop_gen_2023

# FST

I started with putting all bam files in a single folder

Then created two text files as pop1 and pop2 containing list of samples for each of the populations ( Males and females in this case )


pop1
```txt
F_Ghana_WZ_BJE4687_combined__sorted.bam_rg_rh.bam
F_IvoryCoast_xen228_combined__sorted.bam_rg_rh.bam
F_Nigeria_EUA0331_combined__sorted.bam_rg_rh.bam
F_Nigeria_EUA0333_combined__sorted.bam_rg_rh.bam
F_SierraLeone_AMNH17272_combined__sorted.bam_rg_rh.bam
F_SierraLeone_AMNH17274_combined__sorted.bam_rg_rh.bam
all_ROM19161_sorted.bam
XT10_WZ_no_adapt._sorted.bam_rg_rh.bam
XT11_WW_trim_no_adapt_scafconcat_sorted.bam_rg_rh.bam
```


pop2
```txt
M_Ghana_WY_BJE4362_combined__sorted.bam_rg_rh.bam
M_Ghana_ZY_BJE4360_combined__sorted.bam_rg_rh.bam
M_Nigeria_EUA0334_combined__sorted.bam_rg_rh.bam
M_Nigeria_EUA0335_combined__sorted.bam_rg_rh.bam
M_SierraLeone_AMNH17271_combined__sorted.bam_rg_rh.bam
M_SierraLeone_AMNH17273_combined__sorted.bam_rg_rh.bam
XT1_ZY_no_adapt._sorted.bam_rg_rh.bam
XT7_WY_no_adapt__sorted.bam_rg_rh.bam
```

Then used the following script to cal FST changing essential fields

**removed "-r Chr7:1-50000000" to cal FST for the wholw gwnome **

cal_fst.sh

```bash
#!/bin/sh
#SBATCH --job-name=fst
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --mem=512gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

# Load modules
module load nixpkgs/16.09  gcc/7.3.0
module load angsd
module load gsl/2.5
module load htslib

# Define directory containing BAM files
bamdir=./

# Define output directory for Fst results
outdir=./FST_outs_chr7_only

# Define reference genome
refgenome=./reference_genome/XENTR_10.0_genome_scafconcat_goodnamez.fasta

# Loop over all BAM files in the directory
for bamfile in ${bamdir}/*.bam; do
    # Extract the sample name from the file name
    sample=$(basename ${bamfile} .bam)


#this is with 2pops
#first calculate per pop saf for each populatoin
angsd -b pop1 -doMaf 0 -doCounts 1 -setMinDepthInd 5 -setMaxDepth 100 -minMapQ 20 -anc ${refgenome} -out ${outdir}/${sample}_pop1 -dosaf 1 -gl 1 -r Chr7:1-50000000
angsd -b pop2 -doMaf 0 -doCounts 1 -setMinDepthInd 5 -setMaxDepth 100 -minMapQ 20 -anc ${refgenome} -out ${outdir}/${sample}_pop2 -dosaf 1 -gl 1 -r Chr7:1-50000000

#calculate the 2dsfs prior
realSFS ${outdir}/${sample}_pop1.saf.idx ${outdir}/${sample}_pop2.saf.idx >pop1.pop2.ml

#prepare the fst for easy window analysis etc
realSFS fst index ${outdir}/${sample}_pop1.saf.idx ${outdir}/${sample}_pop2.saf.idx -sfs pop1.pop2.ml -fstout here

#get the global estimate
realSFS fst stats here.fst.idx

#below is not tested that much, but seems to work
realSFS fst stats2 here.fst.idx -win 5000 -step 5000 >slidingwindow

done
```
To plot all the populations together,

First I downloaded and renamed slidinwindow files as 'ghana_inds,niger_inds' etc so the following R script can prepare columns when connecting data frames as needed.

then I did put them all in a seperate folder and ran the following R script

***** IF you re run this script, remove any other files other than R script and pop files so the pattern matches the file selection****

```R
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))

library(ggplot2)

# Delete previous plot before starting if it is in the same folder*****************


#remove scientific notation
options(scipen=999)

pop_name<-""

file_list<-grep(list.files(path="./"), pattern='*.R', invert=TRUE, value=TRUE)

# combine files adding pop and sex from file name
df <- do.call(rbind, lapply(file_list, function(x) cbind(read.table(x), file_name=x)))

# Split name column into firstname and last name
require(stringr)
df[c('pop', 'sex')] <- str_split_fixed(df$file_name, '_', 2)

# Subset Rows by column value
df<-df[df$region == 'Chr7',]

#change names

names(df)<-c("Chr","mid_pos","Nsites","Fst","File_name","pop","sex")

#list populations

pop_list<-unique(df$pop)

p<-ggplot(df,aes(x=mid_pos/1000000,y=Fst,col=Fst))+
  geom_point(alpha=)+
  geom_smooth()+
  #labs(title = expression(paste(italic(F[ST]),"  ","Ghana only","    -      Chr7")))+
  #ylim(-0.05, 0.25)+
  theme_bw()+
  xlab("Chrom_position")+
  ylab(expression(paste(italic(F[ST]))))+
  theme_classic()
  
#chgange order here
p+facet_wrap(~factor(pop,levels=c("all","notad","sierra","ghana","tads","niger")),ncol=1,labeller = as_labeller(c(all='All Individuals', notad='No tadpoles', sierra='Sierra Leone', ghana='Ghana', tads = 'Tadpoles',niger = 'Nigeria')))


ggsave(paste("20M_",pop_name,"_no_point_plot.pdf",sep = ""))
```
# Admixture

Made a sub folder inside previous folder with bam files and started Admixture there

first created the text file bam file list

bam.filelist
```txt
../F_Ghana_WZ_BJE4687_combined__sorted.bam_rg_rh.bam
../F_IvoryCoast_xen228_combined__sorted.bam_rg_rh.bam
../F_Nigeria_EUA0331_combined__sorted.bam_rg_rh.bam
../F_Nigeria_EUA0333_combined__sorted.bam_rg_rh.bam
../F_SierraLeone_AMNH17272_combined__sorted.bam_rg_rh.bam
../F_SierraLeone_AMNH17274_combined__sorted.bam_rg_rh.bam
../all_ROM19161_sorted.bam
../XT10_WZ_no_adapt._sorted.bam_rg_rh.bam
../XT11_WW_trim_no_adapt_scafconcat_sorted.bam_rg_rh.bam
../M_Ghana_WY_BJE4362_combined__sorted.bam_rg_rh.bam
../M_Ghana_ZY_BJE4360_combined__sorted.bam_rg_rh.bam
../M_Nigeria_EUA0334_combined__sorted.bam_rg_rh.bam
../M_Nigeria_EUA0335_combined__sorted.bam_rg_rh.bam
../M_SierraLeone_AMNH17271_combined__sorted.bam_rg_rh.bam
../M_SierraLeone_AMNH17273_combined__sorted.bam_rg_rh.bam
../XT1_ZY_no_adapt._sorted.bam_rg_rh.bam
../XT7_WY_no_adapt__sorted.bam_rg_rh.bam
```

then prepared the files for Admixture analysis with

prep_admix.sh

```bash
#!/bin/sh
#SBATCH --job-name=fst
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=512gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

# Load modules
module load nixpkgs/16.09  gcc/7.3.0
module load angsd
module load gsl/2.5
module load htslib

angsd -GL 1 -out genolike -nThreads 10 -doGlf 2 -doMajorMinor 1 -doMaf 0 -doCounts 1 -setMinDepthInd 5 -setMaxDepth 100 -minMapQ 20 -bam bam.filelist
```
cal Admix

```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --mem=256gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben
#SBATCH --array=2-5

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

# Load modules
module load nixpkgs/16.09  gcc/7.3.0
module load angsd
module load gsl/2.5
module load htslib

NGSadmix -likes genolike.beagle.gz -K ${SLURM_ARRAY_TASK_ID} -P 4 -o myoutfiles${SLURM_ARRAY_TASK_ID} -minMaf 0.05
```
then downloaded the out files and ran the following R script
```R
#!/usr/bin/Rscript

# Usage: plotADMIXTURE.r -p <prefix> -i <info file, 2-column file with ind name and population/species name> 
#                        -k <max K value> -l <comma-separated list of populations/species in the order to be plotted>
# This R script makes barplots for K=2 and all other K values until max K (specified with -k). It labels the individuals 
# and splits them into populations or species according to the individual and population/species names in the 2-column file specified with -i.
# The order of populations/species follows the list of populations/species given with -l.
# Usage example: plotADMIXTURE.r -p fileXY -i file.ind.pop.txt -k 4 -pop pop1,pop2,pop3
# In this example, the script would use the files fileXY.2.Q, fileXY.3.Q, fileXY.4.Q to make barplots for the three populations.
# file.ind.pop.txt should contain one line for each individual in the same order as in the admixture files e.g.
# ind1 pop1
# ind2 pop1
# ind3 pop2
# ind4 pop3

# Author: Joana Meier, September 2019

# I am customizing this script. Making it easy to use with R studio-TS
# ***keeping this script in current directory for clarity. But setting the working directory to the directory with all output files***


#set working directory to all_outputs inside the current path
setwd(paste(dirname(rstudioapi::getSourceEditorContext()$path),"/outfiles",sep=""))

#set default values here so it can run on R studio without parsing options
#These will come into action if you do not define these options like you do in bash

#change prefix here
p_input<-"myoutfiles"
#change sample list file here
i_input<-paste(dirname(rstudioapi::getSourceEditorContext()$path),"/pop_info.txt",sep="")
# change maximum k value to plot here
k_input<-5
# change minimum k value to plot here
m_input<-2
#add a list of populations seperated by commas here. This should be exactly similiar to the populations in your sample list file. plots will be created according to this population order
# you will have to change this every time you edit sample list
l_input<-"Liberia,Sierra_Leone,Ivory_coast,Ghana,Lab_tads,Nigeria"
#add the location and file name for the plots here. I am setting this to the directory I am creating in the next line
o_input<-paste(dirname(rstudioapi::getSourceEditorContext()$path),"/plot_outs/plot",sep="")

#create a directory for plot output if it doesn't already exist in the directory with the script
dir.create(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "plot_outs"))

# Read in the arguments
library("optparse")
option_list = list(
  make_option(c("-p", "--prefix"), type="character", default=p_input, 
              help="prefix name (with path if not in the current directory)", metavar="character"),
  make_option(c("-i", "--infofile"), type="character", default=i_input, 
              help="info text file containing for each individual the population/species information", metavar="character"),
  make_option(c("-k", "--maxK"), type="integer", default=k_input, 
              help="maximum K value", metavar="integer"),
  make_option(c("-m", "--minK"), type="integer", default=m_input, 
              help="minimum K value", metavar="integer"),
  make_option(c("-l", "--populations"), type="character", default=l_input, 
              help="comma-separated list of populations/species in the order to be plotted", metavar="character"),
  make_option(c("-o", "--outPrefix"), type="character", default=o_input, 
              help="output prefix (default: name provided with prefix)", metavar="character")
) 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check that all required arguments are provided
if (is.null(opt$prefix)){
  print_help(opt_parser)
  stop("Please provide the prefix", call.=FALSE)
}else if (is.null(opt$infofile)){
  print_help(opt_parser)
  stop("Please provide the info file", call.=FALSE)
}else if (is.null(opt$maxK)){
  print_help(opt_parser)
  stop("Please provide the maximum K value to plot", call.=FALSE)
}else if (is.null(opt$populations)){
  print_help(opt_parser)
  stop("Please provide a comma-separated list of populations/species", call.=FALSE)
}

# If no output prefix is given, use the input prefix
if(opt$outPrefix=="default") opt$outPrefix=opt$prefix

# Assign the first argument to prefix
prefix=opt$prefix

# Get individual names in the correct order
labels<-read.table(opt$infofile)

# Name the columns
names(labels)<-c("ind","pop")

# Add a column with population indices to order the barplots
# Use the order of populations provided as the fourth argument (list separated by commas)
labels$n<-factor(labels$pop,levels=unlist(strsplit(opt$populations,",")))
levels(labels$n)<-c(1:length(levels(labels$n)))
labels$n<-as.integer(as.character(labels$n))

# read in the different admixture output files
minK=opt$minK
maxK=opt$maxK
tbl<-lapply(minK:maxK, function(x) read.table(paste0(prefix,x,".qopt")))

# Prepare spaces to separate the populations/species
rep<-as.vector(table(labels$n))
spaces<-0
for(i in 1:length(rep)){spaces=c(spaces,rep(0,rep[i]-1),0.15)}
spaces<-spaces[-length(spaces)]

#change the space between individuals
spaces<-replace(spaces, spaces==0, 0.02)

#change the space between populations- only change pop_space

pop_space<-0.10
spaces<-replace(spaces, spaces==0.15, pop_space)


# Plot the cluster assignments as a single bar for each individual for each K as a separate row
tiff(file=paste0(opt$outPrefix,".tiff"),width = 2000, height = 1300,res=200)
# change 3rd value in "oma" if you wanna add labels to gain enough space. Use mai to adjust space between different k valued plots
par(mfrow=c(maxK-1,1),mar=c(0,1,0.1,1),oma=c(2,4,4,1),mgp=c(0,0.2,0),mai=c(0.05,0,0,0),xaxs="i",cex.lab=1.2,cex.axis=0.8)

#Create a list for plots
plot_list<-list()
#assign colors to populations at once

col1<-"red"
col2<-"green"
col3<-"orange"
col4<-"skyblue"
col5<-"darkblue"
col6<-"yellow"
col7<-"pink"
col8<-"purple"
col9<-"grey"
col10<-"black"
col11<-"forestgreen"
col12<-"brown"

#create color palettes for each k value
col_palette_k2<-c(col1,col2)
col_palette_k3<-c(col1,col2,col3)
col_palette_k4<-c(col3,col1,col4,col2)
col_palette_k5<-c(col3,col2,col5,col1,col4)


col_palette_k6<-c(col1,col2,col6,col5,col4,col3)
col_palette_k7<-c(col1,col7,col5,col6,col4,col3,col2)
col_palette_k8<-c(col5,col8,col1,col7,col6,col3,col4,col2)
col_palette_k9<-c(col8,col6,col7,col2,col5,col4,col3,col1,col9)
col_palette_k10<-c(col8,col7,col4,col10,col2,col9,col6,col1,col3,col5)
col_palette_k11<-c(col6,col10,col2,col4,col1,col7,col9,col8,col3,col11,col5)
col_palette_k12<-c(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12)

# make a variable with full sample names 
full_sample_names<-labels$ind[order(labels$n)]

# get simpler names  for those complex sample names
library(plyr)
#first renaming the sample names which does not follow the specific format manually**DO NOT USE"_" HERE AS SPLIT WILL USE THIS IN NEXT LINE
simple_sample_list_changing<-mapvalues(full_sample_names, from = c("all_ROM19161_sorted.bam","XT10_WZ_no_adapt._sorted.bam","XT11_WW_trim_no_adapt_scafconcat_sorted.bam","XT1_ZY_no_adapt._sorted.bam","XT7_WY_no_adapt__sorted.bam"), 
                                       to = c("o_Liberia1","F_Lab1","F_Lab2","M_Lab3","M_Lab4"))

#convert facter list into chrs to rename
s_list_chr<-as.character(simple_sample_list_changing)


#then remove the parts after the first"_" from other samples
shortened_sample_list<-sapply(strsplit(s_list_chr,split = "_"),`[`, 2)


# following plots are written in a way you can just copy paste only changing k_val for the different number of 'k's
#paste whats inside * marks for different k values and then change k_val
#**********
# Plot k=2
# change only here
k_val<-2




bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
axis(3,at=bp,labels=shortened_sample_list,las=2,tick=F,cex=0.6)




#***********

#commwnt out this section to automate coloring. Keep commented to manually assign colours for each pop with the next section( you can do it once by selecting lines and ctrl+shift+c)
# Plot higher K values
#  if(maxK>minK)lapply(2:(maxK-1), function(x) barplot(t(as.matrix(tbl[[x]][order(labels$n),])), col=rainbow(n=x+1),xaxt="n", border=NA,ylab=paste0("K=",x+1),yaxt="n",space=spaces))
#  axis(1,at=c(which(spaces==0.15),bp[length(bp)])-diff(c(1,which(spaces==pop_space),bp[length(bp)]))/2,
#      labels=unlist(strsplit(opt$populations,",")))
# dev.off()

#**********
# Plot k
# change only here
k_val<-3




bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
#axis(1,at=c(which(spaces==0.15),bp[length(bp)])-diff(c(1,which(spaces==pop_space),bp[length(bp)]))/2,
#     labels = FALSE)




#***********

#**********
# Plot k
# change only here
k_val<-4




bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
#axis(1,at=c(which(spaces==0.15),bp[length(bp)])-diff(c(1,which(spaces==pop_space),bp[length(bp)]))/2,
#     llabels = FALSE)



#***********

#**********
# Plot k
# change only here
k_val<-5

#uncomment and change values of this and line after bp to manually change label placing
label_points<-c(0.50,3.10,5.80,8.00,11.80,15.72)

bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
#axis(1,at=c(which(spaces==pop_space),bp[length(bp)])-diff(c(1,which(spaces==pop_space),bp[length(bp)]))/2,
#using following line to use manual measures
axis(1,at=label_points,
     labels=unlist(strsplit(opt$populations,",")))



# #***********
# #**********
# # Plot k
# # change only here
# k_val<-12
# 
# 
# 
# 
# bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
# axis(1,at=c(which(spaces==0.15),bp[length(bp)])-diff(c(1,which(spaces==0.15),bp[length(bp)]))/2,
#      labels=unlist(strsplit(opt$populations,",")))

dev.off()

```

# Calculate depth

made a new directory inside the directory with bam files and ran this script

```bash
#!/bin/sh
#SBATCH --job-name=depth
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=128gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL


module load bwa
module load samtools/1.10
for i in ./../*.bam ; do samtools depth -a -r Chr7:1-20000000 $i > ./$i"_depth" ; done
```
Then I calculated windowed depth with the following script

```R
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))

library(reshape) # to rename columns
library(data.table) # to make sliding window dataframe
library(zoo) # to apply rolling function for sliding window
library(ggplot2)

file_list<-grep(list.files(path="./"), pattern='*bam_depth', value=TRUE)

for (filex in file_list) {
  

#upload data to dataframe, rename headers, make locus continuous, create subsets
depth <- read.table(filex, sep="\t", header=T)

# give new column names
names(depth)<-c('Chr','pos','depth')

#install.packages("data.table")
library(data.table)
#install.packages("zoo")
library(zoo)

Xdepth<-subset(depth, select = c("Chr", "pos","depth"))

#genome coverage as sliding window
Xdepth.average<-setDT(Xdepth)[, .(
  window.start = rollapply(pos, width=5000, by=5000, FUN=min, align="left", partial=TRUE),
  window.end = rollapply(pos, width=5000, by=5000, FUN=max, align="left", partial=TRUE),
  coverage = rollapply(depth, width=5000, by=5000, FUN=mean, align="left", partial=TRUE)
), .(Chr)]

write.table(Xdepth.average, file=paste(filex,".windowed",sep = ''), quote=FALSE, sep='\t', col.names = NA)

}

```
Then plotted data with

```R
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))

library(ggplot2)

# Delete previous plot before starting if it is in the same folder*****************


#remove scientific notation
options(scipen=999)

pop_name<-""

file_list<-grep(list.files(path="./"), pattern='\\windowed$', value=TRUE)

print(file_list)

# combine files adding pop and sex from file name
df <- do.call(rbind, lapply(file_list, function(x) cbind(read.table(x), file_name=x)))

# give new column names
names(df)<-c("Chr","Start","End","Depth","file_name")

# Split name column into firstname and last name
require(stringr)
df[c('sex', 'pop','ind','info')] <- str_split_fixed(df$file_name, '_', 4)

# Subset Rows by column value
df<-df[df$Chr == 'Chr7',]

# Subset Rows by column value
df_females<-df[df$sex == 'F',]
df_males<-df[df$sex == 'M',]

#create an empty column for means/males
df_males$ind_mean<-"NA"

#add mean column males
male_means<-tapply(df_males$Depth,df_males$ind, mean)
male_mean_df<-as.data.frame(male_means)
# DataFrame and put name of function rn
male_mean_df <- tibble::rownames_to_column(male_mean_df, "ind")

#conditionally fill the column comparing values
for (i in 1:length(male_mean_df$ind)) {
  for (j in 1:length(df_males$Depth)) {
    if (male_mean_df$ind[i]==df_males$ind[j]) {
      df_males$ind_mean[j]<-male_mean_df$male_means[i]
    }
    j=j+1
  }
}

#add a column with std depth
df_males$std_depth<-as.numeric(df_males$Depth)/as.numeric(df_males$ind_mean)

#create an empty column for means/females
df_females$ind_mean<-"NA"

#add mean column females
male_means<-tapply(df_females$Depth,df_females$ind, mean)
male_mean_df<-as.data.frame(male_means)
# DataFrame and put name of function rn
male_mean_df <- tibble::rownames_to_column(male_mean_df, "ind")

#conditionally fill the column comparing values
for (i in 1:length(male_mean_df$ind)) {
  for (j in 1:length(df_females$Depth)) {
    if (male_mean_df$ind[i]==df_females$ind[j]) {
      df_females$ind_mean[j]<-male_mean_df$male_means[i]
    }
    j=j+1
  }
}

#add a column with std depth
df_females$std_depth<-as.numeric(df_females$Depth)/as.numeric(df_females$ind_mean)

#list populations

pop_list<-unique(df$pop)

print(pop_list)

p<-ggplot(df_males,aes(x=Start/1000000,y=std_depth,color=ind))+
  #geom_point(alpha=0.1)+
  geom_smooth()+
  scale_color_manual(values = c("BJE4360"="blue","BJE4362"='blue',"EUA0334"='blue',"EUA0335"='blue',"AMNH17271"='blue',"AMNH17273"='blue',"XT1"='blue',"XT7"='blue',"BJE4687"='red',"xen228"='red',"all"='red',"GermSeq"='red',"EUA0331"='red',"EUA0333"='red',"JBL052"='red',"AMNH17272"='red',"AMNH17274"='red',"XT10"='red',"XT11"='red' ))+
  geom_smooth(data = df_females,aes(x=Start/1000000,y=std_depth,color=ind))+
  #scale_color_manual(labels = c("Males", "Females"), values = c("blue", "red"))+
  labs(title = pop_name)+
  #ylim(-0.05, 0.25)+
  xlim(1,20)+
  theme_bw()+
  xlab("Chrom_position - Chr 7 (mb)")+
  ylab("Depth")+
  theme_classic()+ 
  theme(legend.position="none")

#chgange order here
p+facet_wrap(~factor(pop,levels=c("Liberia","SierraLeone","IvoryCoast","Ghana","tad","Nigeria","scaffs","mello")),ncol=1,labeller = as_labeller(c(Liberia='Liberia',SierraLeone='Sierra Leone',IvoryCoast='Ivory Coast',Ghana='Ghana',tad='Tadpoles',Nigeria='Nigeria',scaffs='Scaffold data',mello='Mellotropicalis data')))


ggsave(paste("20M_",pop_name,"_depth_plot.pdf",sep = ""))
```

# Nucleotide diversity

created a bam file list for each of the different populations.

bam.filelist
```bash
../all_ROM19161_sorted.bam
../F_Ghana_WZ_BJE4687_combined__sorted.bam_rg_rh.bam
../F_IvoryCoast_xen228_combined__sorted.bam_rg_rh.bam
../F_Nigeria_EUA0331_combined__sorted.bam_rg_rh.bam
../F_Nigeria_EUA0333_combined__sorted.bam_rg_rh.bam
../F_SierraLeone_AMNH17272_combined__sorted.bam_rg_rh.bam
../F_SierraLeone_AMNH17274_combined__sorted.bam_rg_rh.bam
../JBL052_concatscafs_sorted.bam_rg_rh.bam
../mello_GermSeq_sorted.bam_rg_rh.bam
../M_Ghana_WY_BJE4362_combined__sorted.bam_rg_rh.bam
../M_Ghana_ZY_BJE4360_combined__sorted.bam_rg_rh.bam
../M_Nigeria_EUA0334_combined__sorted.bam_rg_rh.bam
../M_Nigeria_EUA0335_combined__sorted.bam_rg_rh.bam
../M_SierraLeone_AMNH17271_combined__sorted.bam_rg_rh.bam
../M_SierraLeone_AMNH17273_combined__sorted.bam_rg_rh.bam
../XT10_WZ_no_adapt._sorted.bam_rg_rh.bam
../XT11_WW_trim_no_adapt_scafconcat_sorted.bam_rg_rh.bam
../XT1_ZY_no_adapt._sorted.bam_rg_rh.bam
../XT7_WY_no_adapt__sorted.bam_rg_rh.bam
```
then calculated tP(π) with

cal_pi_for_inds.sh
```bash
#!/bin/sh
#SBATCH --job-name=fst
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=64gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben
#SBATCH --array=1-19

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

# Load modules
module load nixpkgs/16.09  gcc/7.3.0
module load angsd
module load gsl/2.5
module load htslib

angsd -bam ind${SLURM_ARRAY_TASK_ID} -doMaf 0 -doCounts 1 -setMinDepthInd 5 -setMaxDepth 100 -minMapQ 20 -doSaf 1 -anc ../../reference_genome/XENTR_10.0_genome_scafconcat_goodnamez.fasta -GL 1 -P 24 -out out_${SLURM_ARRAY_TASK_ID}

realSFS out_${SLURM_ARRAY_TASK_ID}.saf.idx -P 24 > out_${SLURM_ARRAY_TASK_ID}.sfs
realSFS saf2theta out_${SLURM_ARRAY_TASK_ID}.saf.idx -sfs out_${SLURM_ARRAY_TASK_ID}.sfs -outname out_${SLURM_ARRAY_TASK_ID}
thetaStat do_stat out_${SLURM_ARRAY_TASK_ID}.thetas.idx -win 5000 -step 5000  -outnames theta.thetasWindow_${SLURM_ARRAY_TASK_ID}.gz
```
After downloading the files, 

calculated moving average by windows

```R
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))

library(reshape) # to rename columns
library(data.table) # to make sliding window dataframe
library(zoo) # to apply rolling function for sliding window
library(ggplot2)

#change window size here ************
win_size<-5000

#change moving win size here ***********
mov_win_size<-5000

file_list<-grep(list.files(path="./Individual_NDs",full.names = TRUE), pattern='*', value=TRUE)

print(file_list)

# create a directory for moving avg files if it doesn't already exist

if (!dir.exists("moving_avg")){
  dir.create("moving_avg")
}else{
  print("dir exists")
}

for (filex in file_list) {
  

#upload data to dataframe, rename headers, make locus continuous, create subsets
depth <- read.table(filex, sep="\t", header=T)

# give new column names
names(depth)<-c("indexStart_indexStop_firstPos_withData_lastPos_withData_WinStart,WinStop","Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")

#*******dividing pi by no of sites to get proper pi

depth$tP<-depth$tP/depth$nSites

#**************************************************

#install.packages("data.table")
library(data.table)
#install.packages("zoo")
library(zoo)

Xdepth<-subset(depth, select = c("Chr", "WinCenter","tP"))

#genome coverage as sliding window
#Xdepth.average<-setDT(Xdepth)[, .(
#  window.start = rollapply(pos, width=win_size, by=win_size, FUN=min, align="left", partial=TRUE),
#  window.end = rollapply(pos, width=win_size, by=win_size, FUN=max, align="left", partial=TRUE),
#  coverage = rollapply(depth, width=win_size, by=win_size, FUN=mean, align="left", partial=TRUE)
#), .(Chr)]


# cal needed number of windows for moving average
mov_avg_no<-length(Xdepth$WinCenter)-mov_win_size+1

#subset needed dataset

moving_avg_df<-Xdepth[c(1:mov_avg_no),c(1:3)]

#add moving avg to dataframe
moving_avg_df$moving_avg<-rollmean(Xdepth$tP, k = mov_win_size)

#extract file name to save
saving_name<-tail(strsplit(filex, "/")[[1]],n=1)

#write.table(Xdepth.average, file=paste("./windowed_depth/",saving_name,".windowed",sep = ''), quote=FALSE, sep='\t', col.names = NA)

write.table(moving_avg_df, file=paste("./moving_avg/",saving_name,"_movingavg",sep = ''), quote=FALSE, sep='\t', col.names = NA)

}

```
Then plot

```R
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))

library(ggplot2)

# Delete previous plot before starting if it is in the same folder*****************


#remove scientific notation
options(scipen=999)

pop_name<-""

file_list<-grep(list.files(path="./moving_avg",full.names = TRUE), pattern='\\movingavg$', value=TRUE)

print(file_list)

# combine files adding pop and sex from file name
df <- do.call(rbind, lapply(file_list, function(x) cbind(read.table(x), file_name=tail(strsplit(x, "/")[[1]],n=1))))

# give new column names
#names(df)<-c("Chr","Start","End","Depth","file_name")

# Split name column into firstname and last name
require(stringr)
df[c( 'ind','sex','pop','info')] <- str_split_fixed(df$file_name, '_', 4)

# Subset Rows by column value
df<-df[df$Chr == 'Chr7',]

# Subset Rows by column value
df_females<-df[df$sex == 'F',]
df_males<-df[df$sex == 'M',]

# #create an empty column for means/males
# df_males$ind_mean<-"NA"
# 
# #add mean column males
# male_means<-tapply(df_males$moving_avg,df_males$ind, mean)
# male_mean_df<-as.data.frame(male_means)
# # DataFrame and put name of function rn
# male_mean_df <- tibble::rownames_to_column(male_mean_df, "ind")
# 
# #conditionally fill the column comparing values
# for (i in 1:length(male_mean_df$ind)) {
#   for (j in 1:length(df_males$moving_avg)) {
#     if (male_mean_df$ind[i]==df_males$ind[j]) {
#       df_males$ind_mean[j]<-male_mean_df$male_means[i]
#     }
#     j=j+1
#   }
# }
# 
# #add a column with std depth
# df_males$std_depth<-as.numeric(df_males$moving_avg)/as.numeric(df_males$ind_mean)
# 
# #create an empty column for means/females
# df_females$ind_mean<-"NA"
# 
# #add mean column females
# male_means<-tapply(df_females$moving_avg,df_females$ind, mean)
# male_mean_df<-as.data.frame(male_means)
# # DataFrame and put name of function rn
# male_mean_df <- tibble::rownames_to_column(male_mean_df, "ind")
# 
# #conditionally fill the column comparing values
# for (i in 1:length(male_mean_df$ind)) {
#   for (j in 1:length(df_females$moving_avg)) {
#     if (male_mean_df$ind[i]==df_females$ind[j]) {
#       df_females$ind_mean[j]<-male_mean_df$male_means[i]
#     }
#     j=j+1
#   }
# }
# 
# #add a column with std depth
# df_females$std_depth<-as.numeric(df_females$moving_avg)/as.numeric(df_females$ind_mean)

#list populations

pop_list<-unique(df$pop)

print(pop_list)

# #add fake data for absent populations*****
# 
# pops_to_add<-c("Liberia","Ivory","scaffold","mello")
# 
# #extract no of rows the lengnth of absent pops
# 
# fake_df<-head(df,length(pops_to_add))
# fake_df$pop<-pops_to_add
# fake_df$moving_avg<-rep(0,length(pops_to_add))
# 
# df<-rbind(df,fake_df)
# 
# pop_list<-unique(df$pop)
# 
# print(pop_list)

p<-ggplot(df_males,aes(x=WinCenter/1000000,y=moving_avg,color=sex))+
  geom_rect(data=df, inherit.aes=FALSE, aes(xmin=6.5, xmax=9, ymin=min(df$moving_avg),ymax=max(df$moving_avg)), fill="lightgrey", alpha=0.3)+
  #geom_point(alpha=0.1)+
  geom_point(size=0.5)+
  scale_color_manual(values = c("M"="blue","F"='red'))+
  geom_point(data = df_females,size=0.5,aes(x=WinCenter/1000000,y=moving_avg))+
  #scale_color_manual(labels = c("males", "females"), values = c("blue", "red"))+
  labs(title = pop_name)+
  #ylim(-0.05, 0.25)+
  xlim(1,20)+
  theme_bw()+
  xlab("Position - Chr 7 (MB)")+
  ylab(expression(pi))+
  theme_classic()+ 
  theme(legend.position="none")
  #geom_hline(yintercept = 1,color="lightgray")

#chgange order here
p+facet_wrap(~factor(pop,levels=c("Liberia","Sierra","Ivory","Ghana","Tad","Nigeria","Scaffold","Mellotropicalis","Cal")),ncol=1,labeller = as_labeller(c(Liberia='Liberia',Ghana='Ghana',Tad='Tadpoles',Nigeria='Nigeria',Sierra='Sierra Leone',Ivory='Ivory Coast',Scaffold='Reference Genome',Mellotropicalis='Mellotropicalis data',Cal="Calcaratus")))


ggsave(paste("20M_",pop_name,"_moving_avg_plot.pdf",sep = ""),width = 5, height = 10)

```
# PCA

merging all vcfs

```bash
#!/bin/sh
#SBATCH --job-name=abba
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=128gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load bcftools
module load vcftools

for i in *.vcf; do  bgzip -c $äiå > $äiå.gz ; done
for i in *.vcf.gz; do bcftools index -t $äiå ; done

vcf-concat combined_Chr10.g.vcf.gz_Chr10_GenotypedSNPs.vcf.gz_filtered.vcf.gz_selected.vcf.gz combined_Chr1.g.vcf.gz_Chr1_GenotypedSNPs.vcf.gz_filtered.vcf.gz_selected.vcf.gz combined_Chr2.g.vcf.gz_Chr2_GenotypedSNPs.vcf.gz_filtered.vcf.gz_selected.vcf.gz combined_Chr3.g.vcf.gz_Chr3_GenotypedSNPs.vcf.gz_filtered.vcf.gz_selected.vcf.gz combined_Chr4.g.vcf.gz_Chr4_GenotypedSNPs.vcf.gz_filtered.vcf.gz_selected.vcf.gz combined_Chr5.g.vcf.gz_Chr5_GenotypedSNPs.vcf.gz_filtered.vcf.gz_selected.vcf.gz combined_Chr6.g.vcf.gz_Chr6_GenotypedSNPs.vcf.gz_filtered.vcf.gz_selected.vcf.gz combined_Chr7.g.vcf.gz_Chr7_GenotypedSNPs.vcf.gz_filtered.vcf.gz_selected.vcf.gz combined_Chr8.g.vcf.gz_Chr8_GenotypedSNPs.vcf.gz_filtered.vcf.gz_selected.vcf.gz combined_Chr9.g.vcf.gz_Chr9_GenotypedSNPs.vcf.gz_filtered.vcf.gz_selected.vcf.gz combined_Scafs.g.vcf.gz_Scafs_GenotypedSNPs.vcf.gz_filtered.vcf.gz_selected.vcf.gz ö bgzip -c > combined_data_2023_Aug.vcf.gz
```
then removed not needed samples

```bash
vcftools --remove-indv all_calcaratus_sorted.bam --remove-indv mello_GermSeq_sorted.bam --gzvcf combined_data_2023_Aug.vcf.gz --out no_mello_no_lib_whole_genome.vcf.gz --recode
```




