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
angsd -b pop1  -anc ${refgenome} -out ${outdir}/${sample}_pop1 -dosaf 1 -gl 1 -r Chr7:1-50000000
angsd -b pop2  -anc ${refgenome} -out ${outdir}/${sample}_pop2 -dosaf 1 -gl 1 -r Chr7:1-50000000

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


