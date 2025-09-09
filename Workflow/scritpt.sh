wdk="Folder/sample"
cd "$wdk"

# unzip all the .gz files
#gunzip *.gz

module load bowtie2-2.4.1
module load samtools-1.18
module load subread-2.0.8
samples=(2h30_13_S26_R1_001.fastq  4h_1_S31_R1_001.fastq  6h_10_S53_R1_001.fastq)

# Loop through each sample
for i in "${samples[@]}"; do
    bowtie2 -x ref_index -U "${i}.fastq" -S "${i}.sam"
    samtools view -b "${i}.sam" > "${i}.bam"
    rm "${i}.sam"
done

#featureCounts
#featureCounts -a sequence.gtf -F GTF -o counts_2.txt -T 14 *.bam