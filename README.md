# Transcriptomic_analysis_of_*Saccharolobus_islandicus*
## 1.Data source
The experimental design of the study focused on understanding the gene expression of Saccharolobus islandicus relative to its different cell-cycle stages.  Cultures of S. islandicus were made to begin the cell cycle at once by a 6-hour treatment with acetic acid. After the removal of the acetic acid, the cells resumed their cell cycle. Total RNA was extracted from 3 samples at three specific time points after synchronization:
- sample 1:  2 hours and 30 minutes (in M-G1 phase)
- sample 2: 4 hours (in S phase)
- sample 3: 6 hours( in G2 phase)
## 2.Obtaining the data
The reads (Sequenced with Illumina NextSeq 2000) were obtaiend from the NCBI database along with a ```.fasta``` reference and ```gunzip *.fastq.gz``` was used to unzip the file. 
## 3.Workflow
This workflow explains the sequence processing of subsamples
### 3.1 Mapping reads to reference 
```module load bowtie2-2.4.1, module load samtools-1.18,module load subread-2.0.8```  was ran to load the tool into the enivironmnet.
##### Loop through each sample
```samples=(2h30_13_S26_R1_001.fastq  4h_1_S31_R1_001.fastq  6h_10_S53_R1_001.fastq)
for i in "${samples[@]}"; do
    bowtie2 -x ref_index -U "${i}.fastq" -S "${i}.sam"
    samtools view -b "${i}.sam" > "${i}.bam"
    rm "${i}.sam"
done
```
### 3.2 FeatureCounts for quantifying reads
```
featureCounts -a sequence.gtf -F GTF -o counts_2.txt -T 14 *.bam
```



