# Transcriptomic_analysis_of_Saccharolobus_islandicus
## 1.Data source
The experimental design of the study focused on understanding the gene expression of Saccharolobus islandicus relative to its different cell-cycle stages.  Cultures of S. islandicus were made to begin the cell cycle at once by a 6-hour treatment with acetic acid. After the removal of the acetic acid, the cells resumed their cell cycle. Total RNA was extracted from 3 samples at three specific time points after synchronization:
- sample 1:  2 hours and 30 minutes (in M-G1 phase)
- sample 2: 4 hours (in S phase)
- sample 3: 6 hours( in G2 phase)
## 2.Obtaining the data
The reads (Sequenced with Illumina NextSeq 2000) were obtaiend from the NCBI database and ```gunzip *.fastq.gz``` was used to unzip the file. 
## 3.Workflow
This workflow exaplains the hybrid de novo assembly of *Pseudomonas aeruginosa* YK01 (script stored in workflow/ )
### 3.1 Quality check
The tool ```fastqc``` was used to check the read quality . ```module load fastqc-0.11.7```  was ran to load the tool into the enivironmnet and the tool was ran:
<pre> fastqc SRR30916324_1.fastq
fastqc SRR30916324_2.fastq</pre>
Output is an ```.html``` file that you can view in your web browser. 
Trimming of low quality was then performed with ```trimmomatic-0.36``` for the paired-end reads:
<pre> java -jar $TRIMMOMATIC PE SRR30916324_1.fastq SRR30916324_2.fastq SRR30916324_1_paired.fastq SRR30916324_1_unpaired.fastq SRR30916324_2_paired.fastq SRR30916324_2_unpaired.fastq SLIDINGWINDOW:4:28 MINLEN:50 </pre>
Output is a trimmed ```.fastq``` file.
### 3.2 Hybrid de novo assembly
In this step, long reads were assembled and refined using short reads to construct the genome
```module load unicycler```  was ran to load the tool into the enivironmnet. ```Unicycler``` is the assmbely tool which reqiured ```module load racon``` and ```module load samtools-1.7 ``` in this instance.
<pre>unicycler -1 SRR30916324_1_paired.fastq -2 SRR30916324_2_paired.fastq -l SRR30916323.fastq -o Assembled_output --threads 8 --no_pilon</pre>
Output is a complete genome as a .fasta file
### 3.3 Checking completeness of the genome
Using  ```CheckM``` we are able to analyis the genome for Completeness, Contamination and Strain heterogeneity.```module load checkm```  was ran to load the tool into the enivironmnet.
<pre>checkm lineage_wf -t 8 -x fasta Assembled_output/ checkm_output/</pre>
Output is  hmmer.analyze.txt (for completeness/contamination estimation)  hmmer.tree.txt (for phylogenetic tree marker genes)




