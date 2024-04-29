source /apps/profiles/modules_asax.sh.dyn
module load star/2.7.6a
module load zlib/1.2.7 


#indexing
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir /scratch/aubclsc0336/Lizard_index \
--genomeFastaFiles /home/aubclsc0336/SceUnd1.0_top24.fasta \
--sjdbGTFfile /home/aubclsc0336/SceUnd1.0_top24.gtf \
--sjdbOverhang 149

#mapping
#you need to create a file STAR within the PracticeRNAseq_Full_Script to run it
STAR --genomeDir /scratch/aubclsc0336/Lizard_index/ \
--runThreadN 6 \
--readFilesIn  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825925_1_paired.fastq  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825925_2_paired.fastq \
--outFileNamePrefix  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/STAR/SRR6825925 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \

STAR --genomeDir /scratch/aubclsc0336/Lizard_index/ \
--runThreadN 6 \
--readFilesIn  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825926_1_paired.fastq  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825926_2_paired.fastq \
--outFileNamePrefix  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/STAR/SRR6825926 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \

STAR --genomeDir /scratch/aubclsc0336/Lizard_index/ \
--runThreadN 6 \
--readFilesIn  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825928_1_paired.fastq  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825928_2_paired.fastq \
--outFileNamePrefix  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/STAR/SRR6825928 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \

STAR --genomeDir /scratch/aubclsc0336/Lizard_index/ \
--runThreadN 6 \
--readFilesIn  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825930_1_paired.fastq  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825930_2_paired.fastq \
--outFileNamePrefix  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/STAR/SRR6825930 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \

STAR --genomeDir /scratch/aubclsc0336/Lizard_index/ \
--runThreadN 6 \
--readFilesIn  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825932_1_paired.fastq  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825932_2_paired.fastq \
--outFileNamePrefix  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/STAR/SRR6825932 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \

STAR --genomeDir /scratch/aubclsc0336/Lizard_index/ \
--runThreadN 6 \
--readFilesIn  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825933_1_paired.fastq  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825933_2_paired.fastq \
--outFileNamePrefix  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/STAR/SRR6825933 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \

STAR --genomeDir /scratch/aubclsc0336/Lizard_index/ \
--runThreadN 6 \
--readFilesIn  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825934_1_paired.fastq  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825934_2_paired.fastq \
--outFileNamePrefix  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/STAR/SRR6825934 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \

STAR --genomeDir /scratch/aubclsc0336/Lizard_index/ \
--runThreadN 6 \
--readFilesIn  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825935_1_paired.fastq  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825935_2_paired.fastq \
--outFileNamePrefix  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/STAR/SRR6825935 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \

STAR --genomeDir /scratch/aubclsc0336/Lizard_index/ \
--runThreadN 6 \
--readFilesIn  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825938_1_paired.fastq  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825938_2_paired.fastq \
--outFileNamePrefix  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/STAR/SRR6825938 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \

STAR --genomeDir /scratch/aubclsc0336/Lizard_index/ \
--runThreadN 6 \
--readFilesIn  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825939_1_paired.fastq  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825939_2_paired.fastq \
--outFileNamePrefix  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/STAR/SRR6825939 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \

STAR --genomeDir /scratch/aubclsc0336/Lizard_index/ \
--runThreadN 6 \
--readFilesIn  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825940_1_paired.fastq  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825940_2_paired.fastq \
--outFileNamePrefix  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/STAR/SRR6825940 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \

STAR --genomeDir /scratch/aubclsc0336/Lizard_index/ \
--runThreadN 6 \
--readFilesIn  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825941_1_paired.fastq  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/CleanData/SRR6825941_2_paired.fastq \
--outFileNamePrefix  /scratch/aubclsc0336/PracticeRNAseq_Full_Script/STAR/SRR6825941 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \




