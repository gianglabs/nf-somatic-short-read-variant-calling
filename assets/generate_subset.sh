work/d4/48c702608dc99a24c089cba2e2387c/HCC1395N.1_aligned.bam
work/67/418a16d043223572001f915ac1363c/HCC1395T.1_aligned.bam

# raw fastq
samtools view -b -F 4 HCC1395N_merged.bam chr22:18000000-20000000 > subset.bam
samtools sort -n -o subset_name.bam subset.bam
samtools fastq \
  -1 HCC1395N_subset_R1.fastq.gz \
  -2 HCC1395N_subset_R2.fastq.gz \
  -0 /dev/null \
  -s /dev/null \
  -n subset_name.bam


samtools view -b -F 4 HCC1395T_merged.bam chr22:18000000-20000000 > subset.bam
samtools sort -n -o subset_name.bam subset.bam
samtools fastq \
  -1 HCC1395T_subset_R1.fastq.gz \
  -2 HCC1395T_subset_R2.fastq.gz \
  -0 /dev/null \
  -s /dev/null \
  -n subset_name.bam

# subset genome
# download only chr22
aws s3 ls --no-sign-request s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/Chromosomes/chr22.fasta .
# index
samtools faidx chr22.fasta
# get fasta dict
picard CreateSequenceDictionary \
    R=chr22.fasta \
    O=chr22.dict