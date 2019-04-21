https://en.wikipedia.org/wiki/List_of_sequenced_plant_genomes#Press_releases_announcing_sequencing
this link was useful in selecting an organism with low quality genome (Corchorus capsularis CVL-1).

1- download samples from SRA
'''
mkdir /home/manar/ngs1_project/plant_samples && cd /home/manar/ngs1_project/plant_samples
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1790382/SRR1790382.sra
'''
2- extract fastq from sra data 
'''
fastq-dump --outdir . --gzip --split-3 SRR1790382.sra
'''
3- dwonload gtf reference file 
'''
mkdir /home/manar/ngs1_project/plant_refernce && cd /home/manar/ngs1_project/plant_refernce
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-43/gff3/corchorus_capsularis/Corchorus_capsularis.CCACVL1_1.0.43.gff3.gz
'''
4- download fasta file refrence (ftp://ftp.ensemblgenomes.org/pub/plants/release-43/fasta/corchorus_capsularis/dna/)
'''
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-43/fasta/corchorus_capsularis/dna/Corchorus_capsularis.CCACVL1_1.0.dna.nonchromosomal.fa.gz
gunzip Corchorus_capsularis.CCACVL1_1.0.dna.nonchromosomal.fa.gz
'''

I was confused between different fasta files to use as a refernec. To choose the fasta file which I can use for alignment,
I searched and found that I should avoid  any file with "rm" in the name but I can use the files with "sm" in the name are "soft masked" or without masking. 
the answer available in these websites:
http://genomespot.blogspot.com/2015/06/mapping-ngs-data-which-genome-version.html
https://bioinformatics.stackexchange.com/questions/540/what-ensembl-genome-version-should-i-use-for-alignments-e-g-toplevel-fa-vs-p

5- fastqc 
'''
mkdir /home/manar/ngs1_project/plant_samples/fastqc/ && cd /home/manar/ngs1_project/plant_samples/fastqc 
cp /home/manar/ngs1_project/plant_samples/SRR1790382_1.fastq.gz .
cp /home/manar/ngs1_project/plant_samples/SRR1790382_2.fastq.gz .
for f in /home/manar/ngs1_project/plant_samples/fastqc/*.fastq.gz;do fastqc -t 1 -f fastq -noextract $f;done
'''

6- Trimming
'''
mkdir ~/ngs1_project/plant_trimmed && cd ~/ngs1_project/plant_trimmed
 
f1="/home/manar/ngs1_project/plant_samples/SRR1790382_1.fastq.gz";
f2="/home/manar/ngs1_project/plant_samples/SRR1790382_2.fastq.gz";
newf1="/home/manar/ngs1_project/plant_trimmed/plant_sample_r1.pe.trim.fastq.gz"; 
newf2="/home/manar/ngs1_project/plant_trimmed/plant_sample_r2.pe.trim.fastq.gz"; 
newf1U="/home/manar/ngs1_project/plant_trimmed/plant_sample_r1.se.trim.fastq.gz";
newf2U="/home/manar/ngs1_project/plant_trimmed/plant_sample_r2.se.trim.fastq.gz";
adap="/home/manar/anaconda3/envs/ngs1/share/trimmomatic-0.39-0/adapters"; 
trimmomatic PE -threads 1 -phred33 -trimlog trimLogFile -summary statsSummaryFile  $f1 $f2 $newf1 $newf1U $newf2 $newf2U ILLUMINACLIP:$adap/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:15 MINLEN:36
'''

7- HISAT allignment 
index your genome
'''
mkdir /home/manar/ngs1_project/plant_hisat_alignment/hisatIndex && cd /home/manar/ngs1_project/plant_hisat_alignment/hisatIndex

ln -s /home/manar/ngs1_project/plant_refernce/Corchorus_capsularis.CCACVL1_1.0.dna.nonchromosomal.fa .
ln -s /home/manar/ngs1_project/plant_refernce/Corchorus_capsularis.CCACVL1_1.0.43.gff3 .
hisat2_extract_splice_sites.py /home/manar/ngs1_project/plant_refernce/Corchorus_capsularis.CCACVL1_1.0.43.gff3 > splicesites.tsv
hisat2_extract_exons.py /home/manar/ngs1_project/plant_refernce/Corchorus_capsularis.CCACVL1_1.0.43.gff3 > exons.tsv
hisat2-build -p 1 --ss splicesites.tsv --exon exons.tsv /home/manar/ngs1_project/plant_refernce/Corchorus_capsularis.CCACVL1_1.0.dna.nonchromosomal.fa Corchorus_capsularis.CCACVL1_1.0.dna.nonchromosomal_indexed
'''

sequence alignment
'''
cd ~/ngs1_project/plant_hisat_alignment
R1="/home/manar/ngs1_project/plant_trimmed/plant_sample_r1.pe.trim.fastq.gz"
R2="/home/manar/ngs1_project/plant_trimmed/plant_sample_r2.pe.trim.fastq.gz"
hisat2 -p 1 -x hisatIndex/Corchorus_capsularis.CCACVL1_1.0.dna.nonchromosomal_indexed --rna-strandness RF -1 $R1 -2 $R2 -S plant_hisat_aligment.sam
'''

stats 

'''
samtools flagstat /home/manar/ngs1_project/plant_hisat_alignment/plant_hisat_aligment.sam > plant_haisat_alignment_stats.out

'''

Convert the SAM files from Haisat aligment into BAM file 
'''
samtools view -bS /home/manar/ngs1_project/plant_hisat_alignment/plant_hisat_aligment.sam > /home/manar/ngs1_project/plant_hisat_alignment/plant_hisat_aligment.bam
samtools sort /home/manar/ngs1_project/plant_hisat_alignment/plant_hisat_aligment.bam -o /home/manar/ngs1_project/plant_hisat_alignment/plant_hisat_aligment.bam_sorted.bam

'''
Bam parser development
1- To develope Bam parser, I should have study well AAM\BAM specification through this link 
http://chagall.med.cornell.edu/galaxy/references/SAM_BAM_Specification.pdf

2- Bam parsers include ( pysam or pysbm or samtools) I used samtools (http://www.htslib.org/doc/samtools.html)

3- Proper discordant read extraction by FLAG number using samtools
https://github.com/arq5x/lumpy-sv/issues/193
Some websites mentioned -F 1294 others -F 3854 I searched for the explanation by (https://broadinstitute.github.io/picard/explain-flags.html)
I believe that 3854 is the correct from my point of view, Also, I tried both and give thwe same results.

'''
samtools view -b -F 1294 plant_hisat_aligment.bam_sorted.bam > plant1294_sample.discordants_sorted.bam

samtools view -b -F 3854 plant_hisat_aligment.bam_sorted.bam > plant3854_sample.discordants_sorted.bam
'''
