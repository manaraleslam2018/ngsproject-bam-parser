https://en.wikipedia.org/wiki/List_of_sequenced_plant_genomes#Press_releases_announcing_sequencing
this link was useful in selecting an organism.

https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1790382
SRR1790382

download samples
'''
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1790382/SRR1790382.sra
'''
to dwonload gtf 
'''
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-43/gff3/corchorus_capsularis/Corchorus_capsularis.CCACVL1_1.0.43.gff3.gz
'''
download fasta file refrence (ftp://ftp.ensemblgenomes.org/pub/plants/release-43/fasta/corchorus_capsularis/dna/)
'''
ftp://ftp.ensemblgenomes.org/pub/plants/release-43/fasta/corchorus_capsularis/dna/Corchorus_capsularis.CCACVL1_1.0.dna.nonchromosomal.fa.gz
gunzip Corchorus_capsularis.CCACVL1_1.0.dna.nonchromosomal.fa.gz
'''

to choose the fasta file i can use for alignment 

http://genomespot.blogspot.com/2015/06/mapping-ngs-data-which-genome-version.html
https://bioinformatics.stackexchange.com/questions/540/what-ensembl-genome-version-should-i-use-for-alignments-e-g-toplevel-fa-vs-p

Repeat Masking?
The short answer is that for the purpose of mapping NGS reads, you don't want to use a RepeatMasked genome. Repeat-Masking is a process to identify all repetitive parts of the genome such as LINE-1 and Alu and replace these sequences with "N". While this may sound like a good idea, it's not because mapping software are specifically designed to handle reads from these repeats and have sophisticated heuristics to either discard those alignments or find all possible alignments. So select a build that has no masking and avoid any with "rm" in the name (ie: Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa.gz). The files with "sm" in the name are "soft masked" which means instead of the repeat regions being given Ns, they are converted to lower-case. These are OK to use because most major aligners recognise lower-case letters as valid bases.

extract fastq from sra data 
'''
fastq-dump --outdir . --gzip --split-3 SRR1790382.sra
'''
fastqc 
for f in /home/manar/ngs1_project/plant_samples/fastqc/*.fastq.gz;do fastqc -t 1 -f fastq -noextract $f;done

/home/manar/ngs1_project/plant_samples/fastqc/SRR1790382_1.fastq.gz


4- Trimming
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

HISAT allignment 
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
 Proper discordant read extraction - samblaster vs. samtools #193 
https://github.com/arq5x/lumpy-sv/issues/193
'''
samtools view -b -F 1294 plant_hisat_aligment.bam_sorted.bam > plant1294_sample.discordants_sorted.bam

samtools view -b -F 3854 plant_hisat_aligment.bam_sorted.bam > plant3854_sample.discordants_sorted.bam
'''
