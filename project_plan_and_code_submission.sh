I started with take a look in this paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6052005/) to know general information aout dogs their breeds. 
Afterthat I tried to download the sra reads (prefetch SRR7120302) from sra (https://www.ncbi.nlm.nih.gov/sra/SRX4041929[accn]#). 
Also I searched for refrence genome for dogs and download the fasta file from (ftp://ftp.ensembl.org/pub/release-96/fasta/canis_familiaris/dna/)
and download the gtf annotation for them but I don't know which one gtf annotation I should choose (ftp://ftp.ensembl.org/pub/release-96/gtf/canis_familiaris/)

I searchd to know the difference between them I found this (https://www.biostars.org/p/217700/#241695)

The one without 'chr' contains annotations for genes on unplaced or unlocalized contigs, while the one with 'chr' only contains annotation for assembled chromosomes, both of them have no prefix 'chr' in chromosome name and based on their recommendation to download the one without 'chr'  if you do not want to loose information about any annotated gene.

wget ftp://ftp.ensembl.org/pub/release-96/gtf/canis_familiaris/Canis_familiaris.CanFam3.1.96.gtf.gz
gunzip Canis_familiaris.CanFam3.1.96.gtf.gz

the sample download failed because of the big data and low internet connection so I thought to download small sample and and small reference and tried the command on it then apply to large dataset 

https://github.com/drtamermansour/nu-ngs02/blob/master/Crash_variant_calling.md

wget https://de.cyverse.org/dl/d/3CE425D7-ECDE-46B8-AB7F-FAF07048AD42/samples.tar.gz
tar xvzf samples.tar.gz

wget https://de.cyverse.org/dl/d/A9330898-FC54-42A5-B205-B1B2DC0E91AE/dog_chr5.fa.gz
gunzip dog_chr5.fa.gz

cd /home/manar/ngs1_project/samples/
for f in /home/manar/ngs1_project/samples/*.pe.fq.gz;do fastqc -t 1 -f fastq -noextract $f;done



hisat alignment
index your genome
'''

mkdir ~/home/manar/ngs1_project/hisat_alignmet && cd ~/home/manar/ngs1_project/hisat_alignmet 
mkdir hisatIndex && cd hisatIndex

ln -s /home/manar/ngs1_project/reference/dog_chr5.fa .
ln -s /home/manar/ngs1_project/reference/Canis_familiaris.CanFam3.1.96.gtf .

hisat2_extract_splice_sites.py /home/manar/ngs1_project/reference/Canis_familiaris.CanFam3.1.96.gtf > splicesites.tsv
hisat2_extract_exons.py /home/manar/ngs1_project/reference/Canis_familiaris.CanFam3.1.96.gtf > exons.tsv
hisat2-build -p 1 --ss splicesites.tsv --exon exons.tsv /home/manar/ngs1_project/reference/dog_chr5.fa dog_chr5_indexed
'''
sequence alignment

'''
cd ~/ngs1_project/hisat_align

for r in 5 6 
do
R1="/home/manar/ngs1_project/trimmed/BD143_TGACCA_L00${r}_R1_001.pe.trim.fastq.gz"
R2="/home/manar/ngs1_project/trimmed/BD143_TGACCA_L00${r}_R2_001.pe.trim.fastq.gz"
hisat2 -p 1 -x /home/manar/ngs1_project/hisat_alignmet/hisatIndex/dog_chr5_indexed --dta --rna-strandness RF -1 $R1 -2 $R2 -S hisat_aligment${r}.sam
done
'''
/home/manar/ngs1_project/trimmed/BD143_TGACCA_L005_R2_001.pe.trim.fastq.gz
stats

'''
for r in 5 6
do
    samtools flagstat /home/manar/ngs1_project/hisat_alignmet/hisat_aligment${r}.sam > haisat_alignment_sample${r}_stats.out
done
'''


trimming 
'''
mkdir /ngs1_project/trimmed && cd /ngs1_project/trimmed 
for r in 5 6 ;do 
f1="/home/manar/ngs1_project/sampels/samples/BD143_TGACCA_L00${r}_R1_001.pe.fq.gz";
f2="/home/manar/ngs1_project/sampels/samples/BD143_TGACCA_L00${r}_R2_001.pe.fq.gz";
newf1="/home/manar/ngs1_project/trimmed/BD143_TGACCA_L00${r}_R1_001.pe.trim.fastq.gz"; 
newf2="/home/manar/ngs1_project/trimmed/BD143_TGACCA_L00${r}_R2_001.pe.trim.fastq.gz"; 
newf1U="/home/manar/ngs1_project/trimmed/BD143_TGACCA_L00${r}_R1_001.se.trim.fastq.gz";
newf2U="/home/manar/ngs1_project/trimmed/BD143_TGACCA_L00${r}_R2_001.se.trim.fastq.gz";
adap="/home/manar/anaconda3/envs/ngs1/share/trimmomatic-0.39-0/adapters"; 
trimmomatic PE -threads 1 -phred33 -trimlog trimLogFile -summary statsSummaryFile${r} $f1 $f2 $newf1 $newf1U $newf2 $newf2U ILLUMINACLIP:$adap/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:15 MINLEN:36;
done
'''

