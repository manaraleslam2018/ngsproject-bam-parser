Project steps 

1- I started with take an overview in this paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6052005/) to know general information about dogs their breeds and their genome assembly. 

2- I tried to download the sra reads (prefetch SRR7120302) from sra (https://www.ncbi.nlm.nih.gov/sra/SRX4041929[accn]#) but the download failed so I started to work with already dog sequences we used in the lab. 
'''
mkdir /home/manar/ngs1_project && cd /home/manar/ngs1_project
mkdir /home/manar/ngs1_project/samples/ && cd /home/manar/ngs1_project/samples/
wget https://de.cyverse.org/dl/d/3CE425D7-ECDE-46B8-AB7F-FAF07048AD42/samples.tar.gz
tar xvzf samples.tar.gz
'''
3- download the GTF and fasta file reference :

I don't know which one gtf annotation I should choose (ftp://ftp.ensembl.org/pub/release-96/gtf/canis_familiaris/). Therefore I searchd to know the difference between (File:Canis_familiaris.CanFam3.1.96.chr.gtf.gz , Canis_familiaris.CanFam3.1.96.gtf.gz) and I found this (https://www.biostars.org/p/217700/#241695) and 
This link explained that the gtf without 'chr' contains annotations for genes on unplaced or unlocalized contigs, while the one with 'chr' only contains annotation for assembled chromosomes, both of them have no prefix 'chr' in chromosome name and based on their recommendation to download the one without 'chr'  if you do not want to loose information about any annotated gene.

'''
mkdir /home/manar/ngs1_project/reference && cd /home/manar/ngs1_project/reference
wget ftp://ftp.ensembl.org/pub/release-96/gtf/canis_familiaris/Canis_familiaris.CanFam3.1.96.gtf.gz
gunzip Canis_familiaris.CanFam3.1.96.gtf.gz
wget https://de.cyverse.org/dl/d/A9330898-FC54-42A5-B205-B1B2DC0E91AE/dog_chr5.fa.gz
gunzip dog_chr5.fa.gz
'''

4- Run the FASTQC for each read end
'''
mkdir/home/manar/ngs1_project/fastqc && cd /home/manar/ngs1_project/fastqc
cp /home/manar/ngs1_project/sampels/samples/BD143_TGACCA_L005_R1_001.pe.fq.gz .
cp /home/manar/ngs1_project/sampels/samples/BD143_TGACCA_L005_R2_001.pe.fq.gz .
cp /home/manar/ngs1_project/sampels/samples/BD143_TGACCA_L006_R1_001.pe.fq.gz .
cp /home/manar/ngs1_project/sampels/samples/BD143_TGACCA_L006_R2_001.pe.fq.gz .

for f in /home/manar/ngs1_project/fastqc/*.pe.fq.gz;do fastqc -t 1 -f fastq -noextract $f;done
'''
5- trimming 
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

6- hisat alignment
index your genome

'''
mkdir home/manar/ngs1_project/hisat_alignmet && cd /home/manar/ngs1_project/hisat_alignmet 
mkdir hisatIndex && cd hisatIndex

ln -s /home/manar/ngs1_project/reference/dog_chr5.fa .
ln -s /home/manar/ngs1_project/reference/Canis_familiaris.CanFam3.1.96.gtf .

hisat2_extract_splice_sites.py /home/manar/ngs1_project/reference/Canis_familiaris.CanFam3.1.96.gtf > splicesites.tsv
hisat2_extract_exons.py /home/manar/ngs1_project/reference/Canis_familiaris.CanFam3.1.96.gtf > exons.tsv
hisat2-build -p 1 --ss splicesites.tsv --exon exons.tsv /home/manar/ngs1_project/reference/dog_chr5.fa dog_chr5_indexed
'''
sequence alignment
'''
cd ..

for r in 5 6 
do
R1="/home/manar/ngs1_project/trimmed/BD143_TGACCA_L00${r}_R1_001.pe.trim.fastq.gz"
R2="/home/manar/ngs1_project/trimmed/BD143_TGACCA_L00${r}_R2_001.pe.trim.fastq.gz"
hisat2 -p 1 -x /home/manar/ngs1_project/hisat_alignmet/hisatIndex/dog_chr5_indexed --dta --rna-strandness RF -1 $R1 -2 $R2 -S hisat_aligment${r}.sam
done
'''

stats
'''
for r in 5 6
do
    samtools flagstat /home/manar/ngs1_project/hisat_alignmet/hisat_aligment${r}.sam > haisat_alignment_sample${r}_stats.out
done
'''
the stats file showed there is no any supplementary or discorcdant aligment reads so I started to search for new organism to work in.





