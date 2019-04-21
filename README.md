# ngsproject-bam-parser
At first I worked on dog samples (DNA) all steps :
1- Run the FASTQC for each read end
2- trimming
3- I choose hisat alignment because it is splice aware aligner.
4- convert SAM file to BAM file 
I make stat on bam file, I didn't find (supplementry, mismatch, discordant reads) therefore I searched for another organism.

I repeat all steps on plant ( Corchorus capsularis CVL-1).
after alignment I found discordant reads so I made BAM parser by samtools and axtract these discordant reads.
I tried to found the cause to be discordant:
after searching I found the definition of discordant: If the distance between R1 and R2 is outside 500bp (+/- 1SD) or they relative orientation is not Forward/Reverse thy will be called as discordant pairs.
so I made feature count to link between the features and my reads to know if these reads were too short or too long to make the distance out of range. 
also to know if these reads clustered in a specefic region in the genome or dipersed.
