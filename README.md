# suggestions : compiled by Victor Jiménez-Vásquez (vr.jimenez.vs@gmail.com) 

# non_human_reads : a collection of codes to get non-human reads from fastq-files by mapping against human genome (Homo sapiens reference genome GRCh38.p14, https://www.ncbi.nlm.nih.gov/search/all/?term=human%20genome) 

## The code 
```r
#1# preparar las instrucciones generales#
for r1 in *fastq.gz
do
prefix=$(basename $r1 _L001_R1_001.fastq.gz)
r2=${prefix}_L001_R2_001.fastq.gz

#2# instrucciones para generar el archivo .bam#
bwa mem -t 15 human.fasta $r1 $r2 > ${prefix}_uno.sam ;
samtools view -@ 15 -bS -T human.fasta ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -@ 15 -n ${prefix}_unoa.bam -o ${prefix}_dosa.bam ;
samtools fixmate -@ 15 -m ${prefix}_dosa.bam ${prefix}_tresa.bam ;
samtools sort -@ 15 ${prefix}_tresa.bam -o ${prefix}_cuatroa.bam ;
samtools markdup -@ 15 ${prefix}_cuatroa.bam ${prefix}.bam ;
samtools index -@ 15 ${prefix}.bam ;

#3# remover los archivos intermediarios#
rm ${prefix}_uno.bam ${prefix}_uno.sam ${prefix}_unoa.bam ${prefix}_dosa.bam ${prefix}_tresa.bam ${prefix}_cuatroa.bam ;

#4# extraer del archivo .bam, todos los reads alineados con la referencia y generar dos fastq (f y r)
bam2fastq --no-aligned -o ${prefix}_unal#.fastq ${prefix}.bam ;
done ;

#5# mover estos nuevos fastq a una carpeta nueva para prerar el ensamblaje de novo#
mkdir aligned_for_denovo ;
mv *.fastq aligned_for_denovo ;
cd aligned_for_denovo ;

#6# generar archivos fastq solo con reads pareados#
for r1 in *fastq
do
prefix=$(basename $r1 _unal_1.fastq)
r2=${prefix}_unal_2.fastq
fastq_pair $r1 $r2 ; 
rm *.single.fq ;
done ;

#7# realizar en ensamblaje de novo con spades#
for r1 in *fq
do
prefix=$(basename $r1 _unal_1.fastq.paired.fq)
r2=${prefix}_unal_2.fastq.paired.fq
spades --pe1-1 $r1 --pe1-2 $r2 --careful -t 15 -m 15 -o ${prefix}_spades ;
mv ${prefix}_spades/scaffolds.fasta ${prefix}_spades/${prefix}_spades_scaffolds.fasta ;
mv ${prefix}_spades/${prefix}_spades_scaffolds.fasta . ;
done ;
rmdir *fq_spades ;
mv *scaffolds.fasta .. ;
exit
```
