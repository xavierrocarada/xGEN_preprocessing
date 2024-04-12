# merging fastq files from different lanes

```
#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=00:30:00
#SBATCH --mem=16GB

cat LP122_1_S154_L001_R1_001.fastq.gz LP122_1_S154_L002_R1_001.fastq.gz LP122_1_S52_L001_R1_001.fastq.gz LP122_1_S52_L002_R1_001.fastq.gz LP122_1_S52_L003_R1_001.fastq.gz LP122_1_S52_L004_R1_001.fastq.gz > LP122_1.R1.fastq.gz
cat LP122_1_S154_L001_R2_001.fastq.gz LP122_1_S154_L002_R2_001.fastq.gz LP122_1_S52_L001_R2_001.fastq.gz LP122_1_S52_L002_R2_001.fastq.gz LP122_1_S52_L003_R2_001.fastq.gz LP122_1_S52_L004_R2_001.fastq.gz > LP122_1.R2.fastq.gz
```

# construct unmapped bam

```
#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=12:00:00
#SBATCH --mem=32GB

module purge
module use /apps/modules/all
module load Java/1.8.0_121
module load picard/2.22.3-Java-1.8.0_121

for d in LP122_* ; do

	if [ -d "$d" ]; then
			
		(cd "$d" && 

    java -jar $EBROOTPICARD/picard.jar FastqToSam \
    FASTQ="$d".R1.fastq.gz \
    FASTQ2="$d".R2.fastq.gz \
    O="$d".bam \
    SM="$d"
    
    )
	fi
done
```

# extract UMIs and add RX tag

Install fgbio in conda.  Setup miniconda as https://github.com/ACAD-UofA/Bioinformatics-Wiki/wiki/Setup-Miniconda
To install new softwares, such as fastp or multiqc:
```
conda activate
conda create -n fastp/multiqc
conda install -c bioconda fastp / conda install -c bioconda multiqc
```
Script:
```
#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=6:00:00
#SBATCH --mem=16GB

conda activate fgbio

for f in *.unmapped.bam ; do

fgbio ExtractUmisFromBam --input=${f} --output=${f}.withUMI.bam --read-structure=8M+T 8M+T --molecular-index-tags=ZA ZB --single-tag=RX

done
```
# align reads


```
#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=72:00:00
#SBATCH --mem=64GB

module purge
module use /apps/modules/all
module load BWA/0.7.15-foss-2017a
module load Java/1.8.0_121
module load picard/2.22.3-Java-1.8.0_121

for f in *.withUMI.bam ; do

java -jar $EBROOTPICARD/picard.jar SamToFastq I=${f} F=/dev/stdout \
INTERLEAVE=true \
  | bwa mem -p -t 4 /hpcfs/users/a1782219/ref/human_g1k_v37_decoy.fasta.gz /dev/stdin \
  | java -jar $EBROOTPICARD/picard.jar MergeBamAlignment \
    UNMAPPED=${f} ALIGNED=/dev/stdin O=${f}.mapped.bam R=/hpcfs/users/a1782219/ref/human_g1k_v37_decoy.fasta.gz \
    SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=1 \
    ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

done
```

# markdup by UMI

```
#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=6:00:00
#SBATCH --mem=32GB

module purge
module use /apps/modules/all
module load Java/1.8.0_121
module load picard/2.22.3-Java-1.8.0_121

for f in *.mapped.bam ; do

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
  I=${f} \
  O=${f}.markdup.bam \
  M=${f}.markdup.metrics.txt \
  BARCODE_TAG=RX 

done
```

# rmdup

```
#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=2:00:00
#SBATCH --mem=16GB

module purge
module use /apps/modules/all
module load SAMtools/1.9-foss-2016b

for f in *.markdup.bam; do 

samtools view -b -u \
        -F `python -c 'print(hex(0x400))'` \
        ${f} \
        | samtools sort -O bam -o ${f}.rmdup.bam

done
```
# correct errors UMI seq
```
#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=6:00:00
#SBATCH --mem=16GB

conda activate fgbio

for f in *.unmapped.bam ; do

fgbio CorrectUmis -i ${f} -o ${f}.fixedumi.bam \
  --max-mismatches=3 \
  --min-distance=1 \
  -M ${f}.metrics.txt \
  -r ${f}.rejected.bam -t RX -u GAGACGAT TTCCAAGG CGCATGAT ACGGAACA CGGCTAAT GCTATCCT TGGACTCT \
  ATCCAGAG CTTAGGAC GTGCCATA TCGCTGTT TTCGTTGG AAGCACTG GTCGAAGA ACCACGAT GATTACCG GCACAACT \
  GCGTCATT GAAGGAAG ACTGAGGT TGAAGACG GTTACGCA AGCGTGTT GATCGAGT TTGCGAAG CTGTTGAC GATGTGTG \
  ACGTTCAG TTGCAGAC CAATGTGG ACGACTTG ACTAGGAG

done

```
