# grape_ref_guided_assembly

Modules:

<pre style="color: silver; background: black;">
module load fastqc/0.11.5
module load Trimmomatic/0.36
module load bcftools/1.6
module load BEDtools/2.27.1
module load bowtie2/2.2.9
module load samtools/1.7
module load GATK/3.7
module load picard/2.9.2
module load fastqc
module load Trimmomatic/0.36
module load MUMmer/4.0.2
module load seqtk/1.2
module load AMOS/3.1.0
module load bamtools/2.5.1
module load SOAP-denovo/2.04</pre>

Non-module programs:
<pre style="color: silver; background: black;">
progPath=/UCHC/PROJECTS/Vitis-genome/guided_assembly/Programs/
progRemovShortSeq=${progPath}/RemoveShortSeq.jar
progGetBlocks=${progPath}/GetBlocks.jar
progFastaToAmos=${progPath}/FastaToAmos.jar
progWriteSoapConfig=${progPath}/WriteSoapConfig.jar
progFastaStats=${progPath}/FastaStats.jar</pre>

Global variables:
<pre style="color: silver; background: black;">
workPathFiles=/UCHC/PROJECTS/Vitis-genome/guided_assembly
ref=/UCHC/PROJECTS/Vitis-genome/guided_assembly/GCF_000003745.3_12X_genomic.fna
refRed=/UCHC/PROJECTS/Vitis-genome/guided_assembly/Vitis_10kb
primerFile=/UCHC/PROJECTS/Vitis-genome/guided_assembly/Programs/AdapterSeq_new.fa
primerFileMP=/UCHC/PROJECTS/Vitis-genome/guided_assembly/Programs/AdapterSeqMP_new.fa

NThreads=32      &#35; set the number of threads of every parallelizable step
maxReadLength=100
kmer=61         &#35;define best K (need to be adapted)
</pre>

<strong>Possible problem variables:</strong>
These variables were unchanged from the source from which this template was taken. Check each step for the occurence of these variables and _how_ they're being used.

<pre style="color: silver; background: black;">
&#35; mate-pair libraries ---------------------
mateName=Aly_sim_Mate
mateLib=(3000 5000 7000 11000 15000)      &#35; set insertion libraries
mateInsLow=(2000 4000 6000 9000 13000)    &#35; lower bound of insertion size
mateInsHigh=(4000 6000 8000 13000 17000)  &#35; upper bound of insertion size
mateLibsd=(400 400 400 400 400)           &#35; sd of insertion size

&#35; list of files with forward reads according to lib array
mateReads1=(Aly_3kb_1.fq Aly_5kb_1.fq Aly_7kb_1.fq Aly_11kb_1.fq Aly_15kb_1.fq)
&#35; list of files with reverse reads according to lib array
mateReads2=(Aly_3kb_2.fq Aly_5kb_2.fq Aly_7kb_2.fq Aly_11kb_2.fq Aly_15kb_2.fq)
&#35; short names of libraries
mateShortNames=(Aly_3kb Aly_5kb Aly_7kb Aly_11kb Aly_15kb)</pre>

Concatenated trimmed reads (paired-end):
<pre style="color: silver; background: black;">
reads1=$( (ls /UCHC/PROJECTS/Vitis-genome/trimmed_seq/all_L1_trimmed.fastq) )
reads2=$( (ls /UCHC/PROJECTS/Vitis-genome/trimmed_seq/all_L2_trimmed.fastq) )</pre>

Prefix:
<pre style="color: silver; background: black;">
shortNames=(Vitis_150)
</pre>

Work path and log:
<pre style="color: silver; background: black;">
workPath=${workPathFiles}/${name}_soap
log=${workPath}/log_${name}_soap.txt
</pre>


<strong>SCRIPT:</strong>
<pre style="color: silver; background: black;">
&#35;!/bin/bash
&#35;SBATCH --job-name=guid_assembly
&#35;SBATCH -N 1
&#35;SBATCH -n 1
&#35;SBATCH -c 32
&#35;SBATCH --partition=general
&#35;SBATCH --mail-type=END
&#35;SBATCH --mem=200G
&#35;SBATCH --mail-user=pjmartinezgarcia1@gmail.com
&#35;SBATCH -o gui_ass_grape_%j.out
&#35;SBATCH -e gui_ass_grape_%j.err


&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;
&#35;  Reference-guided de novo assembly - SOAP
&#35; ====================================================
&#35; Adapted from Heidi Lischer, 2015/2016 by Pedro Martinez
&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;


module load fastqc/0.11.5
module load Trimmomatic/0.36
module load bcftools/1.6
module load BEDtools/2.27.1
module load bowtie2/2.2.9
module load samtools/1.7
module load GATK/3.7
module load picard/2.9.2 
module load fastqc
module load Trimmomatic/0.36
module load MUMmer/4.0.2
module load seqtk/1.2
module load AMOS/3.1.0
module load bamtools/2.5.1
module load SOAP-denovo/2.04

&#35; set variables &#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;
workPathFiles=/UCHC/PROJECTS/Vitis-genome/guided_assembly
ref=/UCHC/PROJECTS/Vitis-genome/guided_assembly/GCF_000003745.3_12X_genomic.fna
refRed=/UCHC/PROJECTS/Vitis-genome/guided_assembly/Vitis_10kb
primerFile=/UCHC/PROJECTS/Vitis-genome/guided_assembly/Programs/AdapterSeq_new.fa
primerFileMP=/UCHC/PROJECTS/Vitis-genome/guided_assembly/Programs/AdapterSeqMP_new.fa
progPath=/isg/shared/apps/

NThreads=32      &#35; set the number of threads of every parallelizable step
maxReadLength=100
kmer=61         &#35;define best K (need to be adapted)


&#35; paired-end libraries -------------------
name=Vitis           &#35; set name of your species
lib=(150)      &#35; set insertion libraries
insLow=(50)     &#35; lower bound of insertion size
insHigh=(300)  &#35; upper bound of insertion size
libsd=(34)       &#35; sd of insertion size

&#35; list of files with forward reads according to lib array
reads1=$( (ls /UCHC/PROJECTS/Vitis-genome/trimmed_seq/all_L1_trimmed.fastq) )
&#35; list of files with rewerse reads according to lib array
reads2=$( (ls /UCHC/PROJECTS/Vitis-genome/trimmed_seq/all_L2_trimmed.fastq) )
&#35; short names of libraries
shortNames=(Vitis_150) 


&#35; mate-pair libraries ---------------------
mateName=Aly_sim_Mate
mateLib=(3000 5000 7000 11000 15000)      &#35; set insertion libraries
mateInsLow=(2000 4000 6000 9000 13000)    &#35; lower bound of insertion size
mateInsHigh=(4000 6000 8000 13000 17000)  &#35; upper bound of insertion size
mateLibsd=(400 400 400 400 400)           &#35; sd of insertion size

&#35; list of files with forward reads according to lib array
mateReads1=(Aly_3kb_1.fq Aly_5kb_1.fq Aly_7kb_1.fq Aly_11kb_1.fq Aly_15kb_1.fq)
&#35; list of files with reverse reads according to lib array
mateReads2=(Aly_3kb_2.fq Aly_5kb_2.fq Aly_7kb_2.fq Aly_11kb_2.fq Aly_15kb_2.fq)
&#35; short names of libraries
mateShortNames=(Aly_3kb Aly_5kb Aly_7kb Aly_11kb Aly_15kb)


&#35; set work path ---------------------------
workPath=${workPathFiles}/${name}_soap
&#35; log file
log=${workPath}/log_${name}_soap.txt


&#35; Programs --------------------------------
progPath=/home/hlischer/Programs
progFastQC=${progPath}/FastQC/fastqc
progTrimmomatic=${progPath}/Trimmomatic-0.32/trimmomatic-0.32.jar
progSamtools=samtools
progVcfutils=${progPath}/bcftools-1.3/vcfutils.pl
progBcftools=${progPath}/bcftools-1.3/bcftools
progBamtools=${progPath}/bamtools/bin/bamtools-2.3.0
progBedtools=${progPath}/bedtools2-2.19.1/bin/bedtools
progPicard=${progPath}/picard-tools-1.109
progBowtie2=${progPath}/bowtie2-2.2.1/bowtie2
progSeqtk=${progPath}/seqtk-master/seqtk
progAmos=${progPath}/amos-3.1.0/bin/
progNucmer=${progPath}/mummer-4.0.0beta2/nucmer
progGatk=${progPath}/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar
progSoapdenovo2=${progPath}/SOAPdenovo2-src-r240
&#35;progGatk=${progPath}/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar
&#35;progSoapdenovo2=${progPath}/SOAPdenovo2-src-r240
&#35;progSoapdenovo2=${progPath}/SOAPdenovo2-src-r240
&#35;progPicard=${progPath}/picard-tools-1.109
&#35;progBowtie2=${progPath}/bowtie2-2.2.1/bowtie2
&#35;progSeqtk=${progPath}/seqtk-master/seqtk
&#35;progAmos=${progPath}/amos-3.1.0/bin/
&#35;progNucmer=${progPath}/mummer-4.0.0beta2/nucmer
&#35;progGatk=${progPath}/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar
&#35;progSoapdenovo2=${progPath}/SOAPdenovo2-src-r240
progRemovShortSeq=${progPath}/RemoveShortSeq.jar
progGetBlocks=${progPath}/GetBlocks.jar
progFastaToAmos=${progPath}/FastaToAmos.jar
progWriteSoapConfig=${progPath}/WriteSoapConfig.jar
progFastaStats=${progPath}/FastaStats.jar
progSplitSeqLowCov=${progPath}/SplitSeqLowCov.jar

&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;



&#35; run pipeline &#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;
mkdir ${workPath}


&#35; 1. Step: quality/adapter trimming and quality check:
&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;
  &#35; quality check ----------
  echo "quality check of raw reads..."
  echo "quality check of raw reads..." > $log
  
  cd $workPathFiles
  fastqcOut=${workPathFiles}/FastQC_het
  mkdir -p ${fastqcOut}
    
  for i in ${!lib[*]}  &#35;for all indexes in the array
  do 
    fastqc -t ${NThreads} -o ${fastqcOut} ${reads1[i]} ${reads2[i]}
  done  
  
  for i in ${!mateLib[*]}  &#35;for all indexes in the array
  do 
    fastqc -t ${NThreads} -o ${fastqcOut} ${mateReads1[i]} ${mateReads2[i]}
  done  
  
  
  &#35; quality/adapter trimming --------
  &#35; - remove Illumina adapters provided in the primer files
  &#35; - remove leading and trailing low quality basses (<3) or N
  &#35; - 4 base sliding window -> remove when average quality is < 15
  &#35; - remove reads which are shorter than 40 bp
  echo "quality/adapter trimming..."
  echo "quality/adapter trimming..." >> $log
  
  trimOut=${workPathFiles}/Trim_het
  mkdir $trimOut
  
  read1TrimPair=()
  read1TrimUnPair=()
  read2TrimPair=()
  read2TrimUnPair=()
  for i in ${!lib[*]}  &#35;for all indexes in the array
  do
    read1TrimPair[i]=${trimOut}/${shortNames[i]}_R1_trimPair.fastq 
    read1TrimUnPair[i]=${trimOut}/${shortNames[i]}_R1_trimUnPair.fastq 
    read2TrimPair[i]=${trimOut}/${shortNames[i]}_R2_trimPair.fastq 
    read2TrimUnPair[i]=${trimOut}/${shortNames[i]}_R2_trimUnPair.fastq 
    
    java -jar /isg/shared/apps/Trimmomatic/0.36/trimmomatic-0.36.jar  PE -threads ${NThreads} ${reads1[i]} ${reads2[i]} ${read1TrimPair[i]} ${read1TrimUnPair[i]} ${read2TrimPair[i]} ${read2TrimUnPair[i]} ILLUMINACLIP:${primerFile}:2:30:7:5:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 >> $log
  done
  
  mateRead1TrimPair=()
  mateRead1TrimUnPair=()
  mateRead2TrimPair=()
  mateRead2TrimUnPair=()
  for i in ${!mateLib[*]}  &#35;for all indexes in the array
  do
    mateRead1TrimPair[i]=${trimOut}/${mateShortNames[i]}_R1_trimPair.fastq 
    mateRead1TrimUnPair[i]=${trimOut}/${mateShortNames[i]}_R1_trimUnPair.fastq 
    mateRead2TrimPair[i]=${trimOut}/${mateShortNames[i]}_R2_trimPair.fastq 
    mateRead2TrimUnPair[i]=${trimOut}/${mateShortNames[i]}_R2_trimUnPair.fastq 
    
    java -jar /isg/shared/apps/Trimmomatic/0.36/trimmomatic-0.36.jar PE -threads ${NThreads} ${mateReads1[i]} ${mateReads2[i]} ${mateRead1TrimPair[i]} ${mateRead1TrimUnPair[i]} ${mateRead2TrimPair[i]} ${mateRead2TrimUnPair[i]} ILLUMINACLIP:${primerFileMP}:2:30:7:5:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 >> $log
  done
  
  
  &#35; quality check ----------
  echo "quality check of trimmed reads..."
  echo "quality check of trimmed reads..." >> $log
  
  for i in ${!lib[*]}  &#35;for all indexes in the array
  do 
    fastqc -t ${NThreads} -o ${fastqcOut} ${read1TrimPair[i]} ${read2TrimPair[i]} ${read1TrimUnPair[i]} ${read2TrimUnPair[i]}
  done
  
  for i in ${!mateLib[*]}  &#35;for all indexes in the array
  do 
    fastqc -t ${NThreads} -o ${fastqcOut} ${mateRead1TrimPair[i]} ${mateRead2TrimPair[i]} ${mateRead1TrimUnPair[i]} ${mateRead2TrimUnPair[i]}
  done
  
  

&#35; 2. Step: map reads against reference 
&#35;          and define blocks and superblocks
&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;
  &#35; prepare Reference ---------
  echo "prepare reference..."
  echo "prepare reference..." >> $log
  
  &#35;remove scaffolds shorter than 10 kb
  java -jar ${progRemovShortSeq} -i $ref -o ${refRed}.fa -length 10000 
  
  &#35;create index files
  ${progSamtools} faidx ${refRed}.fa
  java -jar /isg/shared/apps/picard/picard-tools-2.9.2/picard.jar CreateSequenceDictionary R=${refRed}.fa O=${refRed}.dict
  

  &#35; map reads against reference ----------
  echo "run reference mapping..."
  echo "run reference mapping..." >> $log
  
  cd ${workPath}
  
  &#35;merge unpaired files
  readTrimUnPair=${trimOut}/${name}_trimUnpair_mod.fastq
  cat ${read1TrimUnPair[*]} ${read2TrimUnPair[*]} > ${readTrimUnPair}
  libUnpair=Unpair
  
  &#35; index reference file  
  echo "run reference mapping..."  
  bowtie2-build ${refRed}.fa ${refRed}
  
  mappedAll=()
  unmapped=()
  mapped=()
  mappedFiltered=()
  bowtieFailPair=()
  bowtieFailUnPair1=()
  bowtieFailUnPair2=()
  count=0
  for i in ${!lib[*]}  &#35;for all indexes in the array
  do 
    mappedAll[i]=${workPath}/${shortNames[i]}_all.sorted.bam
    unmapped[i]=${workPath}/${shortNames[i]}_unmapped.sorted.bam
    mapped[i]=${workPath}/${shortNames[i]}.sorted.bam
    mappedFiltered[i]=${workPath}/${shortNames[i]}.sorted.filtered.bam
    bowtieFailPair[i]=${workPath}/${shortNames[i]}_failPair.fastq
    bowtieFailUnPair1[i]=${workPath}/${shortNames[i]}_failUnPairR1.fastq
    bowtieFailUnPair2[i]=${workPath}/${shortNames[i]}_failUnPairR2.fastq    
    (
      bowtie2 --fast-local -p 1 -q --phred33 -I ${insLow[i]} -X ${insHigh[i]} -x ${refRed} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]} -U ${read1TrimUnPair[i]},${read2TrimUnPair[i]} | samtools view -bS - | samtools sort - -T ${shortNames[i]} -o ${mappedAll[i]}
      ${progSamtools} index ${mappedAll[i]}
    
      &#35;filter unmapped reads    
      samtools view -b -F 4 ${mappedAll[i]} > ${mapped[i]}
      samtools index ${mappped[i]}
    
      &#35;get unmapped reads
      samtools view -b -f 4 ${mappedAll[i]} > ${unmapped[i]}
      samtools view -b -f 9 ${unmapped[i]} > ${unmapped[i]%.sorted.bam}_pair.sorted.bam
      java -jar /isg/shared/apps/picard/picard-tools-2.9.2/picard.jar SamToFastq INPUT=${unmapped[i]%.sorted.bam}_pair.sorted.bam FASTQ=${bowtieFailPair[i]%.fastq}.1.fastq SECOND_END_FASTQ=${bowtieFailPair[i]%.fastq}.2.fastq
      rm ${unmapped[i]%.sorted.bam}_pair.sorted.bam
      samtools view -b -F 8 -f 64 ${unmapped[i]} | bamtools convert -format fastq -out ${bowtieFailUnPair1[i]}  
      samtools view -b -F 8 -f 128 ${unmapped[i]} | bamtools convert -format fastq -out ${bowtieFailUnPair2[i]}
    
      bamtools stats -in ${mappedAll[i]} >> $log
      echo "--> ${mappedAll[i]}" >> $log
    
      &#35;filter for mapping quality >=10    
      samtools view -b -q 10 ${mapped[i]} > ${mappedFiltered[i]}
      bamtools stats -in ${mappedFiltered[i]} >> $log
      echo "--> ${mappedFiltered[i]}" >> $log
    
      &#35;check insertion size
      &#35;java -jar /isg/shared/apps/picard/picard-tools-2.9.2/picard.jar CollectInsertSizeMetrics R=${refRed}.fa I=${mapped[i]} O=${mapped[i]%.bam}_insertSize.txt HISTOGRAM_FILE=${mapped[i]%.bam}_insertSizeHist.pdf
    ) &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait
    
&#35;  mateMappedAll=()
&#35;  mateUnmapped=()
&#35;  mateMapped=()
&#35;  mateBowtieFailPair=()
&#35;  mateBowtieFailUnPair1=()
&#35;  mateBowtieFailUnPair2=()
&#35;  count=0
&#35;  for i in ${!mateLib[*]}  &#35;for all indexes in the array
&#35;  do 
&#35;    mateMappedAll[i]=${workPath}/${mateShortNames[i]}_all.sorted.bam
&#35;    mateMapped[i]=${workPath}/${mateShortNames[i]}.sorted.bam
&#35;    mateUnmapped[i]=${workPath}/${mateShortNames[i]}_unmapped.sorted.bam
&#35;    mateBowtieFailPair[i]=${workPath}/${mateShortNames[i]}_failPair.fastq
&#35;    mateBowtieFailUnPair1[i]=${workPath}/${mateShortNames[i]}_failUnPairR1.fastq
&#35;    mateBowtieFailUnPair2[i]=${workPath}/${mateShortNames[i]}_failUnPairR2.fastq      
&#35;    ( 
&#35;      ${progBowtie2} --fast-local -p 1 -q --phred33 -I ${mateInsLow[i]} -X ${mateInsHigh[i]} -x ${refRed} -1 ${mateRead1TrimPair[i]} -2 ${mateRead2TrimPair[i]} --rf | ${progSamtools} view -bS - | ${progSamtools} sort - -T ${mateShortNames[i]} -o ${mateMappedAll[i]}
&#35;      ${progSamtools} index ${mateMappedAll[i]}
&#35;    
&#35;      &#35;filter unmapped reads    
&#35;      ${progSamtools} view -b -F 4 ${mateMappedAll[i]} > ${mateMapped[i]}
&#35;      ${progSamtools} index ${mateMapped[i]}
&#35;    
&#35;      &#35;get unmapped reads
&#35;      ${progSamtools} view -b -f 4 ${mateMappedAll[i]} > ${mateUnmapped[i]}
&#35;      ${progSamtools} view -b -f 9 ${mateUnmapped[i]} > ${mateUnmapped[i]%.sorted.bam}_pair.sorted.bam
&#35;      java -jar ${progPicard}/SamToFastq.jar INPUT=${mateUnmapped[i]%.sorted.bam}_pair.sorted.bam FASTQ=${mateBowtieFailPair[i]%.fastq}.1.fastq SECOND_END_FASTQ=${mateBowtieFailPair[i]%.fastq}.2.fastq
&#35;      rm ${mateUnmapped[i]%.sorted.bam}_pair.sorted.bam
&#35;      ${progSamtools} view -b -F 8 -f 64 ${mateUnmapped[i]} | ${progBamtools} convert -format fastq -out ${mateBowtieFailUnPair1[i]} 
&#35;      ${progSamtools} view -b -F 8 -f 128 ${mateUnmapped[i]} | ${progBamtools} convert -format fastq -out ${mateBowtieFailUnPair2[i]}
&#35;    
&#35;      ${progBamtools} stats -in ${mateMappedAll[i]} >> $log
&#35;      echo "--> ${mateMappedAll[i]}" >> $log
&#35;    
&#35;      &#35;check insertion size
&#35;      &#35;java -jar ${progPicard}/CollectInsertSizeMetrics.jar R=${refRed}.fa I=${mateMapped[i]} O=${mateMapped[i]%.bam}_insertSize.txt HISTOGRAM_FILE=${mateMapped[i]%.bam}_insertSizeHist.pdf
&#35;    ) &
&#35;    let count+=1
&#35;    [[ $((count%${NThreads})) -eq 0 ]] && wait
&#35;  done
&#35;  wait  
  
  &#35;merge alignment files
  mappedMerged=${workPath}/${name}_mate.sorted_mapped
  samtools merge ${mappedMerged}.bam ${mapped[*]} ${mateMapped[*]}
  
  
  &#35; get blocks and superblocks ----------
  echo "get blocks and superblocks..."
  echo "get blocks and superblocks..." >> $log
  
  &#35;get coverage along genome
  covFile=${workPath}/${name}_mate_coverage.txt
  bedtools genomecov -ibam ${mappedMerged}.bam -bga > ${covFile}
  &#35;only proparly paired reads
  samtools view -bf 0x2 ${mappedMerged}.bam | samtools sort - -n | bedtools bamtobed -i - -bedpe | awk '$1 == $4' | cut -f 1,2,6 | sort -k 1,1 | bedtools genomecov -i - -bga -g ${refRed}.fa.fai > ${covFile%.txt}Paired.txt
  
  &#35;get blocks with a minimal coverage of 10 (paired-end) reads and create superblocks of at least 12000bp length and a minimal overlap of 300bp (max. overlap = 3*300bp)
  blocks=${workPath}/blocks.txt
  superblocks=${workPath}/superblocks.txt
  java -jar ${progGetBlocks} -i ${covFile} -paired ${covFile%.txt}Paired.txt -o ${blocks} -oSuper ${superblocks} -mCov 10 -sLength 12000 -sOverlap 300 -maxLength 100000



&#35; 3. Step: do deNovo assembly within superblocks
&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;
  echo "deNovo assembly within superblocks..."
  echo "deNovo assembly within superblocks..." >> $log
  
  cd ${workPath}
  soapRes=${workPath}/soapResults
  mkdir ${soapRes}
  mkdir ${soapRes}/contigs
  mkdir ${soapRes}/scaffolds
  
  &#35;split superbolcks.txt into ${NThreads} files to run it in parallel
  size=$(($(wc -l < ${superblocks})/$((NThreads))+1))
  split -l $size -d ${superblocks} ${superblocks%.txt}_
  
  count=0
  for (( j=0; j<$((NThreads)); j++ ))
  do
    file=${superblocks%.txt}_0${j}
    fileout=${superblocks%.txt}_0${j}_run.sh
    logout=${superblocks%.txt}_0${j}.log
    
    array=(${file//_/ })
    number=${array[${&#35;array[*]}-1]}
    if [ ${number} != "00" ]
    then
      number=`echo $number|sed 's/^0*//'`
    fi
    
    if [[ $number =~ ^[0-9]+$ ]]
    then
      blockNb=$(($number*$size+1))
    else
      blockNb=1
    fi
    
    printf "&#35;"'!'"/bin/bash\n"                 > ${fileout}
    printf "\n"                                >> ${fileout}
    printf "mkdir ${file}_temp\n"              >> ${fileout}
    printf "cd ${file}_temp\n"                 >> ${fileout}
    printf "\n" >> ${fileout}

    printf "blockNb=$blockNb\n" >> ${fileout}
    printf "start=\`date +%%s\`\n" >> ${fileout}
    printf "for block in \$(cat ${file})\n" >> ${fileout}
    printf "do\n" >> ${fileout}
    printf "  echo \$blockNb >> $logout\n" >> ${fileout}
    printf "\n" >> ${fileout}
    printf "  &#35;extract sequence names within specified region\n" >> ${fileout}
    
    seqNames=()
    seqBam=()
    subSeq1=()
    subSeq2=()
    for i in ${!lib[*]}
    do
      seqNames[i]=sequences_${shortNames[i]}
      seqBam[i]=sequences_${shortNames[i]}.bam
      subSeq1[i]=subseq_${shortNames[i]}_R1.fastq
      subSeq2[i]=subseq_${shortNames[i]}_R2.fastq
      printf "  samtools view -b ${mapped[i]} \$block | samtools sort - -T ${shortNames[i]} -no ${seqBam[i]}\n" >> ${fileout}
      printf "  bedtools bamtofastq -i ${seqBam[i]} -fq 1_${subSeq1[i]} -fq2 1_${subSeq2[i]}\n" >> ${fileout}

      &#35;extract paired reads with one pair unmapped
      printf "  samtools view  -b -f 72 ${seqBam[i]} | bamtools convert -format fastq | paste - - - - | sort -k1,1 -t \" \" | tr \"\\\t\" \"\\\n\" > 2_${subSeq1[i]}\n" >> ${fileout}
      printf "  samtools view  -f 72 ${seqBam[i]} | cut -f 1 | awk \'{print \$0\"/2\"}\' > ${seqNames[i]}_R2.txt\n" >> ${fileout}
      printf "  seqtk subseq ${bowtieFailUnPair2[i]} ${seqNames[i]}_R2.txt | paste - - - - | sort -k1,1 -t \" \" | tr \"\\\t\" \"\\\n\" > 2_${subSeq2[i]}\n" >> ${fileout}
  
      printf "  samtools view  -b -f 136 ${seqBam[i]} | bamtools convert -format fastq | paste - - - - | sort -k1,1 -t \" \" | tr \"\\\t\" \"\\\n\" > 3_${subSeq1[i]}\n" >> ${fileout}
      printf "  samtools view  -f 136 ${seqBam[i]} | cut -f 1 | awk \'{print \$0\"/1\"}\' > ${seqNames[i]}_R1.txt\n" >> ${fileout}
      printf "  seqtk subseq ${bowtieFailUnPair1[i]} ${seqNames[i]}_R1.txt | paste - - - - | sort -k1,1 -t \" \" | tr \"\\\t\" \"\\\n\" > 3_${subSeq2[i]}\n" >> ${fileout}
  
      printf "  cat 1_${subSeq1[i]} 2_${subSeq1[i]} 3_${subSeq1[i]} > ${subSeq1[i]}\n" >> ${fileout}
      printf "  cat 1_${subSeq2[i]} 2_${subSeq2[i]} 3_${subSeq2[i]} > ${subSeq2[i]}\n" >> ${fileout}
      printf "\n" >> ${fileout}
    done

    printf "  &#35;deNovo assembly (k41,k51,k61,k71,k81)------\n" >> ${fileout}
    printf "  &#35;SOAP2\n" >> ${fileout}
    for i in ${!lib[*]}
    do
      if [ $i == 0 ]
      then
        libList=${lib[i]}
        forwardReads=${subSeq1[i]}
        reverseReads=${subSeq2[i]}
      else
        libList=${libList},${lib[i]}
        forwardReads=${forwardReads},${subSeq1[i]}
        reverseReads=${reverseReads},${subSeq2[i]}
      fi      
    done
 
    &#35;write config file
    soapConf=soap.config
    printf "  java -jar ${progWriteSoapConfig} -insLength ${libList} -r1 ${forwardReads} -r2 ${reverseReads} -max ${maxReadLength} -ru 3 -rank -o ${soapConf}\n" >> ${fileout}

    for (( i=41; i<=81; i+=10 ))
    do
      printf "  SOAPdenovo-127mer all -s ${soapConf} -K ${i} -o sblock\${blockNb}_${i}_Soap -p 1 -F\n" >> ${fileout}
    done
    printf "  if [ ! -f sblock\${blockNb}_61_Soap.scafSeq ]\n" >> ${fileout}
    printf "  then\n" >> ${fileout}
    printf "    echo \"\$blockNb soap failed\" >> $logout\n" >> ${fileout}
    printf "  fi\n" >> ${fileout}
    printf "  cp *_Soap.contig ${soapRes}/contigs/.\n" >> ${fileout}
    printf "  cp *_Soap.scafSeq ${soapRes}/scaffolds/.\n" >> ${fileout}
    printf "\n" >> ${fileout}
    
    printf "  rm -rf ${file}_temp/*\n" >> ${fileout}
    printf "  ((blockNb++))\n" >> ${fileout}
    printf "done\n" >> ${fileout}
    printf "\n" >> ${fileout}   
    printf "rm -rf ${file}_temp\n" >> ${fileout}
    printf "end=\`date +%%s\`\n" >> ${fileout}
    printf "echo \$((end-start))\n" >> ${fileout}
    printf "\n" >> ${fileout}
    
    chmod +x ${fileout}
    ${fileout} &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait
  
    
  &#35; deNovo assembly of unassembled reads ----------
  echo "deNovo assembly of unassembled reads..."
  echo "deNovo assembly of unassembled reads..." >> $log
    
  unassFolder=${workPath}/Unassembled
  mkdir ${unassFolder}
  cd ${unassFolder}
  
  for i in ${!lib[*]}
    do
      if [ $i == 0 ]
      then
        libList=${lib[i]}
        forwardReads=${bowtieFailPair[i]%.fastq}.1.fastq
        reverseReads=${bowtieFailPair[i]%.fastq}.2.fastq
      else
        libList=${libList},${lib[i]}
        forwardReads=${forwardReads},${bowtieFailPair[i]%.fastq}.1.fastq
        reverseReads=${reverseReads},${bowtieFailPair[i]%.fastq}.2.fastq
      fi      
    done
 
  &#35;write config file
  soapConf=soap.config
  java -jar ${progWriteSoapConfig} -insLength ${libList} -r1 ${forwardReads} -r2 ${reverseReads} -max ${maxReadLength} -ru 3 -rank -o ${soapConf}

  for (( i=81; i>=41; i-=10 ))
  do
    SOAPdenovo-127mer all -s ${soapConf} -K ${i} -o Unass_${j}_Soap -p ${NThreads} -F
  done
  
  
  
&#35; 4. Step: get non-redundant supercontigs
&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35; 
  echo "get supercontigs..."
  echo "get supercontigs..." >> $log
  cd ${workPath}
  
  &#35;merge deNovo assembled superblocks into one FASTA file
  soapContigs=${soapRes}/soapContigs.fa 
  for (( j=41; j<=81; j+=10 ))
  do
    cat ${soapRes}/contigs/*${j}_Soap.contig > ${soapRes}/${j}_soapContigs.fa
  done   
  cat ${soapRes}/*_soapContigs.fa > ${soapContigs}
  
  &#35;merge Unassembled files
  orgUnass=${unassFolder}/Unass_Soap.contig
  cat ${unassFolder}/Unass_*_Soap.contig > $orgUnass
  &#35;remove short seq (<500)
  echo "remove seq < 500 in ${orgUnass}" >> $log
  unass500=${unassFolder}/Unassembled_Soap_500.fa
  java -jar ${progRemovShortSeq} -i ${orgUnass} -o ${unass500} -length 500  >> $log
  
  &#35;merge all files
  superblockSeq=${workPath}/deNovo_Superblocks.fa
  cat ${soapContigs} ${unass500} > ${superblockSeq}
  
  &#35;remove short seq (<200)
  echo "remove seq < 200 in ${superblockSeq}" >> $log
  superblockSeq200=${superblockSeq%.fa}_200.fa
  java -jar ${progRemovShortSeq} -i ${superblockSeq} -o ${superblockSeq200} -length 200 -u  >> $log
 
  &#35;remove redundency with AMOScmp
  amosFolder=${workPath}/AMOScmp
  mkdir ${amosFolder}
  cd ${amosFolder}
  
  &#35;assemble all assembled superblocks with AMOScmp to supercontigs (with the help of reference)
  &#35;changed parameters in AMOScmp: (casm-layout -t 1000 (maximum ignorable trim length), make-consensus -o 10 (minimum overlap base)) 
  superblockSeqAmos=${superblockSeq200%.fa}_Amos.afg 
  java -jar ${progFastaToAmos} -i ${superblockSeq200} -o ${superblockSeqAmos}
  supercontigs=Amos_supercontigs
  amosSupercontigs=${amosFolder}/${supercontigs}.fasta
  
  echo "run AMPScmp..." >> $log
  &#35;${progAmos}AMOScmp -D TGT=${superblockSeqAmos} -D REF=${refRed}.fa ${supercontigs}
  
  &#35; running AMPScmp step by step and use multithread nucmer to spead it up
  &#35;&#35; Building AMOS bank
  echo "  build AMPS bank..." >> $log
  /isg/shared/apps/AMOS/3.1.0/bin/bank-transact -c -z -b ${supercontigs}.bnk -m ${superblockSeqAmos}

  &#35;&#35; Collecting clear range sequences
  echo "  clear range sequences..." >> $log
  /isg/shared/apps/AMOS/3.1.0/bin/dumpreads ${supercontigs}.bnk > ${supercontigs}.seq

  &#35;&#35; Running nucmer
  echo "  run nucmer..." >> $log
  /isg/shared/apps/MUMmer/4.0.2/bin/nucmer --maxmatch --threads=${NThreads} --prefix=${supercontigs} ${refRed}.fa ${supercontigs}.seq

  &#35;&#35; Running layout
  echo "  run layout..." >> $log
  /isg/shared/apps/AMOS/3.1.0/bin/casm-layout -t 1000 -U ${supercontigs}.layout -C ${supercontigs}.conflict -b ${supercontigs}.bnk ${supercontigs}.delta

  &#35;&#35; Running consensus
  echo "  run consensus..." >> $log
  /isg/shared/apps/AMOS/3.1.0/bin/make-consensus -o 10 -B -b ${supercontigs}.bnk

  &#35;&#35; Outputting contigs
  echo "  output contigs..." >> $log
  /isg/shared/apps/AMOS/3.1.0/bin/bank2contig ${supercontigs}.bnk > ${supercontigs}.contig

  &#35;&#35; Outputting fasta
  echo "  output fasta..." >> $log
  /isg/shared/apps/AMOS/3.1.0/bin/bank2fasta -b ${supercontigs}.bnk > ${supercontigs}.fasta
      
 
 
&#35; 5. Step: map reads on supercontigs
&#35;          and de novo assemble unmapped reads
&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35; 
  echo "map reads on supercontigs and correct them..."
  echo "map reads on supercontigs and correct them..." >> $log
  
  &#35;make seqnames unique
  amosSupercontigsUnique=${amosSupercontigs%.fasta}_unique.fa
  java -jar ${progRemovShortSeq} -i ${amosSupercontigs} -o ${amosSupercontigsUnique} -length 1 -u
   
  &#35;get statistics
  echo ${amosSupercontigsUnique} >> $log
  java -jar ${progFastaStats} -i ${amosSupercontigsUnique} -min 200 >> $log
  
  &#35;prepare reference
  bowtie2-build ${amosSupercontigsUnique} ${amosSupercontigsUnique%.fa}
  
  supercontMappedAll=()
  supercontUnmapped=()
  supercontFailPair=()
  supercontMappedFiltered=()
  count=0
  for i in ${!lib[*]}  &#35;for all indexes in the array
  do 
    supercontMappedAll[i]=${amosFolder}/${shortNames[i]}_all.sorted.bam
    supercontUnmapped[i]=${amosFolder}/${shortNames[i]}_unmapped.sorted.bam
    supercontFailPair[i]=${amosFolder}/${shortNames[i]}_failPair.fastq
    supercontMappedFiltered[i]=${amosFolder}/${shortNames[i]}.filtered.sorted.bam
    (
      bowtie2 --sensitive -p 1 -q --phred33 -I ${insLow[i]} -X ${insHigh[i]} -x ${amosSupercontigsUnique%.fa} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]} | samtools view -bS - | samtools sort - -T ${shortNames[i]} -o ${supercontMappedAll[i]}
      samtools index ${supercontMappedAll[i]}
    
      &#35;get unmapped reads
      samtools view -b -f 4 ${supercontMappedAll[i]} > ${supercontUnmapped[i]}
      samtools view -b -f 9 ${supercontUnmapped[i]} > ${supercontUnmapped[i]%.sorted.bam}_pair.sorted.bam
      java -jar /isg/shared/apps/picard/picard-tools-2.9.2/picard.jar SamToFastq  INPUT=${supercontUnmapped[i]%.sorted.bam}_pair.sorted.bam FASTQ=${supercontFailPair[i]%.fastq}.1.fastq SECOND_END_FASTQ=${supercontFailPair[i]%.fastq}.2.fastq
      rm ${supercontUnmapped[i]%.sorted.bam}_pair.sorted.bam
      
      bamtools stats -in ${supercontMappedAll[i]} >> $log
      echo "--> ${supercontMappedAll[i]}" >> $log
     
      &#35;filter for mapping quality >=10    
      samtools view -b -F 4 -q 10 ${supercontMappedAll[i]} > ${supercontMappedFiltered[i]}
      bamtools stats -in ${supercontMappedFiltered[i]} >> $log
      echo "--> ${supercontMappedFiltered[i]}" >> $log
    ) &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait
       

  &#35; deNovo assemble unassembled reads ----------
  echo "deNovo assemble unassembled reads..."
  echo "deNovo assemble unassembled reads..." >> $log
  
  supercontFailUnpairMerged=${amosFolder}/${name}_failUnp.fastq
  cat ${supercontFailUnpair[*]} > ${supercontFailUnpairMerged}
  
  supercontUnassFolder=${amosFolder}/Unassembled
  mkdir ${supercontUnassFolder}
  cd ${supercontUnassFolder}
  
  for i in ${!lib[*]}
    do
      if [ $i == 0 ]
      then
        libList=${lib[i]}
        forwardReads=${supercontFailPair[i]%.fastq}.1.fastq
        reverseReads=${supercontFailPair[i]%.fastq}.2.fastq
      else
        libList=${libList},${lib[i]}
        forwardReads=${forwardReads},${supercontFailPair[i]%.fastq}.1.fastq
        reverseReads=${reverseReads},${supercontFailPair[i]%.fastq}.2.fastq
      fi      
    done
 
  &#35;write config file
  soapConf=soap.config
  java -jar ${progWriteSoapConfig} -insLength ${libList} -r1 ${forwardReads} -r2 ${reverseReads}  -max ${maxReadLength} -ru 3 -rank -o ${soapConf}

  SOAPdenovo-127mer all -s ${soapConf} -K ${kmer} -o Unass_${kmer}_Soap -p ${NThreads} -F
  cd ${supercontUnassFolder}
    
  
  &#35;remove contigs shorter than 200 bp
  supercontSeqUnass=${supercontUnassFolder}/Unass-contigs_200.fa
  java -jar ${progRemovShortSeq} -i Unass_${kmer}_Soap.contig -o ${supercontSeqUnass} -length 100 >> $log
  echo ${supercontSeqUnass} >> $log
  java -jar ${progFastaStats} -i ${supercontSeqUnass} -min 200 >> $log
  
 
 
&#35; 6. Step: map reads to all supercontics and correct them 
&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35; 
  echo "merge contigs..." 
  echo "merge contigs..." >> $log
  
  mergedFolder=${workPath}/merged_corr
  mkdir $mergedFolder
  cd $mergedFolder
  
  merged=${mergedFolder}/${name}_supercontSeq_Unass.fa
  cat ${amosSupercontigsUnique} ${supercontSeqUnass} > $merged
  
  echo "${merged} >> $log
  java -jar ${progFastaStats} -i ${merged} -min 200 >> $log
  
    bowtie2-build ${merged} ${merged%.fa}
  
  mergedMappedAll=()
  mergedMappedFiltered=()
  count=0
  for i in ${!lib[*]}  &#35;for all indexes in the array
  do 
    mergedMappedAll[i]=${mergedFolder}/${shortNames[i]}_all.sorted.bam
    mergedMappedFiltered[i]=${mergedFolder}/${shortNames[i]}.filtered.sorted.bam
    (
      bowtie2 --sensitive -p 1 -q --phred33 -I ${insLow[i]} -X ${insHigh[i]} -x ${merged%.fa} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]} | samtools view -bS - | samtools sort - -T ${shortNames[i]} -o ${mergedMappedAll[i]}
      samtools index ${mergedMappedAll[i]}
        
      bamtools stats -in ${mergedMappedAll[i]} >> $log
      echo "--> ${mergedMappedAll[i]}" >> $log
     
      &#35;filter for mapping quality >=10
      samtools view -b -F 4 -q 10 ${mergedMappedAll[i]} > ${mergedMappedFiltered[i]}
    
      bamtools stats -in ${mergedMappedFiltered[i]} >> $log
      echo "--> ${mergedMappedFiltered[i]}" >> $log
    ) &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait
  
  
  &#35; error correction ----------  
  &#35; add RG header of second file (lost while merging)
  for i in ${!lib[*]}  &#35;for all indexes in the array
  do
    echo -e "@RG\tID:${shortNames[i]}.filtered.sorted\tPL:illumina\tPU:${lib[i]}\tLB:${lib[i]}\tSM:${shortNames[i]}" >> rg
  done
  samtools view -H ${mergedMappedFiltered[1]} | cat - rg > header
      
  &#35;realign reads
  mergedMappedMerged=${mergedFolder}/${name}.filtered_RG.sorted.bam
  samtools merge -r -h header ${mergedMappedMerged} ${mergedMappedFiltered[*]}
  rm rg header
 
  samtools index ${mergedMappedMerged}
  samtools faidx ${merged}
  java -jar /isg/shared/apps/picard/picard-tools-2.9.2/picard.jar CreateSequenceDictionary R=${merged} O=${merged%.fa}.dict
  java -jar /isg/shared/apps/GATK/3.7/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${merged} -I ${mergedMappedMerged} -o target_intervals.list
  mergedMappedMergedReal=${mergedMappedMerged%.bam}_realigned.bam
  java -jar /isg/shared/apps/GATK/3.7/GenomeAnalysisTK.jar -T IndelRealigner -R ${merged} -I ${mergedMappedMerged} -targetIntervals target_intervals.list -o ${mergedMappedMergedReal}

  &#35;get alternative seq
  mergedCorr=${merged%.fa}_corr.fq
  samtools mpileup -uf ${merged} ${mergedMappedMergedReal} | bcftools call -c - | /isg/shared/apps/bcftools/1.6/bin/vcfutils.pl vcf2fq -d 1 > ${mergedCorr}    
  
  &#35;remove start and end N
  mergedCorrWN=${mergedCorr%.fq}WN.fa 
  echo ${mergedCorrWN} >> $log
  java -jar ${progRemovShortSeq} -i ${mergedCorr} -o ${mergedCorrWN} -length 100 -n -fq >> $log
  
  &#35;get statistics
  echo ${mergedCorrWN} >> $log
  java -jar ${progFastaStats} -i ${mergedCorrWN} -min 200 >> $log
  

  &#35; split sequences at places with no coverage ----------
  bowtie2-build ${mergedCorrWN} ${mergedCorrWN%.fa}
  
  mergedCorrMappedAll=()
  mergedCorrMappedFiltered=()
  count=0
  for i in ${!lib[*]}  &#35;for all indexes in the array
  do 
    mergedCorrMappedAll[i]=${mergedFolder}/${shortNames[i]}_corrWN_all.sorted.bam
    mergedCorrMappedFiltered[i]=${mergedFolder}/${shortNames[i]}_corrWN.filtered.sorted.bam
    (
      bowtie2 --sensitive -p 1 -q --phred33 -I ${insLow[i]} -X ${insHigh[i]} -x ${mergedCorrWN%.fa} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]} | samtools view -bS - | samtools sort - -T ${shortNames[i]} -o ${mergedCorrMappedAll[i]}
      samtools index ${mergedCorrMappedAll[i]}
   
      bamtools stats -in ${mergedCorrMappedAll[i]} >> $log
      echo "--> ${mergedCorrMappedAll[i]}" >> $log
    
      &#35;filter for mapping quality >=10
      samtools view -b -F 4 -q 10 ${mergedCorrMappedAll[i]} > ${mergedCorrMappedFiltered[i]}
    
      bamtools stats -in ${mergedCorrMappedFiltered[i]} >> $log
      echo "--> ${mergedCorrMappedFiltered[i]}" >> $log
    ) &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait

  mergedCorrMappedFilteredMerged=${mergedFolder}/${name}_corrWN.filtered.sorted.bam
  samtools merge ${mergedCorrMappedFilteredMerged} ${mergedCorrMappedFiltered[*]}
  bedtools genomecov -ibam ${mergedCorrMappedFilteredMerged} -bga > ${mergedCorrWN%.fa}_filteredCov.txt
  &#35;only proparly paired reads
  samtools faidx $mergedCorrWN
  samtools view -bf 0x2 ${mergedCorrMappedFilteredMerged} | samtools sort - -n | bedtools bamtobed -i - -bedpe | awk '$1 == $4' | cut -f 1,2,6 | sort -k 1,1 | bedtools genomecov -i - -bga -g ${mergedCorrWN}.fai  > ${mergedCorrWN%.fa}_filteredPairedCov.txt
  
  java -jar ${progSplitSeqLowCov} -i ${mergedCorrWN%.fa}_filteredCov.txt -paired ${mergedCorrWN%.fa}_filteredPairedCov.txt -o ${mergedCorrWN%.fa}_filteredNotCov.txt -mCov 1 -fasta ${mergedCorrWN} -fastaOut ${mergedCorrWN%.fa}_splitFiltered.fa >> $log
  echo ${mergedCorrWN%.fa}_splitFiltered.fa >> $log
  java -jar ${progFastaStats} -i ${mergedCorrWN%.fa}_splitFiltered.fa -min 200 >> $log
  


&#35; 7. Step: scaffolding and gap closing
&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;&#35;
  echo "scaffolding..." 
  echo "scaffolding..." >> $log
  
  scafFolder=${mergedFolder}/scaffold_gapClosed
  mkdir $scafFolder
  cd $scafFolder
  
  for i in ${!lib[*]}
  do
    if [ $i == 0 ]
    then
      libList=${lib[i]}
      forwardReads=${read1TrimPair[i]}
      reverseReads=${read2TrimPair[i]}
    else
      libList=${libList},${lib[i]}
      forwardReads=${forwardReads},${read1TrimPair[i]}
      reverseReads=${reverseReads},${read2TrimPair[i]}
    fi      
      mateLibList=${mateLibList},${mateLib[i]}
      mateForwardReads=${mateForwardReads},${mateRead1TrimPair[i]}
      mateReverseReads=${mateReverseReads},${mateRead2TrimPair[i]}
    fi      
  done
  
  &#35;write config file
  soapConf=${scafFolder}/soap.config
  java -jar ${progWriteSoapConfig} -insLength ${libList} -r1 ${forwardReads} -r2 ${reverseReads} -max ${maxReadLength} -ru 2 -mateInsLength ${mateLibList} -mateR1 ${mateForwardReads} -mateR2 ${mateReverseReads} -mateRu 2 -rank -o ${soapConf}
  scafFile=${name}_${kmer}
  SOAPdenovo-127mer -D -c ${mergedCorrWN%.fa}_splitFiltered.fa -K ${kmer} -g ${scafFile} -p ${NThreads}
  SOAPdenovo-127mer map -s ${soapConf} -g ${scafFile} -p ${NThreads}
  SOAPdenovo-127mer scaff -g ${scafFile} -p ${NThreads} -F
   
   
  &#35;remove scaffolds < 200 bp ----------
  scafSeq=${scafFolder}/${name}_scafSeq.fa
  echo ${scafFile}.scafSeq >> $log
  java -jar ${progRemovShortSeq} -i ${scafFile}.scafSeq -o ${scafSeq} -length 200 >> $log
  java -jar ${progRemovShortSeq} -i ${scafSeq} -o ${scafSeq%.fa}_500.fa -length 500 >> $log
  java -jar ${progRemovShortSeq} -i ${scafSeq} -o ${scafSeq%.fa}_1000.fa -length 1000 >> $log
 
  &#35;get statistics
  echo ${scafSeq} >> $log
  java -jar ${progFastaStats} -i ${scafSeq} -min 200 >> $log
  java -jar ${progFastaStats} -i ${scafSeq} -min 500 >> $log
  java -jar ${progFastaStats} -i ${scafSeq} -min 1000 >> $log

  &#35;map reads against scaffolds
  bowtie2-build ${scafSeq} ${scafSeq%.fa}
   
  scafMappedAll=()
  scafMappedFiltered=()
  count=0
  for i in ${!lib[*]}  &#35;for all indexes in the array
  do 
    scafMappedAll[i]=${scafFolder}/${shortNames[i]}_all.sorted.bam
    scafMappedFiltered[i]=${scafFolder}/${shortNames[i]}.filtered.sorted.bam
    (
      bowtie2 --sensitive -p 1 -q --phred33 -I ${insLow[i]} -X ${insHigh[i]} -x ${scafSeq%.fa} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]} | samtools view -bS - | samtools sort - -T ${shortNames[i]} -o ${scafMappedAll[i]}
      samtools index ${scafMappedAll[i]}
    
      bamtools stats -in ${scafMappedAll[i]} >> $log
      echo "--> ${scafMappedAll[i]}" >> $log
    
      &#35;filter for mapping quality >=10    
      samtools view -b -F 4 -q 10 ${scafMappedAll[i]} > ${scafMappedFiltered[i]}
    
      bamtools stats -in ${scafMappedFiltered[i]} >> $log
      echo "--> ${scafMappedFiltered[i]}" >> $log
    ) &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait</strong>
