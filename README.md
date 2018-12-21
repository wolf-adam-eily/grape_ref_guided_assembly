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
progPath=/UCHC/PROJECTS/Vitis-genome/guided_assembly/Programs/
progRemovShortSeq=${progPath}/RemoveShortSeq.jar
progGetBlocks=${progPath}/GetBlocks.jar
progFastaToAmos=${progPath}/FastaToAmos.jar
progWriteSoapConfig=${progPath}/WriteSoapConfig.jar
progFastaStats=${progPath}/FastaStats.jar

Global variables:
<pre style="color: silver; background: black;">
workPathFiles=/UCHC/PROJECTS/Vitis-genome/guided_assembly
ref=/UCHC/PROJECTS/Vitis-genome/guided_assembly/GCF_000003745.3_12X_genomic.fna
refRed=/UCHC/PROJECTS/Vitis-genome/guided_assembly/Vitis_10kb
primerFile=/UCHC/PROJECTS/Vitis-genome/guided_assembly/Programs/AdapterSeq_new.fa
primerFileMP=/UCHC/PROJECTS/Vitis-genome/guided_assembly/Programs/AdapterSeqMP_new.fa

NThreads=32      # set the number of threads of every parallelizable step
maxReadLength=100
kmer=61         #define best K (need to be adapted)
</pre>

<strong>Possible problem variables:</strong>
These variables were unchanged from the source from which this template was taken. Check each step for the occurence of these variables and _how_ they're being used.

<pre style="color: silver; background: black;">
# mate-pair libraries ---------------------
mateName=Aly_sim_Mate
mateLib=(3000 5000 7000 11000 15000)      # set insertion libraries
mateInsLow=(2000 4000 6000 9000 13000)    # lower bound of insertion size
mateInsHigh=(4000 6000 8000 13000 17000)  # upper bound of insertion size
mateLibsd=(400 400 400 400 400)           # sd of insertion size

# list of files with forward reads according to lib array
mateReads1=(Aly_3kb_1.fq Aly_5kb_1.fq Aly_7kb_1.fq Aly_11kb_1.fq Aly_15kb_1.fq)
# list of files with reverse reads according to lib array
mateReads2=(Aly_3kb_2.fq Aly_5kb_2.fq Aly_7kb_2.fq Aly_11kb_2.fq Aly_15kb_2.fq)
# short names of libraries
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

<h2 id="First_Point_Header">Step 1 (trimming reads) (COMPLETED)</h2>
 <pre style="color: silver; background: black;">
  echo "quality check of raw reads..."
  echo "quality check of raw reads..." > $log

  cd $workPathFiles
  fastqcOut=${workPathFiles}/FastQC_het
  mkdir -p ${fastqcOut}

  for i in ${!lib[*]}  #for all indexes in the array
  do
    fastqc -t ${NThreads} -o ${fastqcOut} ${reads1[i]} ${reads2[i]}
  done

  for i in ${!mateLib[*]}  #for all indexes in the array
  do
    fastqc -t ${NThreads} -o ${fastqcOut} ${mateReads1[i]} ${mateReads2[i]}
  done

  # quality/adapter trimming --------
  # - remove Illumina adapters provided in the primer files
  # - remove leading and trailing low quality basses (<3) or N
  # - 4 base sliding window -> remove when average quality is < 15
  # - remove reads which are shorter than 40 bp
  echo "quality/adapter trimming..."
  echo "quality/adapter trimming..." >> $log

  trimOut=${workPathFiles}/Trim_het
  mkdir $trimOut

  read1TrimPair=()
  read1TrimUnPair=()
  read2TrimPair=()
  read2TrimUnPair=()
  for i in ${!lib[*]}  #for all indexes in the array
  do
    read1TrimPair[i]=${trimOut}/${shortNames[i]}_R1_trimPair.fastq
    read1TrimUnPair[i]=${trimOut}/${shortNames[i]}_R1_trimUnPair.fastq
    read2TrimPair[i]=${trimOut}/${shortNames[i]}_R2_trimPair.fastq
    read2TrimUnPair[i]=${trimOut}/${shortNames[i]}_R2_trimUnPair.fastq

    java -jar /isg/shared/apps/Trimmomatic/0.36/trimmomatic-0.36.jar  PE -threads ${NThreads} ${reads1[i]} ${reads2[i]} ${read1TrimPair
[i]} ${read1TrimUnPair[i]} ${read2TrimPair[i]} ${read2TrimUnPair[i]} ILLUMINACLIP:${primerFile}:2:30:7:5:true LEADING:3 TRAILING:3 SLID
INGWINDOW:4:15 MINLEN:40 >> $log
  done

  mateRead1TrimPair=()
  mateRead1TrimUnPair=()
  mateRead2TrimPair=()
  
    for i in ${!mateLib[*]}  #for all indexes in the array
  do
    mateRead1TrimPair[i]=${trimOut}/${mateShortNames[i]}_R1_trimPair.fastq
    mateRead1TrimUnPair[i]=${trimOut}/${mateShortNames[i]}_R1_trimUnPair.fastq
    mateRead2TrimPair[i]=${trimOut}/${mateShortNames[i]}_R2_trimPair.fastq
    mateRead2TrimUnPair[i]=${trimOut}/${mateShortNames[i]}_R2_trimUnPair.fastq

    java -jar /isg/shared/apps/Trimmomatic/0.36/trimmomatic-0.36.jar PE -threads ${NThreads} ${mateReads1[i]} ${mateReads2[i]} ${mat
ad1TrimPair[i]} ${mateRead1TrimUnPair[i]} ${mateRead2TrimPair[i]} ${mateRead2TrimUnPair[i]} ILLUMINACLIP:${primerFileMP}:2:30:7:5:tr
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 >> $log
  done


  # quality check ----------
  echo "quality check of trimmed reads..."
  echo "quality check of trimmed reads..." >> $log

  for i in ${!lib[*]}  #for all indexes in the array
  do
    fastqc -t ${NThreads} -o ${fastqcOut} ${read1TrimPair[i]} ${read2TrimPair[i]} ${read1TrimUnPair[i]} ${read2TrimUnPair[i]}
  done

  for i in ${!mateLib[*]}  #for all indexes in the array
  do
    fastqc -t ${NThreads} -o ${fastqcOut} ${mateRead1TrimPair[i]} ${mateRead2TrimPair[i]} ${mateRead1TrimUnPair[i]} ${mateRead2TrimU
ir[i]}
  done
</pre>

<h2 id="Second_Point_Header">Step 2 (aligning reads to reference) (BEING CHECKED FOR ERRORS)</h2>
<pre style="color: silver; background: black;">
#          and define blocks and superblocks
#######################################################
  # prepare Reference ---------
  echo "prepare reference..."
  echo "prepare reference..." >> $log

  #remove scaffolds shorter than 10 kb
  java -jar ${progRemovShortSeq} -i $ref -o ${refRed}.fa -length 10000

  #create index files
  ${progSamtools} faidx ${refRed}.fa
  java -jar /isg/shared/apps/picard/picard-tools-2.9.2/picard.jar CreateSequenceDictionary R=${refRed}.fa O=${refRed}.dict


  # map reads against reference ----------
  echo "run reference mapping..."
  echo "run reference mapping..." >> $log

  cd ${workPath}

  #merge unpaired files
  readTrimUnPair=${trimOut}/${name}_trimUnpair_mod.fastq
  cat ${read1TrimUnPair[*]} ${read2TrimUnPair[*]} > ${readTrimUnPair}
  libUnpair=Unpair

  # index reference file
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
  for i in ${!lib[*]}  #for all indexes in the array
  do
    mappedAll[i]=${workPath}/${shortNames[i]}_all.sorted.bam
    unmapped[i]=${workPath}/${shortNames[i]}_unmapped.sorted.bam
    mapped[i]=${workPath}/${shortNames[i]}.sorted.bam
    mappedFiltered[i]=${workPath}/${shortNames[i]}.sorted.filtered.bam
    bowtieFailPair[i]=${workPath}/${shortNames[i]}_failPair.fastq
    bowtieFailUnPair1[i]=${workPath}/${shortNames[i]}_failUnPairR1.fastq
    bowtieFailUnPair2[i]=${workPath}/${shortNames[i]}_failUnPairR2.fastq
    (
      bowtie2 --fast-local -p 1 -q --phred33 -I ${insLow[i]} -X ${insHigh[i]} -x ${refRed} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]
} -U ${read1TrimUnPair[i]},${read2TrimUnPair[i]} | samtools view -bS - | samtools sort - -T ${shortNames[i]} -o ${mappedAll[i]}
      ${progSamtools} index ${mappedAll[i]}

      #filter unmapped reads
      samtools view -b -F 4 ${mappedAll[i]} > ${mapped[i]}
      samtools index ${mappped[i]}

      #get unmapped reads
      samtools view -b -f 4 ${mappedAll[i]} > ${unmapped[i]}
      samtools view -b -f 9 ${unmapped[i]} > ${unmapped[i]%.sorted.bam}_pair.sorted.bam
      java -jar /isg/shared/apps/picard/picard-tools-2.9.2/picard.jar SamToFastq INPUT=${unmapped[i]%.sorted.bam}_pair.sorted.bam FASTQ
=${bowtieFailPair[i]%.fastq}.1.fastq SECOND_END_FASTQ=${bowtieFailPair[i]%.fastq}.2.fastq
      rm ${unmapped[i]%.sorted.bam}_pair.sorted.bam
           rm ${unmapped[i]%.sorted.bam}_pair.sorted.bam
     samtools view -b -F 8 -f 64 ${unmapped[i]} | bamtools convert -format fastq -out ${bowtieFailUnPair1[i]}
     samtools view -b -F 8 -f 128 ${unmapped[i]} | bamtools convert -format fastq -out ${bowtieFailUnPair2[i]}

     bamtools stats -in ${mappedAll[i]} >> $log
     echo "--> ${mappedAll[i]}" >> $log

     #filter for mapping quality >=10
     samtools view -b -q 10 ${mapped[i]} > ${mappedFiltered[i]}
     bamtools stats -in ${mappedFiltered[i]} >> $log
     echo "--> ${mappedFiltered[i]}" >> $log

     #check insertion size
     #java -jar /isg/shared/apps/picard/picard-tools-2.9.2/picard.jar CollectInsertSizeMetrics R=${refRed}.fa I=${mapped[i]} O=${mappe
[i]%.bam}_insertSize.txt HISTOGRAM_FILE=${mapped[i]%.bam}_insertSizeHist.pdf
   ) &
   let count+=1
   [[ $((count%${NThreads})) -eq 0 ]] && wait
 done
 wait

  mateMappedAll=()
  mateUnmapped=()
  mateMapped=()
  mateBowtieFailPair=()
  mateBowtieFailUnPair1=()
  mateBowtieFailUnPair2=()
  count=0
  for i in ${!mateLib[*]}  #for all indexes in the array
  do
    mateMappedAll[i]=${workPath}/${mateShortNames[i]}_all.sorted.bam
    mateMapped[i]=${workPath}/${mateShortNames[i]}.sorted.bam
</pre>
