K-mer GWAS pipeline
=====================

These scripts show how to perform k-mer GWAS on maize population in our study.

1. Get k-mer sets for all accessions
=======================================
/kmerGWAS/run_kmersGWAS.sh
```Bash
#! /bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
conda activate david
library_dir="/data05/bxin/softwares/kmerGWAS/" # The path to the directory with the library
base_dir="/data05/bxin/kmerGWAS" # This will be our working directory
THREADS=16 #Max number of threads to use 
KMC="$library_dir/external_programs/kmc_v3"
# Parameters for building the k-mers table
K=31 #k-mer size
THRESHOLD_COUNT=2 #threshold to count k-mers per sample
# Sub-folders to use
samples_dir="$base_dir/samples"
# Create directories
mkdir -p $samples_dir
# Go over all strain names in the phenotype file and create separate k-mers list
line_index=0
while read name value 
do
    if [ $line_index -gt 0 ]
    then
        echo "Working on sample #$line_index : $name"
        sample_dir="$samples_dir/$name"
        mkdir -p $sample_dir # Create each samples directory
        if [ -e $sample_dir/kmers_with_strand ]; then   # check if the kmers_with_strand file exist
            echo "kmers_with_strand exists!"
            continue
        fi
        # 1. Download sequence files
        bam=obs://mem-result/$name/$name\_mem_rd.sort.bam
        prefix=$name\_mem_rd
        cd /evs/ &&/bin/rm -rf /evs/* &&  mkdir $name && cd $name 
        date
        obsutil cp $bam /evs/$name/$prefix.sort.bam  -f  -u  -j=3  -vmd5 -vlength >  obs_bam_download.log  2> speed1  
        obsutil cp $bam /evs/$name/$prefix.sort.bam  -f  -u  -j=3  -vmd5 -vlength >>  obs_bam_download.log 2>> speed1 
        obsutil cp $bam /evs/$name/$prefix.sort.bam  -f  -u  -j=3  -vmd5 -vlength >>  obs_bam_download.log 2>> speed1 
        obsutil cp $bam /evs/$name/$prefix.sort.bam  -f  -u  -j=3  -vmd5 -vlength >>  obs_bam_download.log 
        obsutil cp $bam /evs/$name/$prefix.sort.bam  -f  -u  -j=3  -vmd5 -vlength >>  obs_bam_download.log 
        date
        date
        samtools fastq -@ $THREADS $prefix.sort.bam \
            -1 $prefix\_1.fq \
            -2 $prefix\_2.fq \
            -0 /dev/null -s /dev/null -n  && echo 'fastq finished'
        date
        /bin/rm $prefix.sort.bam
        echo /evs/$name/$prefix\_1.fq > input_files.txt
        echo /evs/$name/$prefix\_2.fq >> input_files.txt   # File with list of fastq files
        # 2. Count k-mers twice using KMC
        export LD_LIBRARY_PATH=/data05/bxin/softwares/kmerGWAS:$LD_LIBRARY_PATH
        # 2.1. KMC with canonization    
        CMD="$KMC -t$THREADS -k$K -ci$THRESHOLD_COUNT @input_files.txt output_kmc_canon ./ 1> kmc_canon.1"
        echo "$CMD"; eval $CMD && echo 'KMC with canonization finished'
        # 2.2. KMC without canonization
        CMD="$KMC -t$THREADS -k$K -ci0 -b @input_files.txt output_kmc_all ./ 1> kmc_all.1"
        echo "$CMD"; eval $CMD && echo 'KMC without canonization finished'
        # 3. create a list of k-mers from the two KMC DBs
        export LD_LIBRARY_PATH=/data05/bin/software/gcc-7.5.0/gcc-build-7.5/bin:$LD_LIBRARY_PATH
        export PATH=/data05/bin/software/gcc-7.5.0/gcc-build-7.5/bin:$PATH
        CMD="$library_dir/bin/kmers_add_strand_information -c output_kmc_canon -n output_kmc_all -k $K -o kmers_with_strand"
        echo "$CMD"; eval $CMD && echo 'kmer list finished'
        # 4. copy result to the folder of the cluster
        cp kmers_with_strand $sample_dir && echo 'copy kmer list finished'
        cp obs_bam_download.log $sample_dir && echo 'copy log file finished'
    fi
    ((line_index++))
done < "./sample_list"
 
/data05/bxin/kmerGWAS/submit_jobs_example.sh
 
#!/bin/sh
####### step 1: create sample list, 13 samples lists are created in total ##############
obsutil ls obs://mem-result/ -limit 9999 | grep _mem_rd.sort.bam | awk -F"/" '{print $4}' > all_samples
LEN=`wc -l all_samples | awk '{print $1}'`
echo $LEN
offset=$((LEN/60))    # each sample list has 300 files
leftovers=$((LEN%60))
for i in $(seq 1 $offset)
do
    # file_list=sample_list.$i
    # from=$(((i-1)*60+1))
    # end=$((i*60))
    # echo -e "accession_id\tphenotype_value" > $file_list
    # sed -n ${from},${end}p all_samples | awk '{print $1"\t1"}' >> $file_list
    # cp run_kmersGWAS.sh run_kmersGWAS.$i.sh
    # sed -i s/sample_list/sample_list.$i/g run_kmersGWAS.$i.sh
    # submit jobs
    qsub -cwd -l vf=32g,p=16 -q cn-evs600sas-6 run_kmersGWAS.$i.sh
Done
```

2. Filter the needed k-mers
=============================
/kmerGWAS/step2_combine_filter.sh
```Bash
#! /bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
conda activate david
export LD_LIBRARY_PATH=/data05/bxin/kmerGWAS/my_lib:$LD_LIBRARY_PATH
library_dir="/data05/bxin/softwares/kmerGWAS/" # The path to the directory with the library
base_dir="/data05/bxin/kmerGWAS" # This will be our working directory
operation="$library_dir/bin/list_kmers_found_in_multiple_samples"
 
# cd /data05/bxin/kmerGWAS/
# mkdir combine_and_filter
cd /data05/bxin/kmerGWAS/combine_and_filter/
# 1 create samples files list
ll /data05/bxin/kmerGWAS/samples | awk '{printf "/data05/bxin/kmerGWAS/samples/%s/kmers_with_strand\t%s\n", $NF,$NF}' > kmers_list_paths.txt
# 2 combine and filter
CMD="$operation -l kmers_list_paths.txt -k 31 --mac 5 -p 0.2 -o kmers_to_use"
echo "$CMD"; eval $CMD && echo 'combine and filter finished'
```

3. Build k-mer table for k-mers and accessions
===================================================
/kmerGWAS/step3_create_kmer_table.sh
```Bash
#! /bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
conda activate david
export LD_LIBRARY_PATH=/data05/bxin/kmerGWAS/my_lib:$LD_LIBRARY_PATH
library_dir="/data05/bxin/softwares/kmerGWAS/" # The path to the directory with the library
base_dir="/data05/bxin/kmerGWAS" # This will be our working directory
operation="$library_dir/bin/build_kmers_table"
cd /data05/bxin/kmerGWAS/combine_and_filter
 
# create kmer table
CMD="$operation -l kmers_list_paths.txt -k 31 -a kmers_to_use -o kmers_table"
echo "$CMD"; eval $CMD && echo 'create kmer table finished'
```