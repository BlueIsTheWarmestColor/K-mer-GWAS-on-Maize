K-mer GWAS pipeline
=====================

These scripts show how to perform k-mer GWAS on maize population in our study.

1.Get k-mer sets for all accessions
-------------------------------------
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

2.Filter the needed k-mers
------------------------------
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

3.Build k-mer table for k-mers and accessions
---------------------------------------------------
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

4.Conduct association study
-------------------------------
/kmerGWAS/gwas_test/gwas_*_maf0.01.sh
```Bash
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
conda activate py_27_simon
export LD_LIBRARY_PATH=/data05/bxin/kmerGWAS/my_lib:$LD_LIBRARY_PATH
library_dir="/data05/bxin/softwares/kmerGWAS/" # The path to the directory with the library
base_dir="/data05/bxin/kmerGWAS" # This will be our working directory
operation="$library_dir/bin/emma_kinship_kmers"
gemma_path="/data05/bxin/softwares/kmerGWAS/external_programs/gemma_0_96"
cd /data05/bxin/kmerGWAS/gwas_test/
 
# run the GWAS
gwas_outdir="/data05/bxin/kmerGWAS/gwas_test/KRN_maf0.01"
CMD="python2.7 $library_dir/kmers_gwas.py --pheno cob_row_num_bluped_pheno_3367_used_sorted.pheno -k 100000 --permutations 1000 --gemma_path $gemma_path --kmers_table /data05/bxin/kmerGWAS/gwas_test/kmers_table -l 31 -p 32 --maf 0.01 --mac 5 --pattern_counter --outdir $gwas_outdir "
echo "$CMD"; eval $CMD
```

5.Generate matching table for significant k-mers and accessions
-----------------------------------------------------------------
/kmerGWAS/gwas_test/*_maf0.01/from_kmer_get_loci/step1_get_kmers_in_accessions/step1_kmer_accession.sh
```Bash
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
conda activate david
export LD_LIBRARY_PATH=/data05/bxin/kmerGWAS/my_lib:$LD_LIBRARY_PATH
library_dir="/data05/bxin/softwares/kmerGWAS/" # The path to the directory with the library
base_dir="/data05/bxin/kmerGWAS" # This will be our working directory
operation="$library_dir/bin/filter_kmers"
 
 
cd /data05/bxin/kmerGWAS/gwas_test/EH_PH_maf0.01/from_kmer_get_loci/step1_get_kmers_in_accessions/
outdir="/data05/bxin/kmerGWAS/gwas_test/EH_PH_maf0.01/from_kmer_get_loci/step1_get_kmers_in_accessions/"
 
cut -f 2 ../../kmers/pass_threshold_5per > ./threshold_5per_kmer_file.txt
sed -i '1d' threshold_5per_kmer_file.txt
sed -i 's/_.*//g' threshold_5per_kmer_file.txt
 
cut -f 2 ../../kmers/pass_threshold_10per > ./threshold_10per_kmer_file.txt
sed -i '1d' threshold_10per_kmer_file.txt
sed -i 's/_.*//g' threshold_10per_kmer_file.txt
 
CMD="$operation -t /data05/bxin/kmerGWAS/kmer_table_new_version_tool/kmers_table -k threshold_5per_kmer_file.txt -o $outdir/EH_PH_5per_kmer_in_accessions"
echo "$CMD"; eval $CMD 
 
CMD="$operation -t /data05/bxin/kmerGWAS/kmer_table_new_version_tool/kmers_table -k threshold_10per_kmer_file.txt -o $outdir/EH_PH_10per_kmer_in_accessions"
echo "$CMD"; eval $CMD 
```

6.Trace back the reads k-mer originated
--------------------------------------------
/kmerGWAS/gwas_test/step2_get_reads_redo.sh
```Bash
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
conda activate david
export LD_LIBRARY_PATH=/data05/bxin/kmerGWAS/my_lib:$LD_LIBRARY_PATH
library_dir="/data05/bxin/softwares/kmerGWAS/" # The path to the directory with the library
base_dir="/data05/bxin/kmerGWAS" # This will be our working directory
operation="/data05/bxin/softwares/fetch_reads_with_kmers/fetch_reads"
THREADS="4"
 
for name in `cat accessions_id`
do 
echo "Working on sample  : $name"
 
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
 
        # 3. create a list of k-mers from the two KMC DBs
        export LD_LIBRARY_PATH=/data05/bin/software/gcc-7.5.0/gcc-build-7.5/bin:$LD_LIBRARY_PATH
        export PATH=/data05/bin/software/gcc-7.5.0/gcc-build-7.5/bin:$PATH
        for pheno in EH_PH PH EH ASI
        do
            pheno_dir=/data05/bxin/kmerGWAS/gwas_test/$pheno\_redo_maf0.01
            echo "Working on pheno  : $pheno_dir"
            per5_kmer_fa="$pheno_dir/from_kmer_get_loci/step1_get_kmers_in_accessions/5per_kmer.fa"
            echo "$per5_kmer_fa"
            if [ -e $per5_kmer_fa ]; then   
                echo "Working on pheno/5per  : $pheno"
                CMD="mkdir $pheno_dir/from_kmer_get_loci/step2_get_reads"
                echo "$CMD"; eval $CMD
                CMD="mkdir $pheno_dir/from_kmer_get_loci/step2_get_reads/5per_kmer"
                echo "$CMD"; eval $CMD
                CMD="$operation /evs/$name/$prefix\_1.fq /evs/$name/$prefix\_2.fq $per5_kmer_fa 31 $pheno_dir/from_kmer_get_loci/step2_get_reads/5per_kmer/$name"
                echo "$CMD"; eval $CMD 
            else
                echo "$per5_kmer_fa didn't exist"
            fi
        done
 
Done
```

7.Assemble the reads significant k-mer originated
--------------------
/kmerGWAS/gwas_test/merge_reads_assembly.sh
```Bash
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
 
/data05/bxin/softwares/SPAdes-3.13.0-Linux/bin/spades.py -1 /data05/bxin/kmerGWAS/gwas_test/pheno_maf0.01/from_kmer_get_loci/step2_get_reads/5per_kmer/merge_R1.fq -2 /data05/bxin/kmerGWAS/gwas_test/pheno_maf0.01/from_kmer_get_loci/step2_get_reads/5per_kmer/merge_R2.fq -t 32 --careful -o /data05/bxin/kmerGWAS/gwas_test/pheno_maf0.01/from_kmer_get_loci/step5_merge_reads_assembly/5per_merge_redo &
wait
```

8.The correation between significant k-mers
---------------------------------------------
#### Step1 convert k-mer table to VCF format
```Bash
python /kmerGWAS/gwas_test/KRN_maf0.01/from_kmer_get_loci/step8_kmer_correlation/step1_kmer_vcf.py
```
#### Step2 calculate correlation by LD
```Bash
bash /kmerGWAS/gwas_test/KRN_maf0.01/from_kmer_get_loci/step8_kmer_correlation/step2_cal_correlation.sh
```

#### Step3 prepare file for heatmap drawing
```Bash
python /kmerGWAS/gwas_test/KRN_maf0.01/from_kmer_get_loci/step8_kmer_correlation/step3_prepare_file.py
```
#### Step4 draw correlation heatmap
```R
draw_heatmap.R
 
library(ggcorrplot)
library(ggthemes)
 
data <- read.table("DTA_input_file",header = T,row.names = 1,sep = "\t")
 
##head(data)
 
ggcorrplot(data,
         hc.order = TRUE, outline.col = "gray",tl.cex = 5,
         lab =TRUE,lab_size = 2)
ggsave("DTA.pdf",limitsize = FALSE, width = 50,height =50)
```

9.Idedntifying the variants and location of significant k-mers
----------------------------------------------------------------
#### Step1 calculate assembly scores of contigs significant k-mers realated
In this step, assembly scores are calculated for later contigs filtering
```Bash
python /kmerGWAS/gwas_test/KRN_maf0.01/from_kmer_get_loci/step5_merge_reads_assembly/5per_merge_redo/assembly_score.py
```

#### Step2 blast k-mers to contigs
Map k-mers to contigs<br>
/kmerGWAS/gwas_test/KRN_maf0.01/from_kmer_get_loci/step5_merge_reads_assembly/5per_merge_redo/blast_kmer_contig.sh
```Bash
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
makeblastdb=/data05/bxin/softwares/ncbi-blast-2.11.0+/bin/makeblastdb
blastn=/data05/bxin/softwares/ncbi-blast-2.11.0+/bin/blastn
TARGET=contigs.fasta
QUERY=../../step1_get_kmers_in_accessions/5per_kmer.fa
 
$makeblastdb  -in $TARGET -input_type fasta -dbtype nucl -title contigs_sequence -parse_seqids -out contigs_sequence
$blastn -num_threads 16 -db  contigs_sequence  -query $QUERY  -out kmer_contigs_blast -evalue 1e-2 -outfmt 6
```
#### Step3 map contigs to the reference genome
/kmerGWAS/gwas_test/KRN_maf0.01/from_kmer_get_loci/step5_merge_reads_assembly/5per_merge_redo/
```Bash
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
~/anaconda3/bin/seqkit seq -m 200 contigs.fasta > contigs_len_greater_200.fasta &
~/anaconda3/bin/seqkit seq -M 200 contigs.fasta > contigs_len_less_200.fasta &
wait
conda activate py_37_simon
 
minimap2 -t 8 -r 1500 -c --MD --eqx /data05/bxin/reference/Mo17.V2.nochr.fa contigs_len_greater_200.fasta > contigs_len_greater_200.paf_mo17v2 &
minimap2 -t 8 -x sr -c --MD --eqx /data05/bxin/reference/Mo17.V2.nochr.fa contigs_len_less_200.fasta > contigs_len_less_200.fasta.paf_mo17v2 
wait
cat contigs_len_greater_200.paf_mo17v2 contigs_len_less_200.fasta.paf_mo17v2 > contigs.paf_mo17v2
minimap2 -t 8 -c --MD --eqx /data05/bxin/reference/Zea_mays.B73_RefGen_v4.dna.toplevel.fa contigs_len_greater_200.fasta > contigs_len_greater_200.paf_B73v4 &
minimap2 -t 8 -x sr -c --MD --eqx /data05/bxin/reference/Zea_mays.B73_RefGen_v4.dna.toplevel.fa contigs_len_less_200.fasta > contigs_len_less_200.fasta.paf_B73v4 
wait
cat contigs_len_greater_200.paf_B73v4 contigs_len_less_200.fasta.paf_B73v4 >contigs.paf_B73v4
```
#### Step4 draw syntenic plot
```Bash
python /kmerGWAS/gwas_test/KRN_maf0.01/from_kmer_get_loci/step6_draw_pic/colinearity_pic/draw_pic.py
```
