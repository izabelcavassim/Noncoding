# Noncoding
DFE of noncoding regions 


# Data
30X 1000 genomes


# Pruning YRU population (01_subsampling_new.sh)
```bash
#!/bin/bash
for i in /u/project/klohmuel/DataRepository/Human/Variants/VCF/*.vcf.gz; do 
	filename="$(basename -- $i)"; 
	filename2="Yoruba"$filename; 
	echo $filename2;
	echo $filename2;
	qsub -v i=$i -v filename2=$filename2 01_subsampling_submission_new.sh; 
done 
```
``` bash
#!/bin/bash
#$ -l h_rt=23:59:59,h_data=16G
#$ -pe shared 2
#$ -cwd
#$ -V
#$ -o log.$JOB_ID.out
#$ -e log.$JOB_ID.err
#$ -N filter
/u/local/Modules/default/init/modules.sh
. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
conda activate noncoding
i=$i
echo $i
filename2=$filename2
echo $filename2
bcftools view -S /u/home/m/mica20/project-kirk-bigdata/noncoding_project/data/yoruba_names.txt $i --force-samples | bgzip -c > '/u/scratch/m/mica20/'$filename2
```

# Including SNPs only
``` bash
#!/bin/bash
#qsub -cwd -V -N yYRIsubsetVCFSNPs -l highp,time=60:00:00,h_data=4G -t 1-22:1 -M mica20 -m bea 02_SNPs_only.sh  
. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
conda activate noncoding

out_dir='/u/scratch/m/mica20/SNPs_only' #scratch directory
in_dir='/u/scratch/m/mica20' #scratch directory

CHROM=${SGE_TASK_ID}
OUT_PREFIX='YRI_subset_n108_SNPs'

vcf_in=${in_dir}/'YorubaCCDG_13607_B01_GRM_WGS_2019-02-19_chr'${CHROM}'.recalibrated_variants.vcf.gz'
vcf_out=${out_dir}/${OUT_PREFIX}'_chr'$CHROM'.recalibrated_variants.vcf.gz'

vcftools --gzvcf ${vcf_in} --remove-indels --recode --recode-INFO-all --out ${vcf_out}
bgzip ${vcf_out}
```

# SNP annotation
``` bash
#!/bin/bash
#qsub -cwd -V -N SNP_annotation -l highp,time=60:00:00,h_data=16G -t 1-22:1 -M mica20 -m bea 03_SNP_annotation_test.sh
   
/u/local/Modules/default/init/modules.sh
. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
conda activate noncoding

in_dir='/u/scratch/m/mica20/SNPs_only'
out_dir='/u/scratch/m/mica20/SNPs_only'
#java -Xmx8g -jar /com/extra/SnpEff/LATEST/jar-bin/snpEff.jar -c /com/extra/SnpEff/LATEST/etc/snpEff.config -geneId GRCh38.82 ../../merge/merged_humans_chimp_GATK_conv.vcf > merged_humans_chimp_GATK_res.vcf

# this is necessary to be able to connect to Host on the cluster to download the annotation database
#export _JAVA_OPTIONS="-Dhttp.proxyHost=in -Dhttp.proxyPort=3128 -Dhttps.proxyHost=in -Dhttps.proxyPort=3128"
java -Djava.net.useSystemProxies=true ...

#check which databases are available
#java -jar /u/home/m/mica20/project-kirk-bigdata/noncoding_project/software/snpEff/snpEff.jar databases

#download an annotation database
#java -jar /u/home/m/mica20/project-kirk-bigdata/noncoding_project/software/snpEff/snpEff.jar download -v GRCh38.99

# loop through chromosomes
CHROM=${SGE_TASK_ID}
OUT_PREFIX='YRI_subset_n108_SNPs'

vcf_in=${in_dir}/'YRI_subset_n108_SNPs_chr'${CHROM}'.recalibrated_variants.vcf.gz.recode.vcf' 
vcf_out=${out_dir}/${OUT_PREFIX}'_chr'${CHROM}'.recalibrated_variants.SNPeff.vcf'
stats_out=${out_dir}/'snpEffs_stats_chr'${CHROM}'.txt'

# run the code!
java -Xmx8g -jar /u/home/m/mica20/project-kirk-bigdata/noncoding_project/software/snpEff/snpEff.jar -v GRCh38.99 $vcf_in -stats $stats_out > $vcf_out

# zip the files
bgzip ${vcf_out}
```

# Filtering low quality SNPs



