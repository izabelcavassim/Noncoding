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


# SNP annotation


# Filtering low quality SNPs



