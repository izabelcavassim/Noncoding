# Noncoding
DFE of noncoding regions 


# Data
30X 1000 genomes


# Subsetting only YRU population
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
bash script 01_subsampling_submission_new.sh: 
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
Annotation of SNPs was done with [snpEff](https://pcingola.github.io/SnpEff/)
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
# Syn and Non-synonymous SNPs
``` bash
#!/bin/bash
#qsub -cwd -V -N Filtering_non_syn_SNPs -l highp,time=60:00:00,h_data=16G -t 1-22:1 -M eplau -m bea 04_SNP_annotation_syn_non.sh
/u/local/Modules/default/init/modules.sh
. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
conda activate noncoding

in_dir='/u/scratch/m/mica20/SNPs_only'
out_dir_syn='/u/scratch/m/mica20/SNPs_only/Syn_snps'
out_dir_nonsyn='/u/scratch/m/mica20/SNPs_only/Non_Syn_snps'

# this is necessary to be able to connect to Host on the cluster to download the annotation database
java -Djava.net.useSystemProxies=true ...

# loop through chromosomes
CHROM=${SGE_TASK_ID}
OUT_PREFIX='YRI_subset_n108_SNPs'

vcf_in=${in_dir}/${OUT_PREFIX}'_chr'${CHROM}'.recalibrated_variants.SNPeff.vcf' 

# out vcf files
vcf_out_syn=${out_dir_syn}/${OUT_PREFIX}'_chr'${CHROM}'.recalibrated_variants.SNPeff.Syn.vcf'
vcf_out_nonsyn=${out_dir_nonsyn}/${OUT_PREFIX}'_chr'${CHROM}'.recalibrated_variants.SNPeff.Non.Syn.vcf'

# run the code: synonymous variants
java -Xmx8g -jar /u/home/m/mica20/project-kirk-bigdata/noncoding_project/software/snpEff/SnpSift.jar filter "ANN[0].EFFECT has 'synonymous_variant' " $vcf_in > $vcf_out_syn
# zip the files
bgzip ${vcf_out_syn}

# run the code: nonsynonymous variants
java -Xmx8g -jar /u/home/m/mica20/project-kirk-bigdata/noncoding_project/software/snpEff/SnpSift.jar filter "ANN[0].EFFECT has 'missense_variant' | ANN[0].EFFECT has 'stop_gained' | ANN[0].EFFECT has 'stop_lost' " $vcf_in  > $vcf_out_nonsyn
# zip the files
bgzip ${vcf_out_nonsyn}
```
# Finding coding regions
``` bash
### Finding protein coding regions of Hg38

wget -O - "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz" |\
gunzip -c | grep 'transcript_type "protein_coding"' |\
awk '($3=="exon") {printf("%s\t%s\t%s\n",$1,int($4)-1,$5);}' |\
sort -T . -k1,1 -k2,2n | bedtools merge > coding_region_hg38_gencode_v36_version_2.bed
```
# Restricting HMM state regions to noncoding regions
As explained [here](https://bedtools.readthedocs.io/en/latest/content/tools/subtract.html), one can use bedtools for that.
``` bash
for i in *states_hg38.bed; do bedtools subtract -a $i -b coding_region_hg38_gencode_v36_version_2.bed > subset_only_noncoding_$i; done
```

# Sorting HMM state regions
``` bash
for i in subset_noncoding*; do sort -k1,1 -k2,2n $i > sorted_$i; done
``` 

# Finding the closest quiescent position to the other HMM state
``` bash
for i in sorted*; do closestBed -a $i -b sorted_subset_noncoding_only_quescient_states_hg38.bed > closest_quescient_states_$i; done
```

# Filtering out duplicated positions
``` bash
for i in closest_quescient_states_*; do sort -k1,1 -k2,2n -k3,3n -u $i > unique_$i; done
```

# Mutation rate
``` python
# Extracted denovo mutations from Sasani et al 2019, and overlapped their denovo mutations to the specific hidden state regions

def computing_number_of_shared_mutations(mutations_df, hidden_states_coords_df, middle_distance_threshold, k, signatures_file, name_hidden_state):
    """


    Args:

    Returns:
        :class:`~pandas.DataFrame`. Table with the relative position (to the exon center)
        the number of mutations observed and expected in that position.

    """
    
    # Convert mutations into bed to intersect values
    mutations_bed = pybedtools.BedTool.from_dataframe(mutations_df)
    
    # Convert the hidden states into bed a bed file
    hidden_states_coords_bed = pybedtools.BedTool.from_dataframe(hidden_states_coords_df)
    
    # Filter mutations by the full region of interest
    my_bed = hidden_states_coords_bed.intersect(mutations_bed, wao=True, output='intersection_of_a_and_b.bed')
    
    mutations_in_range = pd.read_table(my_bed.fn, names = ['chrm_hs','start_hs', 'end_hs', 'chrm_dnm', 'start_dnm', 'end_dnm', 'ref_dnm', 'alt_dnm', 'sample_id_dnm', 'family_id_dnm', 'class_dnm', 'type_dnm', 'mut_dnm', 'col14', 'intersection_hits'],  sep="\s+", index_col=False, low_memory = False)
    mutations_in_range_new = mutations_in_range[mutations_in_range['intersection_hits'] != 0]
    print("The total amount of mutations at the middle exon-centered sequences is " +
          str(len(mutations_in_range_new.index)))
    
    # Saving the data
    mutations_in_range_new.to_csv('results/intersection_dnm_' + name_hidden_state + '.tsv', sep="\t", header=True,
                              index=False)
``` 
