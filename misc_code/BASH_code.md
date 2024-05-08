# Code for the main bioinformatic processing using snakemake
### Raphael Eisenhofer

## Get the raw data
git clone https://github.com/EisenRa/2023_BeHo_scrutinising_assembly.git && cd 2023_BeHo_scrutinising_assembly

## Preprocessing of data
I used the EHI_bioinformatics preprocessing pipeline: https://github.com/earthhologenome/EHI_bioinformatics.git
You'll just need to have snakemake installed to run the pipeline (dependencies are handled through mamba). Also note that you may want to change the workload manager (we used SLURM on our server).

```
git clone https://github.com/earthhologenome/EHI_bioinformatics.git

cd EHI_bioinformatics

mv 

snakemake \
-s 0_Code/1_Preprocess_QC.snakefile \
-j 20 \
--cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} -c {threads} -v" \
--use-conda \
--conda-frontend mamba \
--latency-wait 600

```

## Subsampling of preprocessed reads
Our samples (like any many studies) did not have even sequencing coverage. Therefore, we subsampled them with this target in mind: (image)[figures/subsampling_target.png]


#### Additionally, rename the samples so they're easier to interpret:
```
for r in `seq 1 2`;
    do while read sample file; 
        do echo mv "$file"*_M_"$r".fastq.gz "$sample"_M_"$r".fq.gz;
    done < mapping.tsv;
done

```

#### The subsampling code:
```
# Install BBmap via conda if you don't have it:
conda create -n bbmap-39.01 bbmap=39.01
conda activate bbmap-39.01

mkdir subsampled

dos2unix data/bases.tsv

while read id bases;
    do reformat.sh \
    in="$id"_M_1.fq.gz \
    in2="$id"_M_2.fq.gz \
    samplebasestarget=$bases \
    out="$id"_subsam_M_1.fq.gz \
    out2="$id"_subsam_M_2.fq.gz;
done < data/target_bases.tsv

conda deactivate

# get read counts from subsampled files
for i in ind_assembly/2_Reads/4_Host_removed/*_1.fq.gz; 
    do echo $(basename ${i/_subsam_M_1.fq.gz/}) >> sample_names.tsv && zcat $i | grep '@' | wc -l >> counts.tsv;
done

paste sample_names.tsv counts.tsv > subsampled_counts.tsv

```

# Assembly/binning/annotation
I again used the EHI_bioinformatics pipeline for the individual and coassemblies. We have a few different types of assemblies:

## Individual assembly (as it says on the tin):
```
git clone https://github.com/earthhologenome/EHI_bioinformatics.git && mv EHI_bioinformatics subsampled/ind_assembly

mkdir subsampled/ind_assembly/2_Reads/4_Host_removed/

snakemake \
-s 0_Code/2____________.snakefile \
-j 20 --cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} -c {threads} -v" \
--use-conda \```
--conda-frontend conda \
--conda-prefix /projects/mjolnir1/people/ncl550/0_software \
--latency-wait 600

snakemake \
-s 0_Code/3_MAG_dereplication_annotation_mapping.snakefile \
-j 20 --cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} -c {threads} -v" \
--use-conda \
--conda-frontend conda \
--conda-prefix /projects/mjolnir1/people/ncl550/0_software \
--latency-wait 600

#Copy and rename relevant outputs
for i in 3_Outputs/10_Final_tables/*.txt; do cp $i ../../data/tables/ind_$(basename $i); done
cp 3_Outputs/8_GTDB-tk/gtdbtk_combined_summary.tsv ../../data/tables/ind_gtdbtk_combined_summary.tsv
cp 3_Outputs/assembly_summary.tsv ../../data/tables/assembly_summary/ind_assembly_summary.tsv

```

## Coassemblies (all of group's reads assembled/binned)

#### Longitudinal coassembly (coassembly by mouse across all treatments; n=5 per mouse):
```
#Setup sym links + assembly groups:
git clone https://github.com/earthhologenome/EHI_bioinformatics.git && mv EHI_bioinformatics subsampled/coassembly_long

mkdir subsampled/coassembly_long/2_Reads/4_Host_removed/

cut -f4 data/metadata.tsv | uniq | sed '1d;' > data/mice.tsv

while read mouse; 
    do mkdir subsampled/coassembly_long/2_Reads/4_Host_removed/$mouse && ln -s `pwd`/subsampled/ind_assembly/2_Reads/4_Host_removed/"$mouse"* subsampled/coassembly_long/2_Reads/4_Host_removed/"$mouse"/; 
done < data/mice.tsv

cd subsampled/coassembly_long/

#Launch pipeline
snakemake \
-s 0_Code/2_Coassembly_Binning.snakefile \
-j 20 \
--cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} -c {threads} -v" \
--use-conda \
--conda-frontend conda \
--conda-prefix /projects/mjolnir1/people/ncl550/0_software \
--latency-wait 600

#setup for dereplication pipeline
cd 3_Outputs/5_Refined_Bins && mkdir All_metawrap_70_10_bins
cp */metawrap_70_10_bins/* All_metawrap_70_10_bins/
cp */*metawrap_70_10_bins.stats .
cd ../../
ln -s `pwd`/../ind_assembly/2_Reads/4_Host_removed/* 2_Reads/4_Host_removed/

snakemake \
-s 0_Code/3_MAG_dereplication_annotation_mapping.snakefile \
-j 20 --cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} -c {threads} -v" \
--use-conda \
--conda-frontend conda \
--conda-prefix /projects/mjolnir1/people/ncl550/0_software \
--latency-wait 600

#Copy and rename relevant outputs
for i in 3_Outputs/10_Final_tables/*.txt; do cp $i ../../data/tables/long_$(basename $i); done
cp 3_Outputs/8_GTDB-tk/gtdbtk_combined_summary.tsv ../../data/tables/long_gtdbtk_combined_summary.tsv
cp 3_Outputs/*summary.tsv ../../data/tables/assembly_summary/

```



#### cage_treatment assembly (coassembly by cage at a given treatment):
```
#Setup sym links + assembly groups:
git clone https://github.com/earthhologenome/EHI_bioinformatics.git && mv EHI_bioinformatics subsampled/coassembly_cage_treat

mkdir subsampled/coassembly_cage_treat/2_Reads/4_Host_removed

cut -f5,9 data/metadata.tsv | tr '\t' '_' | sed '1d;' | sort | uniq > data/cage_treat_folders.tsv

while read cage; 
    do mkdir subsampled/coassembly_cage_treat/2_Reads/4_Host_removed/$cage;
done < data/cage_treat_folders.tsv

cut -f5,9 data/metadata.tsv | sed '1d;' | sort | uniq > data/cage_treat.tsv

while read cage treat; 
    do ln -s `pwd`/subsampled/EHI_bioinformatics/2_Reads/4_Host_removed/"$cage"*"$treat"* subsampled/coassembly_cage_treat/2_Reads/4_Host_removed/"$cage"_"$treat"/; 
done < data/cage_treat.tsv

cd subsampled/coassembly_cage_treat/

#Launch pipeline
snakemake \
-s 0_Code/2_Coassembly_Binning.snakefile \
-j 20 \
--cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} -c {threads} -v" \
--use-conda \
--conda-frontend conda \
--conda-prefix /projects/mjolnir1/people/ncl550/0_software \
--latency-wait 600 

#setup for dereplication pipeline
cd 3_Outputs/5_Refined_Bins && mkdir All_metawrap_70_10_bins
cp */metawrap_70_10_bins/* All_metawrap_70_10_bins/
cp */*metawrap_70_10_bins.stats .
cd ../../
ln -s `pwd`/../ind_assembly/2_Reads/4_Host_removed/* 2_Reads/4_Host_removed/

snakemake \
-s 0_Code/3_MAG_dereplication_annotation_mapping.snakefile \
-j 20 --cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} -c {threads} -v" \
--use-conda \
--conda-frontend conda \
--conda-prefix /projects/mjolnir1/people/ncl550/0_software \
--latency-wait 600

#Copy and rename relevant outputs
for i in 3_Outputs/10_Final_tables/*.txt; do cp $i ../../data/tables/cage_treat_$(basename $i); done
cp 3_Outputs/8_GTDB-tk/gtdbtk_combined_summary.tsv ../../data/tables/cage_treat_gtdbtk_combined_summary.tsv
cp 3_Outputs/*summary.tsv ../../data/tables/assembly_summary/

```




#### treatment assembly (coassembly of all samples at a given treatment):
```
#Setup sym links + assembly groups:
git clone https://github.com/earthhologenome/EHI_bioinformatics.git && mv EHI_bioinformatics subsampled/coassembly_treat

mkdir subsampled/coassembly_treat/2_Reads/4_Host_removed

cut -f9 data/metadata.tsv | sed '1d;' | sort | uniq > data/treatments.tsv

while read treatment; 
    do mkdir subsampled/coassembly_treat/2_Reads/4_Host_removed/$treatment && ln -s `pwd`/subsampled/ind_assembly/2_Reads/4_Host_removed/*"$treatment"* subsampled/coassembly_treat/2_Reads/4_Host_removed/"$treatment"/; 
done < data/treatments.tsv

cd subsampled/coassembly_treat

snakemake \
-s 0_Code/2_Coassembly_Binning.snakefile \
-j 20 \
--cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} -c {threads} -v" \
--use-conda \
--conda-frontend conda \
--conda-prefix /projects/mjolnir1/people/ncl550/0_software \
--latency-wait 600 

#setup for dereplication pipeline
cd 3_Outputs/5_Refined_Bins && mkdir All_metawrap_70_10_bins
cp */metawrap_70_10_bins/* All_metawrap_70_10_bins/
cp */*metawrap_70_10_bins.stats .
cd ../../
ln -s `pwd`/../ind_assembly/2_Reads/4_Host_removed/* 2_Reads/4_Host_removed/

snakemake \
-s 0_Code/3_MAG_dereplication_annotation_mapping.snakefile \
-j 20 --cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} -c {threads} -v" \
--use-conda \
--conda-frontend conda \
--conda-prefix /projects/mjolnir1/people/ncl550/0_software \
--latency-wait 600

#Copy and rename relevant outputs
for i in 3_Outputs/10_Final_tables/*.txt; do cp $i ../../data/tables/treat_$(basename $i); done
cp 3_Outputs/8_GTDB-tk/gtdbtk_combined_summary.tsv ../../data/tables/treat_gtdbtk_combined_summary.tsv
cp 3_Outputs/*summary.tsv ../../data/tables/assembly_summary/

```

## Cobinning (Individual assemblies -> binning with all samples in group)
Using individual assemblies from above with a different snakefile.

#### Cobinning by treatment
```
#Setup binning group input:
git clone https://github.com/earthhologenome/EHI_bioinformatics.git && mv EHI_bioinformatics subsampled/cobinning_treat

cd subsampled/cobinning_treat

while read treatment; 
    do for assembly in ../ind_assembly/3_Outputs/2_Assemblies/*"$treatment"*.fasta; 
        do for reads in ../ind_assembly/2_Reads/4_Host_removed/*"$treatment"*_1.fq.gz; 
            do echo $(basename $assembly) >> assemblies.tsv && echo $(basename ${reads/_1.fq.gz/}) >> reads.tsv; 
            done; 
        done; 
    done < ../../data/treatments.tsv

paste assemblies.tsv reads.tsv > valid_combinations_temp.tsv

#Filter by only a subset of mice (to save computation)
while read mouse;
    do grep "$mouse".._subsam_contigs valid_combinations_temp.tsv;
    done < ../../data/mice_subset.tsv > valid_combinations.tsv

snakemake -s 0_Code/cobinning.smk -j 120 --cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} -c {threads} -v" --use-conda --conda-frontend conda --conda-prefix /projects/mjolnir1/people/ncl550/0_software --latency-wait 60 -n

```

#### Cobinning by cage_treatment (cobinning by cage at a given treatment):
```
#Setup binning group input:
git clone https://github.com/earthhologenome/EHI_bioinformatics.git && mv EHI_bioinformatics subsampled/cobinning_cage_treat

cd subsampled/cobinning_cage_treat

while read cage treatment; 
    do for assembly in ../ind_assembly/3_Outputs/2_Assemblies/*"$cage"*"$treatment"*.fasta;
        do for reads in ../ind_assembly/2_Reads/4_Host_removed/"$cage"*"$treatment"_subsam_M_1.fq.gz; 
            do echo $(basename $assembly) >> assemblies.tsv && echo $(basename ${reads/_1.fq.gz/}) >> reads.tsv; 
            done; 
        done; 
    done < ../../data/cage_treat.tsv

paste assemblies.tsv reads.tsv > valid_combinations.tsv

snakemake -s 0_Code/cobinning.smk -j 120 --cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} -c {threads} -v" --use-conda --conda-frontend conda --conda-prefix /projects/mjolnir1/people/ncl550/0_software --latency-wait 60 -n

```

#### Longitudinal cobinning (cobinning by mouse across all treatments; n=5 per mouse)
```
#Setup binning group input:
git clone https://github.com/earthhologenome/EHI_bioinformatics.git && mv EHI_bioinformatics subsampled/cobinning_long

cd subsampled/cobinning_long

while read mouse; 
    do for assembly in ../ind_assembly/3_Outputs/2_Assemblies/"$mouse"*.fasta; 
        do for reads in ../ind_assembly/2_Reads/4_Host_removed/"$mouse"*_1.fq.gz; 
            do echo $(basename $assembly) >> assemblies.tsv && echo $(basename ${reads/_1.fq.gz/}) >> reads.tsv; 
            done; 
        done; 
    done < ../../data/mice.tsv

paste assemblies.tsv reads.tsv > valid_combinations.tsv

snakemake -s 0_Code/cobinning.smk -j 120 --cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} -c {threads} -v" --use-conda --conda-frontend conda --conda-prefix /projects/mjolnir1/people/ncl550/0_software --latency-wait 60

```

## Multisplit (Indiviual assemblies concatenated by group and binned with vamb)
We're using Vamb here for this one. Individual assemblies as above, with a new snakefile.

# download checkm-db
mkdir checkm_db && cd checkm_db
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvzf checkm_data_2015_01_16.tar.gz && rm checkm_data_2015_01_16.tar.gz && cd ../


#### Multisplit binning by treatment
```
#Setup binning group input:
git clone https://github.com/earthhologenome/EHI_bioinformatics.git && mv EHI_bioinformatics subsampled/vamb_treat

cd subsampled/vamb_treat

while read treatment; 
    do for assembly in ../ind_assembly/3_Outputs/2_Assemblies/*"$treatment"*.fasta; 
        do echo $(basename $assembly) >> assemblies.tsv 
        done; 
    done < ../../data/treatments.tsv

while read treatment; 
    do for reads in ../ind_assembly/2_Reads/4_Host_removed/*"$treatment"*_1.fq.gz; 
        do echo $(basename ${reads/_1.fq.gz/}) >> reads.tsv && echo $treatment >> treatment.tsv 
        done; 
    done < ../../data/treatments.tsv

paste assemblies.tsv reads.tsv treatment.tsv > valid_combinations.tsv

snakemake -s 0_Code/vamb_multisplit.smk -j 120 --cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} -c {threads} -v" --use-conda --conda-frontend conda --conda-prefix /projects/mjolnir1/people/ncl550/0_software --latency-wait 60

```

#### Multisplit binning by cage_treatment
```
#Setup binning group input:
git clone https://github.com/earthhologenome/EHI_bioinformatics.git && mv EHI_bioinformatics subsampled/vamb_cage_treat

cd subsampled/vamb_cage_treat

while read cage treatment; 
    do for assembly in ../ind_assembly/3_Outputs/2_Assemblies/*"$cage"*"$treatment"*.fasta;
        do echo $(basename $assembly) >> assemblies.tsv; 
        done; 
    done < ../../data/cage_treat.tsv

while read cage treatment; 
    do for reads in ../ind_assembly/2_Reads/4_Host_removed/"$cage"*"$treatment"_subsam_M_1.fq.gz; 
        do echo $(basename ${reads/_1.fq.gz/}) >> reads.tsv && echo "$cage"_"$treatment" >> cage_treat.tsv;
        done;
    done < ../../data/cage_treat.tsv

paste assemblies.tsv reads.tsv cage_treat.tsv > valid_combinations.tsv 

snakemake -s 0_Code/vamb_multisplit.smk -j 120 --cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} -c {threads} -v" --use-conda --conda-frontend conda --conda-prefix /projects/mjolnir1/people/ncl550/0_software --latency-wait 60

```


#### Multisplit binning by longitudinal
```
#Setup binning group input:
git clone https://github.com/earthhologenome/EHI_bioinformatics.git && mv EHI_bioinformatics subsampled/vamb_long

cd subsampled/vamb_long

while read mouse; 
    do for assembly in ../ind_assembly/3_Outputs/2_Assemblies/"$mouse"*.fasta; 
        do echo $(basename $assembly) >> assemblies.tsv; 
        done; 
    done < ../../data/mice.tsv

while read mouse; 
    do for reads in ../ind_assembly/2_Reads/4_Host_removed/"$mouse"*_1.fq.gz;
        do echo $(basename ${reads/_1.fq.gz/}) >> reads.tsv && echo $mouse >> mouse.tsv;
        done
    done < ../../data/mice.tsv

paste assemblies.tsv reads.tsv mouse.tsv > valid_combinations.tsv

snakemake -s 0_Code/vamb_multisplit.smk -j 120 --cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} -c {threads} -v" --use-conda --conda-frontend conda --conda-prefix /projects/mjolnir1/people/ncl550/0_software --latency-wait 60

```


#### Dereplicate all MAGs from all assembly methods + map
```
git clone https://github.com/earthhologenome/EHI_bioinformatics.git && mv EHI_bioinformatics subsampled/all_dereplicated

#Setup sym links for reads
cd subsampled/all_dereplicated/ && mkdir 2_Reads/4_Host_removed

for i in `pwd`/../ind_assembly/2_Reads/4_Host_removed/*.gz; do ln -s $i 2_Reads/4_Host_removed/; done

#Setup sym links for bins + .stats files
mkdir 3_Outputs && mkdir 3_Outputs/5_Refined_Bins && mkdir 3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins

for i in ../ind_assembly/3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/*.fa.gz;
    do mv $i ${i/bins\//bins\/individual_}; 
done

ln -s `pwd`/../ind_assembly/3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/*.gz 3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/

for i in ../coassembly_cage_treat/3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/*.fa.gz;
    do mv $i ${i/bins\//bins\/coassembly_cage_treat_}; 
done

for i in ../coassembly_treat/3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/*.fa.gz;
    do mv $i ${i/bins\//bins\/coassembly_treat_}; 
done

for i in ../coassembly_long/3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/*.fa.gz;
    do mv $i ${i/bins\//bins\/coassembly_long_}; 
done

ln -s `pwd`/../coassembly_*/3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/*.gz 3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/

for i in ../cobinning_cage_treat/3_Outputs/drep/All_metawrap_70_10_bins/*.fa.gz;
    do mv $i ${i/bins\//bins\/cobinning_cage_treat_}; 
done

for i in ../cobinning_treat/3_Outputs/drep/All_metawrap_70_10_bins/*.fa.gz;
    do mv $i ${i/bins\//bins\/cobinning_treat_}; 
done

for i in ../cobinning_long/3_Outputs/drep/All_metawrap_70_10_bins/*.fa.gz;
    do mv $i ${i/bins\//bins\/cobinning_long_}; 
done

ln -s `pwd`/../cobinning_*/3_Outputs/drep/All_metawrap_70_10_bins/*.gz 3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/

for i in ../vamb_long/3_Outputs/drep/All_bins/*.fa;
    do mv $i ${i/bins\//bins\/vamb_long_}; 
done

for i in ../vamb_treat/3_Outputs/drep/All_bins/*.fa;
    do mv $i ${i/bins\//bins\/vamb_treat_}; 
done

for i in ../vamb_cage_treat/3_Outputs/drep/All_bins/*.fa;
    do mv $i ${i/bins\//bins\/vamb_cage_treat_}; 
done

for i in `pwd`/../vamb_*/3_Outputs/drep/All_bins/*.fa;
    do ln -s $i 3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/;
done

```

## Stats files
```
cat ../ind_assembly/3_Outputs/5_Refined_Bins/genome_info.csv | tr ',' '\t' | sed 's@3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/@individual_@g' > 3_Outputs/5_Refined_Bins/individual.stats

cat ../coassembly_long/3_Outputs/5_Refined_Bins/genome_info.csv | tr ',' '\t' | sed 's@3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/@coassembly_long_@g' > 3_Outputs/5_Refined_Bins/coassembly_long.stats

cat ../coassembly_treat/3_Outputs/5_Refined_Bins/genome_info.csv | tr ',' '\t' | sed 's@3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/@coassembly_treat_@g' > 3_Outputs/5_Refined_Bins/coassembly_treat.stats

cat ../coassembly_cage_treat/3_Outputs/5_Refined_Bins/genome_info.csv | tr ',' '\t' | sed 's@3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/@coassembly_cage_treat_@g' > 3_Outputs/5_Refined_Bins/coassembly_cage_treat.stats

cat ../cobinning_long/3_Outputs/drep/genome_info.csv | tr ',' '\t' | sed 's@3_Outputs/drep/All_metawrap_70_10_bins/@cobinning_long_@g' > 3_Outputs/5_Refined_Bins/cobinning_long.stats

cat ../cobinning_treat/3_Outputs/drep/genome_info.csv | tr ',' '\t' | sed 's@3_Outputs/drep/All_metawrap_70_10_bins/@cobinning_treat_@g' > 3_Outputs/5_Refined_Bins/cobinning_treat.stats

cat ../cobinning_cage_treat/3_Outputs/drep/genome_info.csv | tr ',' '\t' | sed 's@3_Outputs/drep/All_metawrap_70_10_bins/@cobinning_cage_treat_@g' > 3_Outputs/5_Refined_Bins/cobinning_cage_treat.stats

cat ../vamb_long/3_Outputs/drep/genome_info.csv | tr ',' '\t' | sed 's@3_Outputs/drep/All_bins/@vamb_long_@g' > 3_Outputs/5_Refined_Bins/vamb_long.stats

cat ../vamb_treat/3_Outputs/drep/genome_info.csv | tr ',' '\t' | sed 's@3_Outputs/drep/All_bins/@vamb_treat_@g' > 3_Outputs/5_Refined_Bins/vamb_treat.stats

cat ../vamb_cage_treat/3_Outputs/drep/genome_info.csv | tr ',' '\t' | sed 's@3_Outputs/drep/All_bins/@vamb_cage_treat_@g' > 3_Outputs/5_Refined_Bins/vamb_cage_treat.stats


## Run the pipeline

snakemake \
-s 0_Code/3_MAG_dereplication_annotation_mapping.snakefile \
-j 20 --cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} -c {threads} -v" \
--use-conda \
--conda-frontend conda \
--conda-prefix /projects/mjolnir1/people/ncl550/0_software \
--latency-wait 600 -n

cd ../
```

## Run DRAM on each treatment's dereplicated genomes

```
#setup directory
git clone https://github.com/earthhologenome/EHI_bioinformatics.git && mv EHI_bioinformatics dram
mkdir dram/3_Outputs/ && mkdir dram/3_Outputs/5_Refined_Bins && mkdir dram/3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins

#create sym links

for i in vamb_*; do for mag in $i/3_Outputs/drep/98/dereplicated_genomes/*.fa.gz; do ln -s `pwd`/$mag dram/3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/"$i"_$(basename $mag); done; done

for i in cobinning_*; do for mag in $i/3_Outputs/drep/98/dereplicated_genomes/*.fa.gz; do ln -s `pwd`/$mag dram/3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/"$i"_$(basename $mag); done; done

for i in coassembly_*; do for mag in $i/3_Outputs/7_Dereplication/dereplicated_genomes/*.fa.gz; do ln -s `pwd`/$mag dram/3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/"$i"_$(basename $mag); done; done

for i in ind_assembly/3_Outputs/7_Dereplication/dereplicated_genomes/*.fa.gz; do ln -s `pwd`/$i dram/3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/individual_$(basename $i); done

for i in all_dereplicated/3_Outputs/7_Dereplication/dereplicated_genomes/*.fa.gz; do ln -s `pwd`/$i dram/3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/ALL_dreped_$(basename $i); done

cd dram

# run snakemake
snakemake -s 0_Code/DRAM.snakefile -j 40 --cluster "sbatch --mem {resources.mem_gb}G -c {threads} --time {resources.time} -v" --use-conda --conda-frontend conda --conda-prefix /projects/mjolnir1/people/ncl550/0_software --latency-wait 600

cd ../

```


#### Run GUNC on each set of dereplicated MAGs
```
#!/bin/sh
#SBATCH -c 16 --mem 72G --time 08:00:00 # number of cores

gunc run \
    --file_suffix .fa.gz \
    --input_dir dram/3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/ \
    -t 16 \
    -r gunc_db/gunc_db_progenomes2.1.dmnd


```


## rename files

```
mkdir tables

# get ncontigs from bins:
for i in ind_assembly/3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/*.gz; 
    do echo $(basename $i) >> names.tsv && zcat $i | grep '>' | wc -l >> n_contigs.tsv; 
done

for i in coassembly_*/3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/*.gz; 
    do echo $(basename $i) >> names.tsv && zcat $i | grep '>' | wc -l >> n_contigs.tsv; 
done

for i in cobinning_*/3_Outputs/drep/All_metawrap_70_10_bins/*.gz; 
    do echo $(basename ${i/}) >> names.tsv && zcat $i | grep '>' | wc -l >> n_contigs.tsv; 
done

paste names.tsv n_contigs.tsv > ncontigs.tsv

cp ind_assembly/3_Outputs/assembly_summary.tsv tables/

for i in coassembly_cage_treat/3_Outputs/*_summary.tsv; 
    do cp $i tables/$(basename ${i/coassembly_summary.tsv/coassembly_cage_treat.tsv});
done

for i in coassembly_long/3_Outputs/*_summary.tsv; 
    do cp $i tables/$(basename ${i/coassembly_summary.tsv/coassembly_long.tsv});
done

for i in coassembly_treat/3_Outputs/*_summary.tsv; 
    do cp $i tables/$(basename ${i/coassembly_summary.tsv/coassembly_treat.tsv});
done

for i in cobinning_cage_treat/3_Outputs/refinement/*/*_metawrap_70_10_bins.stats;
    do cp $i tables/$(basename ${i/_subsam_contigs.fasta_metawrap/_cobinning_cage_treat_metawrap});
done

for i in cobinning_treat/3_Outputs/refinement/*/*_metawrap_70_10_bins.stats;
    do cp $i tables/$(basename ${i/_subsam_contigs.fasta_metawrap/_cobinning_treat_metawrap});
done

for i in cobinning_long/3_Outputs/refinement/*/*_metawrap_70_10_bins.stats;
    do cp $i tables/$(basename ${i/_subsam_contigs.fasta_metawrap/_cobinning_long_metawrap});
done

cp ind_assembly/3_Outputs/5_Refined_Bins/All_bins.stats tables/individual_metawrap_70_10_bins.stats

for i in vamb_long/3_Outputs/checkm/*.tsv;
    do cp $i tables/$(basename ${i/_checkm/_vamb_long_checkm});
done

for i in vamb_cage_treat/3_Outputs/checkm/*.tsv;
    do cp $i tables/$(basename ${i/_checkm/_vamb_cage_treat_checkm});
done

for i in vamb_treat/3_Outputs/checkm/*.tsv;
    do cp $i tables/$(basename ${i/_checkm/_vamb_treat_checkm});
done


# Get bin info from vamb bins (requires bbmap to be installed)

mkdir vamb_stats

for i in vamb_treat/3_Outputs/vamb/*/vamb_bins/*.fa; do
    stats.sh addname=true format=3 in=$i out=vamb_stats/$(basename ${i/.fa/.tsv});
done

for i in vamb_stats/*.tsv; 
    do sed '1d;' $i >> vamb_stats_temp.tsv;
done

head -1 vamb_stats/CD_vae_10.tsv >> vamb_stats_header.tsv

cat vamb_stats_header.tsv vamb_stats_temp.tsv > vamb_mag_stats.tsv
gzip vamb_mag_stats.tsv

mv vamb_mag_stats.tsv tables/

cp coassembly_*/3_Outputs/5_Refined_Bins/*.stats tables/

cp cobinning_*/3_Outputs/refinement/*/*_metawrap_70_10_bins.stats tables/

for i in cobinning*;
    do cp "$i"/3_Outputs/gtdbtk/gtdbtk_combined_summary.tsv tables/$(basename "$i"_gtdbtk_combined_summary.tsv);
done

for i in vamb*;
    do cp "$i"/3_Outputs/gtdbtk/gtdbtk_combined_summary.tsv tables/$(basename "$i"_gtdbtk_combined_summary.tsv);
done

for i in coassembly*;
    do cp "$i"/3_Outputs/8_GTDB-tk/gtdbtk_combined_summary.tsv tables/$(basename "$i"_gtdbtk_combined_summary.tsv);
done

cp ind_assembly/3_Outputs/8_GTDB-tk/gtdbtk_combined_summary.tsv tables/individual_gtdbtk_combined_summary.tsv

gzip tables/*


# rename count tables:
mkdir count_tables

for i in *; do cp $i/3_Outputs/10_Final_tables/mapping_rate.txt count_tables/"$i"_mapping.tsv; done

for i in *; do cp $i/3_Outputs/coverm/count_table.tsv count_tables/"$i"_count_table.tsv; done


gzip count_tables/*


# rename full_trees:
for i in vamb*; 
    do cp $i/3_Outputs/gtdbtk/classify/gtdbtk.bac120.classify.tree full_trees/"$i"_full_tree.tree;
done

for i in cobinning*; 
    do cp $i/3_Outputs/gtdbtk/classify/gtdbtk.bac120.classify.tree full_trees/"$i"_full_tree.tree;
done

for i in coassembly*; 
    do cp $i/3_Outputs/8_GTDB-tk/classify/gtdbtk.bac120.classify.tree full_trees/"$i"_full_tree.tree;
done

cp ind_assembly/3_Outputs/8_GTDB-tk/classify/gtdbtk.bac120.classify.tree full_trees/individual_full_tree.tree

cp all_dereplicated/3_Outputs/8_GTDB-tk/classify/gtdbtk.bac120.classify.tree full_trees/all_dereplicated_full_tree.tree

gzip full_trees/*

```

# GUNC
```
#!/bin/sh
#SBATCH -c 16 --mem 96G --time 24:00:00 # number of cores

gunc run \
    --file_suffix .fa.gz \
    -o GUNC_out \
    --input_dir dram/3_Outputs/5_Refined_Bins/All_metawrap_70_10_bins/ \
    -t 16 \
    -r gunc_db/gunc_db_progenomes2.1.dmnd
```

# Collect runtime stats from snakemake
```
# Run in each snakemake directory:
snakemake -s 0_Code/*.snakefile --report --{approach}_report.html

# Some hacky BASH-fu to extract the values from the .html files into a usable format:
for i in *.html; 
    do grep 'runtimes_spec' $i | tr { \\n | grep 'runtime' | sed '1d;' | grep -v 'field' | grep -v 'label' | grep -v 'vega' | cut -f2,4 -d ' ' | sed 's/"//g' | sed 's/,//g' | sed 's/]//g' | sed 's/}//g' > ${i/.html/.tsv};
done

# Snakemake with benchmark was run for the individual and coassembly pipelines, so collect them into a table format too:
for i in ind_assembly/3_Outputs/0_Logs/*MAG_mapping.benchmark.tsv; do echo "individual" >> ind_strategy.tsv && echo "mag_catalogue_mapping" >> ind_job.tsv && sed '1d;' $i | cut -f1 >> ind_time.tsv; done

paste ind_job.tsv ind_time.tsv ind_strategy.tsv > ind_benchmark.tsv

for i in coassembly_treat/3_Outputs/0_Logs/*MAG_mapping.benchmark.tsv; do echo "coassembly_treat" >> coassembly_strategy.tsv && echo "mag_catalogue_mapping" >> coassembly_job.tsv && sed '1d;' $i | cut -f1 >> coassembly_time.tsv; done

for i in coassembly_cage_treat/3_Outputs/0_Logs/*MAG_mapping.benchmark.tsv; do echo "coassembly_cage_treat" >> coassembly_strategy.tsv && echo "mag_catalogue_mapping" >> coassembly_job.tsv && sed '1d;' $i | cut -f1 >> coassembly_time.tsv; done

for i in coassembly_long/3_Outputs/0_Logs/*MAG_mapping.benchmark.tsv; do echo "coassembly_long" >> coassembly_strategy.tsv && echo "mag_catalogue_mapping" >> coassembly_job.tsv && sed '1d;' $i | cut -f1 >> coassembly_time.tsv; done

for i in coassembly_treat/3_Outputs/0_Logs/*coassembly_mapping.benchmark.tsv; do echo "coassembly_treat" >> coassembly_strategy.tsv && echo "assembly_mapping" >> coassembly_job.tsv && sed '1d;' $i | cut -f1 >> coassembly_time.tsv; done

for i in coassembly_cage_treat/3_Outputs/0_Logs/*coassembly_mapping.benchmark.tsv; do echo "coassembly_cage_treat" >> coassembly_strategy.tsv && echo "assembly_mapping" >> coassembly_job.tsv && sed '1d;' $i | cut -f1 >> coassembly_time.tsv; done

for i in coassembly_long/3_Outputs/0_Logs/*coassembly_mapping.benchmark.tsv; do echo "coassembly_long" >> coassembly_strategy.tsv && echo "assembly_mapping" >> coassembly_job.tsv && sed '1d;' $i | cut -f1 >> coassembly_time.tsv; done

for i in coassembly_treat/3_Outputs/0_Logs/*coassembly.benchmark.tsv; do echo "coassembly_treat" >> coassembly_strategy.tsv && echo "Assembly" >> coassembly_job.tsv && sed '1d;' $i | cut -f1 >> coassembly_time.tsv; done

for i in coassembly_cage_treat/3_Outputs/0_Logs/*coassembly.benchmark.tsv; do echo "coassembly_cage_treat" >> coassembly_strategy.tsv && echo "Assembly" >> coassembly_job.tsv && sed '1d;' $i | cut -f1 >> coassembly_time.tsv; done

for i in coassembly_long/3_Outputs/0_Logs/*coassembly.benchmark.tsv; do echo "coassembly_long" >> coassembly_strategy.tsv && echo "Assembly" >> coassembly_job.tsv && sed '1d;' $i | cut -f1 >> coassembly_time.tsv; done

paste coassembly_job.tsv coassembly_time.tsv coassembly_strategy.tsv > coassembly_benchmark.tsv

```


# Run SingleM microbial fraction
```
#singlem version 0.16, metapackage: S3.2.1.GTDB_r214.metapackage_20231006.smpkg

for i in ind_assembly/2_Reads/4_Host_removed/*M_1.fq.gz;
do singlem pipe \
        -1 $i \
        -2 ${i/_1.fq/_2.fq} \
        --threads 12 \
        --archive-otu-table ${i/M_1.fq.gz/pipe.json} \
        --otu-table ${i/M_1.fq.gz/pipe.tsv} \
        --taxonomic-profile ${i/M_1.fq.gz/pipe_profile.tsv}
done

for i in ind_assembly/2_Reads/4_Host_removed/*pipe_profile.tsv;
do singlem microbial_fraction \
        -p $i \
        -1 ${i/pipe_profile.tsv/M_1.fq.gz} \
        -2 ${i/pipe_profile.tsv/M_2.fq.gz} \
        --output-tsv ${i/pipe_profile.tsv/read_fraction.tsv} \
        --output-per-taxon-read-fractions ${i/pipe_profile.tsv/pertaxon_rf.tsv}
done

```