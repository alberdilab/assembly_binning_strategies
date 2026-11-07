configfile: "0_Code/configs/2_Assembly_Binning_config.yaml"

import pandas as pd

colnames = ['assembly', 'sample', 'treatment']
df = pd.read_csv("valid_combinations.tsv", sep="\t", names=colnames, header=None)

valid_combinations = set(
    (row["assembly"], row["sample"], row["treatment"]) for _, row in df.iterrows()
)

# Function to filter combinations based on assembly
def filter_combinations(treatment):
    return [combo for combo in valid_combinations if combo[2] == treatment]

#Define wildcards for later
ASSEM = list(set(df.iloc[0:, 0].tolist()))
SAM = list(set(df.iloc[0:, 1].tolist()))
TREAT = list(set(df.iloc[0:, 2].tolist()))

# Define desired outputs
rule all:
    input:
        expand(
            os.path.join(
                "3_Outputs/combined_contigs/", "{combo[2]}/{combo[2]}_combined_contigs.fasta.gz"
            ),
            combo=valid_combinations,
        ),
        expand(
            os.path.join(
                "3_Outputs/checkm/", "{combo[2]}_checkm.tsv"
            ),
            combo=valid_combinations,
        ),
        count_table=os.path.join(
            "3_Outputs/coverm/", "count_table.tsv"
        )

# Job rules
################################################################################
### Combine contigs from different individual assemblies
rule combine_contigs:
    input:
        contigs=lambda wildcards: expand(
            os.path.join(
                "../ind_assembly/3_Outputs/2_Assemblies/", "{assembly}"
            ),
            assembly=[combo[0] for combo in filter_combinations(wildcards.treatment)],
        ),
    output:
        os.path.join("3_Outputs/combined_contigs/" "{treatment}/{treatment}_combined_contigs.fasta.gz")
    # conda:
    #     f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 1
    resources:
        mem_gb=16,
        time='00:30:00'
    shell:
        """
        source activate /projects/mjolnir1/people/ncl550/0_software/miniconda/envs/vamb4.1.3

        python /projects/mjolnir1/people/ncl550/1_projects/binners/vamb/vamb/src/concatenate.py \
        {output} {input.contigs}

        """

################################################################################
### Index multisplit contigs
rule multisplit_index:
    input:
        os.path.join(
            "3_Outputs/combined_contigs/" "{treatment}/{treatment}_combined_contigs.fasta.gz"
        )
    output:
        os.path.join(
            "3_Outputs/combined_contigs/" "{treatment}/{treatment}_combined_contigs.fasta.gz.rev.2.bt2l"
        )
    conda:
        "conda_envs/2_Assembly_Binning.yaml"
    threads: 16
    resources:
        mem_gb=64,
        time="04:00:00",
    shell:
        """
        # Index the coassembly
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {input} {input} 
        """

rule mapping:
    input:
        contigs=os.path.join(
            "3_Outputs/combined_contigs/" "{treatment}/{treatment}_combined_contigs.fasta.gz"
            ),
        index=os.path.join(
            "3_Outputs/combined_contigs/" "{treatment}/{treatment}_combined_contigs.fasta.gz.rev.2.bt2l"
            ),
        r1=os.path.join(
            "../ind_assembly/2_Reads/4_Host_removed/", "{sample}_1.fq.gz"
            ),
        r2=os.path.join(
            "../ind_assembly/2_Reads/4_Host_removed/", "{sample}_2.fq.gz"
            )
    output:
        temp(os.path.join(
            "3_Outputs/bams/", "{treatment}/{treatment}_{sample}.bam"
            ))
    wildcard_constraints:
        sample="(.+)_subsam_[^_]+"
    threads: 8
    resources:
        mem_gb="48",
        time='02:00:00'
    conda:
        "conda_envs/2_Assembly_Binning.yaml"
    shell:
        """
        # Map reads to assembly using Bowtie2
        bowtie2 \
            --time \
            --threads {threads} \
            -x {input.contigs} \
            -1 {input.r1} \
            -2 {input.r2} \
        | samtools sort -@ {threads} -o {output}
        """

################################################################################
### Bin contigs using vamb
rule vamb_multisplit:
    input:
        bam=lambda wildcards: expand(
            os.path.join(
                "3_Outputs/bams/", "{treatment}/{treatment}_{sample}.bam"
            ),
            treatment=wildcards.treatment,
            sample=[combo[1] for combo in filter_combinations(wildcards.treatment)],
        ),
        contigs=os.path.join(
            "3_Outputs/combined_contigs/" "{treatment}/{treatment}_combined_contigs.fasta.gz"
            ),
    output:
        os.path.join(
             "3_Outputs/vamb/", "{treatment}/vamb_complete"
            ),
    params:
        outdir=os.path.join("3_Outputs/vamb/" + "{treatment}")
    threads: 16
    # conda:
    #     f"{config['codedir']}/conda_envs/vamb.yaml"
    resources:
        mem_gb=64,
        time="12:00:00"
    shell:
        """
            mkdir -p 3_Outputs/vamb
            rm -rf {params.outdir}
            source activate /projects/mjolnir1/people/ncl550/0_software/miniconda/envs/vamb4.1.3
            # Run vamb
            vamb \
            --outdir {params.outdir} \
            --fasta {input.contigs} \
            --bamfiles {input.bam} \
            -o C \
            --minfasta 200000

            # rename bins from .fna to .fa
            mkdir {params.outdir}/vamb_bins

            for i in {params.outdir}/bins/*/*.fna;
                do mv $i {params.outdir}/vamb_bins/$(basename ${{i/.fna/.fa}});
            done

            # Create output for the next rule
            touch {output}
        """

################################################################################
### checkm
rule checkm:
    input:
        os.path.join(
             "3_Outputs/vamb/", "{treatment}/vamb_complete"
            ),
    output:
        os.path.join(
            "3_Outputs/checkm/", "{treatment}_checkm.tsv"
            )
    wildcard_constraints:
        sample="(.+)_subsam_[^_]+"
    threads: 16
    resources:
        mem_gb="90",
        time='08:00:00'
    params:
        vambdir=os.path.join("3_Outputs/vamb/" + "{treatment}")
    conda:
        "conda_envs/checkm.yaml"
    shell:
        """
        #rename bins
        for i in {params.vambdir}/vamb_bins/*.fa;
            do mv $i {params.vambdir}/vamb_bins/{wildcards.treatment}_"$(basename $i)";
        done

        # Set checkm-db path
        export checkmdb=../checkm_db
        printf $checkmdb | checkm data setRoot

        mkdir -p {params.vambdir}/tmpdir

        # run checkm
        checkm lineage_wf \
            -x fa \
            {params.vambdir}/vamb_bins \
            {params.vambdir}/checkm \
            --tab_table \
            --tmpdir {params.vambdir}/tmpdir \
            -f {output} \
            -t {threads} \
            --pplacer_threads 2

        """

################################################################################
### dRep
rule dRep:
    input:
        lambda wildcards: expand(
            os.path.join(
                "3_Outputs/checkm/", "{treatment}_checkm.tsv",
            ),
            treatment=TREAT
        )
    output:
        os.path.join(
            "3_Outputs/drep/", "mag_catalogue.fa.gz"
        )
    wildcard_constraints:
        sample="(.+)_subsam_[^_]+"
    threads: 16
    resources:
        mem_gb="128",
        time='08:00:00'
    params:
        vambdir="3_Outputs/vamb/"
    # conda:
    #     "conda_envs/drep.yaml"
    shell:
        """
        source activate /projects/mjolnir1/people/ncl550/0_software/e4fd3b3a5e8607bdeda760e6c6c0149b_
        module load bbmap/39.01
        # Copy MAGs into a shared folder
        mkdir -p 3_Outputs/drep/ && mkdir -p 3_Outputs/drep/All_bins
        cp {params.vambdir}/*/vamb_bins/*.fa 3_Outputs/drep/All_bins/

        # Parse/collate metawrap stats files for compatibility with dRep genomeinfo:
        echo -e "genome,completeness,contamination" > 3_Outputs/drep/header.txt
        for i in {input};
            do sed '1d;' $i | cut -f 1,12,13 --output-delimiter ',' >> 3_Outputs/drep/bin_info.txt;
        done
        sed -i'' 's@^@3_Outputs/drep/All_bins/@g' 3_Outputs/drep/bin_info.txt
        sed -i'' 's/,/.fa,/' 3_Outputs/drep/bin_info.txt
        cat 3_Outputs/drep/header.txt 3_Outputs/drep/bin_info.txt > 3_Outputs/drep/genome_info.csv

        # run dRep
        dRep dereplicate \
            3_Outputs/drep/98 \
            -p {threads} \
            -comp 70 \
            -con 10 \
            -sa 0.98 \
            --genomeInfo 3_Outputs/drep/genome_info.csv \
            -g 3_Outputs/drep/All_bins/*.fa

        # Reformat MAG headers for CoverM
        for mag in 3_Outputs/drep/98/dereplicated_genomes/*.fa;
            do rename.sh \
                in=$mag \
                out=3_Outputs/drep/98/dereplicated_genomes/$(basename ${{mag/.fa/.fa.gz}}) \
                zl=9 \
                prefix=$(basename ${{mag/.fa/^}});
        done

        # Concatentate the renamed, dereplicated MAGs
        cat 3_Outputs/drep/98/dereplicated_genomes/*.fa.gz > 3_Outputs/drep/mag_catalogue.fa.gz

        """

################################################################################
### Annotate dereplicated MAGs with gtdb-tk taxonomy:
rule gtdbtk:
    input:
        "3_Outputs/drep/98/figures/Primary_clustering_dendrogram.pdf"
    output:
        "3_Outputs/gtdbtk/classify/gtdbtk.bac120.summary.tsv"
    params:
        GTDB_data = expand("{GTDB_data}", GTDB_data=config['GTDB_data']),
        outdir = "3_Outputs/gtdbtk/",
        bins = "3_Outputs/drep/98/dereplicated_genomes"
    conda:
        "conda_envs/3_GTDB-tk.yaml"
    threads:
        16
    resources:
        mem_gb=456,
        time='08:00:00'
    benchmark:
        "3_Outputs/0_Logs/gtdbtk.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/gtdbtk.log"
    message:
        "Running gtdb-tk to taxonomically annotate MAGs"
    shell:
        """
        # Specify path to reference data:
        export GTDBTK_DATA_PATH={params.GTDB_data}

        mkdir {params.outdir}/tmp

        # Run GTDB-tk:
        gtdbtk classify_wf \
        --genome_dir {params.bins} \
        --extension ".fa.gz" \
        --out_dir {params.outdir} \
        --tmpdir {params.outdir}/tmp \
        --cpus {threads} \
        --skip_ani_screen \
        --full_tree

        # Create a merged summary output for DRAM:
        if [ -s "{params.outdir}/classify/gtdbtk.ar53.summary.tsv" ]
        then
        sed '1d;' {params.outdir}/classify/gtdbtk.ar53.summary.tsv > {params.outdir}/ar53.tsv
        cat {output} {params.outdir}/ar53.tsv > {params.outdir}/gtdbtk_combined_summary.tsv
        rm {params.outdir}/ar53.tsv

        # Otherwise, just use the bacterial summary (if no archaeal bins)
        else
        cat {output} > {params.outdir}/gtdbtk_combined_summary.tsv
        fi
        """

rule index_mags:
    input:
        os.path.join(
            "3_Outputs/drep/", "mag_catalogue.fa.gz"
        )
    output:
        os.path.join(
            "3_Outputs/drep/", "mag_catalogue.fa.gz.rev.2.bt2l"
        )
    wildcard_constraints:
        sample="(.+)_subsam_[^_]+"
    threads: 8
    resources:
        mem_gb="48",
        time='00:45:00'
    conda:
        "conda_envs/2_Assembly_Binning.yaml"
    shell:
        """
        # Index MAG catalogue with Bowtie2
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {input} {input}

        """

rule mag_catalogue_mapping:
    input:
        mags=os.path.join(
            "3_Outputs/drep/", "mag_catalogue.fa.gz"
        ),
        index=os.path.join(
            "3_Outputs/drep/", "mag_catalogue.fa.gz.rev.2.bt2l"
        ),
        r1=os.path.join(
            "../ind_assembly/2_Reads/4_Host_removed/", "{sample}_1.fq.gz"
        ),
        r2=os.path.join(
            "../ind_assembly/2_Reads/4_Host_removed/", "{sample}_2.fq.gz"
        ),
        gtdbtk="3_Outputs/gtdbtk/classify/gtdbtk.bac120.summary.tsv"
    output:
        os.path.join(
            "3_Outputs/mag_catalogue_bams/", "{sample}.bam"
            )
    wildcard_constraints:
        sample="(.+)_subsam_[^_]+"
    threads: 8
    resources:
        mem_gb="48",
        time='01:00:00'
    conda:
        "conda_envs/2_Assembly_Binning.yaml"
    shell:
        """
        # Map reads to assembly using Bowtie2
        bowtie2 \
            --time \
            --threads {threads} \
            -x {input.mags} \
            -1 {input.r1} \
            -2 {input.r2} \
        | samtools sort -@ {threads} -o {output}
        """

rule coverm_mags:
    input:
        lambda wildcards: expand(
            os.path.join(
                "3_Outputs/mag_catalogue_bams/", "{sample}.bam",
            ),
            sample=SAM
        )
    output:
        count_table=os.path.join(
            "3_Outputs/coverm/", "count_table.tsv"
        ),
    wildcard_constraints:
        sample="(.+)_subsam_[^_]+"
    threads: 8
    resources:
        mem_gb="48",
        time='04:00:00'
    conda:
        "conda_envs/2_Assembly_Binning.yaml"
    shell:
        """
        coverm genome \
            -b {input} \
            -s ^ \
            -m count covered_fraction length relative_abundance \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output.count_table}

        """