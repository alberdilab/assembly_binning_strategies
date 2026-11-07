configfile: "0_Code/configs/2_Assembly_Binning_config.yaml"

import pandas as pd

colnames = ['assembly', 'sample']
df = pd.read_csv("valid_combinations.tsv", sep="\t", names=colnames, header=None)

valid_combinations = set(
    (row["assembly"], row["sample"]) for _, row in df.iterrows()
)

#Define wildcards for later
ASSEM = list(set(df.iloc[0:, 0].tolist()))
SAM = list(set(df.iloc[0:, 1].tolist()))

# Function to filter combinations based on assembly
def filter_combinations(assembly):
    return [combo for combo in valid_combinations if combo[0] == assembly]

rule all:
    input:
        expand(
            os.path.join(
                "3_Outputs/refinement/{combo[0]}_refinement/",
                "{combo[0]}_metawrap_70_10_bins.stats",
            ),
            combo=valid_combinations,
        ),
        os.path.join(
            "3_Outputs/coverm/", "count_table.tsv"
        ),

rule mapping:
    input:
        contigs=os.path.join(
            "../ind_assembly/3_Outputs/2_Assemblies/", "{assembly}"
            ),
        r1=os.path.join(
            "../ind_assembly/2_Reads/4_Host_removed/", "{sample}_1.fq.gz"
            ),
        r2=os.path.join(
            "../ind_assembly/2_Reads/4_Host_removed/", "{sample}_2.fq.gz"
            )
    output:
        temp(os.path.join(
            "3_Outputs/bams/", "{assembly}/{assembly}_{sample}.bam"
            ))
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
        # Map reads to assembly using Bowtie2
        bowtie2 \
            --time \
            --threads {threads} \
            -x {input.contigs} \
            -1 {input.r1} \
            -2 {input.r2} \
        | samtools sort -@ {threads} -o {output}
        """

rule binning:
    input:
        bam=lambda wildcards: expand(
            os.path.join(
                "3_Outputs/bams/", "{assembly}/{assembly}_{sample}.bam"
            ),
            assembly=wildcards.assembly,
            sample=[combo[1] for combo in filter_combinations(wildcards.assembly)],
        ),
        contigs=os.path.join(
            "../ind_assembly/3_Outputs/2_Assemblies/", "{assembly}"
        )
    output:
        os.path.join(
            "3_Outputs/binning/", "{assembly}/{assembly}_binning_complete"
        )
    params:
        outdir=os.path.join("3_Outputs/binning" + "/{assembly}_binning")
    threads: 16
    resources:
        mem_gb="72",
        time='08:00:00'
    conda:
        "conda_envs/2_MetaWRAP.yaml"
    shell:
        """
        rm -rf {params.outdir}
        # Create dummy fq/assembly files to trick metaWRAP into running without mapping
        mkdir -p {params.outdir}
        mkdir -p {params.outdir}/work_files

        touch {params.outdir}/work_files/assembly.fa.bwt

        for bam in {input.bam}; do echo "@" > {params.outdir}/work_files/$(basename ${{bam/.bam/_1.fastq}}); done
        for bam in {input.bam}; do echo "@" > {params.outdir}/work_files/$(basename ${{bam/.bam/_2.fastq}}); done

        #Symlink BAMs for metaWRAP
        for bam in {input.bam}; do ln -sf `pwd`/$bam {params.outdir}/work_files/$(basename $bam); done

        # Run metaWRAP binning
        metawrap binning -o {params.outdir} \
            -t {threads} \
            -m {resources.mem_gb} \
            -a {input.contigs} \
            -l 1500 \
            --metabat2 \
            --maxbin2 \
            --concoct \
        {params.outdir}/work_files/*_1.fastq {params.outdir}/work_files/*_2.fastq

        # Create output for the next rule
        touch {output}
        """

rule refinement:
    input:
        os.path.join(
            "3_Outputs/binning/", "{assembly}/{assembly}_binning_complete"
        )
    output:
        stats=os.path.join(
            "3_Outputs/refinement/{assembly}_refinement/",
            "{assembly}_metawrap_70_10_bins.stats",
        ),
        contigmap=os.path.join(
            "3_Outputs/refinement/{assembly}_refinement/",
            "{assembly}_metawrap_70_10_bins.contigs",
        )
    params:
        binning=os.path.join("3_Outputs/binning" + "/{assembly}_binning"),
        outdir=os.path.join("3_Outputs/refinement" + "/{assembly}_refinement"),
        stats_dir=directory("3_Outputs/stats/{assembly}_stats/")
    threads: 6
    resources:
        mem_gb="172",
        time='06:00:00'
    conda:
        "conda_envs/2_MetaWRAP.yaml"
    shell:
        """
        # Setup checkM path
        export checkmdb={config[checkmdb]}
        printf $checkmdb | checkm data setRoot

        rm -rf {params.outdir}

        metawrap bin_refinement \
            -m {resources.mem_gb} \
            -t {threads} \
            -o {params.outdir} \
            -A {params.binning}/concoct_bins/ \
            -B {params.binning}/maxbin2_bins/ \
            -C {params.binning}/metabat2_bins/ \
            -c 70 \
            -x 10

        # Rename output files, and sort metawrap by bin name
        head -1 {params.outdir}/metawrap_70_10_bins.stats > {params.outdir}/mw_colnames.tsv
        sed '1d;' {params.outdir}/metawrap_70_10_bins.stats | sort -k1,1 -t$'\t' > {params.outdir}/mw_sorted.tsv
        cat {params.outdir}/mw_colnames.tsv {params.outdir}/mw_sorted.tsv > {params.outdir}/mw_sorted_col.tsv
        mv {params.outdir}/mw_sorted_col.tsv {output.stats}
        mv {params.outdir}/metawrap_70_10_bins.contigs {output.contigmap}
        sed -i'' '2,$s/bin/{wildcards.assembly}_bin/g' {output.stats}
        sed -i'' 's/bin/{wildcards.assembly}_bin/g' {output.contigmap}

        # Rename metawrap bins to match coassembly group:
        for bin in {params.outdir}/metawrap_70_10_bins/*.fa;
            do mv $bin ${{bin/bin./{wildcards.assembly}_bin.}};
        done

        # Compress output bins
        pigz -p {threads} {params.outdir}/metawrap_70_10_bins/*.fa

        #Print the number of MAGs to a file for combining with the assembly report
        mkdir -p {params.stats_dir}
        ls -l {params.outdir}/metawrap_70_10_bins/*.fa.gz | wc -l > {params.stats_dir}/{wildcards.assembly}_bins.tsv

        """

################################################################################
### dRep
rule dRep:
    input:
        lambda wildcards: expand(
            os.path.join(
                "3_Outputs/refinement/{assembly}_refinement/",
                "{assembly}_metawrap_70_10_bins.stats",
            ),
            assembly=ASSEM
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
        refinedir="3_Outputs/refinement"
    # conda:
    #     "conda_envs/drep.yaml"
    shell:
        """
        source activate /projects/mjolnir1/people/ncl550/0_software/e4fd3b3a5e8607bdeda760e6c6c0149b_
        module load bbmap/39.01
        # Copy MAGs into a shared folder
        mkdir -p 3_Outputs/drep/ && mkdir -p 3_Outputs/drep/All_metawrap_70_10_bins
        cp {params.refinedir}/*/metawrap_70_10_bins/*.fa.gz 3_Outputs/drep/All_metawrap_70_10_bins/

        # Parse/collate metawrap stats files for compatibility with dRep genomeinfo:
        echo -e "genome,completeness,contamination" > 3_Outputs/drep/header.txt
        for i in {params.refinedir}/*/C*_metawrap_70_10_bins.stats;
            do sed '1d;' $i | cut -f 1,2,3 --output-delimiter ',' >> 3_Outputs/drep/bin_info.txt;
                done
        sed -i'' 's@^@3_Outputs/drep/All_metawrap_70_10_bins/@g' 3_Outputs/drep/bin_info.txt
        sed -i'' 's/,/.fa.gz,/' 3_Outputs/drep/bin_info.txt
        cat 3_Outputs/drep/header.txt 3_Outputs/drep/bin_info.txt > 3_Outputs/drep/genome_info.csv

        # run dRep
        dRep dereplicate \
            3_Outputs/drep/98 \
            -p {threads} \
            -sa 0.98 \
            --genomeInfo 3_Outputs/drep/genome_info.csv \
            -g 3_Outputs/drep/All_metawrap_70_10_bins/*.fa.gz

        # Reformat MAG headers for CoverM
        for mag in 3_Outputs/drep/98/dereplicated_genomes/*.fa.gz;
            do rename.sh \
                in=$mag \
                out=3_Outputs/drep/98/dereplicated_genomes/$(basename ${{mag/.fa.gz/_renamed.fa.gz}}) \
                zl=9 \
                prefix=$(basename ${{mag/.fa.gz/^}});
        done

        # Concatentate the renamed, dereplicated MAGs
        cat 3_Outputs/drep/98/dereplicated_genomes/*_renamed.fa.gz > {output}

        # clean up
        rm 3_Outputs/drep/98/dereplicated_genomes/*_renamed.fa.gz

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
        --extension "renamed.fa.gz" \
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
        temp(os.path.join(
            "3_Outputs/mag_catalogue_bams/", "{sample}.bam"
            )
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