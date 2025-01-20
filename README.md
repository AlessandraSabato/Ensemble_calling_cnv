# Ensemble_calling_cnv
Here's an overview of how to get started with a Nextflow workflow that implements an ensemble calling, annotation, and visualization pipeline. This pipeline will process genomic data, including the steps of variant calling, merging multiple callsets, annotating the variants, and visualizing the results.

## 1.1 Install Nextflow and Prerequisites

Before starting, ensure you have **Nextflow** and all the necessary tools installed. You will also need to set up a computational environment that supports Nextflow execution.

## 1.2 Configuring the nextflow.config File

The **nextflow.config** file is crucial for setting up parameters, resource allocations, and execution profiles. Below is an outline of the essential configurations and their meanings.

# 2. Configuration of nextflow.config

This file allows you to define paths, memory allocations, environment variables, and execution profiles for running the workflow.

## 2.1 Configuration Parameters

In the `params` section of the **nextflow.config** file, you will define the following:

params {
    // Path to the reference genome file (FASTA format)
    genome = "/path/to/reference/genome.fa"  
    
    // Path to the directory containing BAM files
    input = "/path/to/bam_files"  
    
    // List of directories containing .VCF.GZ and .VCF.GZ.TBI files
    directories = ["/path/to/vcf_directory1", "/path/to/vcf_directory2"]  
    
    // Pattern for selecting BAM files from the input directory
    samples = "*.bam*"  
    
    // Pattern for selecting VCF files from the directories
    ensamples = "*.vcf*"  
    
    // Memory configuration for each tool
    ins_mem = '24GB'  // Memory for INSurVeyor
    sva_mem = '24GB'  // Memory for SvABA
    clu_mem = '16GB'  // Memory for SurVClusterer
    sur_mem = '24GB'  // Memory for SURVIVOR
    tru_mem = '32GB'  // Memory for Truvari
    annotsv_mem = '20GB'  // Memory for AnnotSV annotation
    knotsv_mem = '20GB'  // Memory for KnotSV annotation (HTML and XLS)
    vep_mem = '20GB'  // Memory for VEP annotation
}
2.1.1 Memory Parameters:

    The memory (sur_mem, tru_mem, etc.) for each tool is adjusted to ensure the processes have enough resources based on your input files and the complexity of the tasks.
    Added annotsv_mem, knotsv_mem, and vep_mem for memory requirements for the AnnotSV, KnotSV, and VEP annotation processes, respectively.

2.1.2 Input and Output File Patterns:

    The input and samples parameters are set for the BAM files to match your workflow setup.
    The directories and ensamples parameters are for selecting VCF files from different directories, ensuring flexibility in file selection.

## 2.2 Environment Variables and Execution Profiles

You will define environment variables, such as paths to libraries and executables, and specify execution profiles for different environments (SLURM or Local).

env {
    LIBRARY_PATH = "/path/to/library"  
    PATH = "/path/to/bin:$PATH"  
    JAVA_HOME = "/path/to/java"  
}

profiles {
    // SLURM profile (for cluster execution)
    slurm {
        executor.name = 'slurm'
        queue = 'batch'
        memory = '32GB'
        cpus = 4
        time = '24:00:00'
    }

    // Local profile (for local execution)
    local {
        executor.name = 'local'
        memory = '16GB'
        cpus = 2
    }
}

# 3. Workflow Execution

Once you’ve configured your **nextflow.config**, you can begin executing your workflow.

## 3.1 Create Directories and Set Up Environment

Before running the workflow, set up the necessary directories:

mkdir nextflow_running
cd nextflow running

## 3.2 Running the Workflow on SLURM Cluster

To run the workflow on a **SLURM**-based cluster, use the following command:

nextflow run main.nf -profile slurm

## 3.3 Running the Workflow Locally

For local execution:

nextflow run main.nf -profile local

## 3.4 Resuming the Workflow

To resume the workflow from a previous point in case of failure:

nextflow run main.nf -profile slurm -resume

# 4. Workflow Processes and Parameters

Each process in the workflow is associated with specific tools and configurations for variant calling, merging, annotation, and visualization.

## 4.1 survclusterer_ensemble (SurVClusterer)

This process clusters structural variants using **SurVClusterer**, which groups similar variants into meaningful clusters.


clusterer -d 1000 -t 4 file.txt ${params.genome} -o output_file.vcf


Parameters:
-d 1000   # Clusters variants within 1000 base pairs.
-t 4      # Uses 4 threads for parallel processing.

## 4.2 survivor_ensemble (SURVIVOR)

**SURVIVOR** merges multiple VCF files into a single output and filters variants based on various parameters.

SURVIVOR merge file.txt 1000 1 1 -1 -1 -1 output_file.vcf

Parameters:
Merges VCF files, with variants within 1000 base pairs being merged.

Explanation: This command ensures that only variants with sufficient support are retained and combines multiple VCF files.

## 4.3 truvari_ensemble (Truvari)

**Truvari** collapses and merges structural variants into a concise VCF representation.

truvari collapse -i bcftools.merge.vcf.gz -o truvari.vcf.gz -c output_file.vcf.gz

Explanation: **Truvari** reduces redundancy and improves structural variant analysis by merging similar variants.

## 4.4 annotsv_annotation (AnnotSV)

**AnnotSV** annotates structural variants with functional information, helping researchers understand their biological relevance.

${ANNOTSV_DIR} -genomeBuild GRCh38 -SVinputFile truvari.vcf -outputFile annotsv -outputDir AnnotSV_output -SVinputInfo 1

Explanation: **AnnotSV** adds biological context, such as gene impact and disease associations, to the structural variants.

## 4.5 knotsv_annotation_html and knotsv_annotation_xslm (KnotSV)

**KnotSV** generates HTML and Excel reports with visualizations of annotated structural variants.

perl /home/ale_sab/knotAnnotSV/knotAnnotSV2XL.pl --configFile /home/ale_sab/knotAnnotSV/config_AnnotSV.yaml --annotSVfile ${annot_out} --outDir KnotSV_dir --outPrefix knotsv_output --genomeBuild GRCh38 --LOEUFcolorRange 1 --geneCountThreshold 40




Parameters:
--LOEUFcolorRange 1  # Specifies color range for Loss of Function (LOEUF) annotations.

Explanation: **KnotSV** outputs detailed visual reports that aid in the interpretation of structural variants.

## 4.6 vep_annotation (VEP)

**VEP** annotates VCF files with detailed gene information, including predicted functional consequences of variants.

singularity exec /home/ale_sab/project_cnv_vr/ensembl-vep_latest.sif vep -i cleaned_input.vcf --format vcf --output_file ${vep_annotated_vcf} --vcf --everything --assembly GRCh38 --symbol --canonical --vcf_info_field ANN --cache --offline --dir_cache ~/.vep --force_overwrite --cache_version 113 --sift b --polyphen b --biotype --hgvs --fasta



Explanation: **VEP** adds in-depth annotations, helping researchers understand the effects of variants on genes, transcripts, and proteins.

















process survivor_ensemble {
    cpus 2
    maxForks 4
    memory params.sur_mem
    publishDir "survivor_vcf"

    input:
        tuple val(core), path(files)

    output:
        tuple path(sur_out), path("${sur_out}.tbi")

    script:
    sur_out = "${core}.surv_rem.vcf.gz"
    """
    #!/bin/bash
    set -euxo pipefail

    # Prepare file list
    [ -f file.txt ] && rm file.txt

    for file in ${files}; do
        [[ \$file == *.vcf.gz ]] && zcat \$file > "\$file.vcf" && echo "\$file.vcf" >> file.txt
    done

    # Merge and process variants using SURVIVOR
    SURVIVOR merge file.txt 1000 1 1 -1 -1 -1 ${core}.survivor.vcf
    bgzip ${core}.survivor.vcf

    # Post-process the merged VCF
    zgrep -v HLA ${core}.survivor.vcf.gz | bgzip | bcftools sort -Oz -o ${core}.survivor_sorted.vcf.gz
    bcftools view -e 'INFO/END < POS' ${core}.survivor_sorted.vcf.gz -o ${sur_out}
    bcftools index -t ${sur_out}
    """
}

process truvari_ensemble {
    cpus 2
    maxForks 4
    memory params.tru_mem
    module 'truvari'
    publishDir "truvari_vcf"

    input:
        path(files)

    output:
        tuple path(tru_out), path("${tru_out}.tbi")

    script:
    tru_out = "tru_coll_sorted.vcf.gz"
    """
    #!/bin/bash
    set -euxo pipefail

    # Merge input VCF files
    vcf_files=""
    for file in ${files}; do
        [[ \$file == *.vcf.gz ]] && vcf_files="\$vcf_files \$file"
    done
    bcftools merge -m none \$vcf_files -Oz -o bcftools.merge.vcf.gz --force-samples
    bcftools index -t bcftools.merge.vcf.gz

    # Collapse variants using Truvari
    truvari collapse -i bcftools.merge.vcf.gz -o truvari.vcf.gz -c tru_collapse.vcf.gz

    # Sort, compress, and index the output
    vcf-sort tru_collapse.vcf.gz | bgzip -c > tru_coll_sorted.vcf.gz
    tabix -p vcf ${tru_out}
    """
}

process annotsv_annotation {
    cpus 2
    maxForks 4
    memory '20GB'
    module 'AnnotSV/3.4'
    publishDir "annotsv_annotations"

    input:
        path(truvari_vcf)

    output:
        tuple path(annotated_vcf), path("${annotated_vcf}.summary")

    script:
    annotated_vcf = "AnnotSV_output/annotsv_annotated.vcf"
    """
    #!/bin/bash
    set -euxo pipefail

    # Prepare AnnotSV input
    gunzip -k -c ${truvari_vcf} > truvari.vcf
    mkdir -p AnnotSV_output

    # Annotate structural variants
    ${ANNOTSV_DIR} -genomeBuild GRCh38 -SVinputFile truvari.vcf -outputFile annotsv -outputDir AnnotSV_output -SVinputInfo 1
    """
}
process knotsv_annotation_html {
    cpus 2
    maxForks 4
    memory '20GB'
    publishDir "knotsv_vcf"
    input:
        tuple path(annot_out), path(summary)
    output:
        tuple path(html_out), path("${html_out}.config")
    script:
    html_out = "KnotSV_dir/knotsv_output.html"
    """
    set -euxo pipefail

    mkdir -p KnotSV_dir
    perl /home/ale_sab/knotAnnotSV/knotAnnotSV.pl --configFile /home/ale_sab/knotAnnotSV/config_AnnotSV.yaml --annotSVfile ${annot_out} --outDir KnotSV_dir --outPrefix knotsv_output --genomeBuild GRCh38 --LOEUFcolorRange 1
    """
}

process knotsv_annotation_xslm {
    cpus 2
    maxForks 4
    memory '20GB'
    publishDir "knotsv_vcf"
    input:
        tuple path(annot_out), path(summary)
    output:
        tuple path(xls_out), path("${xls_out}.config")
    script:
    xls_out = "KnotSV_dir/knotsv_output.xlsx"
    """
    set -euxo pipefail

    mkdir -p KnotSV_dir
    perl /home/ale_sab/knotAnnotSV/knotAnnotSV2XL.pl --configFile /home/ale_sab/knotAnnotSV/config_AnnotSV.yaml --annotSVfile ${annot_out} --outDir KnotSV_dir --outPrefix knotsv_output --genomeBuild GRCh38 --LOEUFcolorRange 1 --geneCountThreshold 40
    """
}
process vep_annotation {
    cpus 2
    maxForks 4
    memory '20GB'
    publishDir "vep_annotated_vcf"

    input:
        path(truvari_vcf)

    output:
        tuple path(vep_annotated_vcf), path("${vep_annotated_vcf}.index")

    script:
    vep_annotated_vcf = "vep_annotated_variants.vcf"
    """
    #!/bin/bash
    set -euxo pipefail

    # Prepare VEP input
    gunzip -k -c ${truvari_vcf} > truvari_input.vcf
    bcftools view truvari_input.vcf | grep -v "chrUn" | grep -v "SVTYPE=NA" > cleaned_input.vcf

    # Annotate using VEP
    singularity exec /home/ale_sab/project_cnv_vr/ensembl-vep_latest.sif vep \
        -i cleaned_input.vcf --format vcf --output_file ${vep_annotated_vcf} --vcf --everything \
        --assembly GRCh38 --symbol --canonical --vcf_info_field ANN --cache --offline \
        --dir_cache ~/.vep --force_overwrite --cache_version 113 --sift b --polyphen b \
        --biotype --hgvs --fasta /home/ale_sab/data/genomes/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa
    """
}

workflow {
    // Input channels
    ins_res = Channel.fromFilePairs(params.input + params.samples) { f -> f.getSimpleName() }
    sva_res = Channel.fromFilePairs(params.input + params.samples) { f -> f.getSimpleName() }
    prev_vcf = Channel.fromPath(params.directories.collect { it + params.ensamples }).map { f -> [f.getSimpleName(), f] }

    // Ensemble processes
    all_res = ins_res.mix(sva_res).map { [it[0], it[1]] }
    clustered_vcf = all_res | survclusterer_ensemble
    survivor_vcf = clustered_vcf | survivor_ensemble | truvari_ensemble

    // Annotation processes
    annotated_vcf = survivor_vcf | annotsv_annotation
    vep_annotated = survivor_vcf | vep_annotation

    // Outputs
    clustered_vcf.view()
    annotated_vcf.view()
    vep_annotated.view()
}












