# Ensemble_calling_cnv

// Define processes for structural variant ensemble and annotation

process survclusterer_ensemble {
    cpus 4
    maxForks 4
    memory params.clu_mem
    publishDir "clusterer_vcf"

    input:
        tuple val(core), path(files)

    output:
        tuple path(clu_out), path("${clu_out}.tbi")

    script:
    clu_out = "${core}.clusterer.vcf.sv.gz"
    """
    #!/bin/bash
    set -euxo pipefail

    # Prepare file list
    [ -f file.txt ] && rm file.txt

    for file in ${files}; do
        [[ \$file == *.vcf.gz ]] && echo \$file >> file.txt
    done

    # Run Clusterer tool
    clusterer -d 1000 -t 4 file.txt ${params.genome} -o ${core}.clusterer.vcf

    # Compress and index output
    bgzip -k ${core}.clusterer.vcf.sv
    tabix ${clu_out}
    """
}

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












