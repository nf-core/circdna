// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/software
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join

// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided as a string i.e. "options.args"
//               where "params.options" is a Groovy Map that MUST be provided via the addParams section of the including workflow.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

params.options = [:]
options        = initOptions(params.options)

process CIRCLEFINDER {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
//    conda (params.enable_conda ? "YOUR-TOOL-HERE" : null)
//    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
//        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
//    } else {
//        container "quay.io/biocontainers/YOUR-TOOL-HERE"
//    }

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/software/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(split), path(concordant)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    path "*.microDNA-JT.txt" optional true
    path "*.circle_finder_exit_log.txt" optional true
    // test
    path "*"
    // TODO nf-core: List additional required output channels/values here

    // path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/software/homer/annotatepeaks/main.nf
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "\\$options.args" variable
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads \\$task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    #!/usr/bin/env bash

    # prefix=${prefix}
    # split=${split}
    # concordant=${concordant}
    
    # Function to output an error if files do not exist
    file_exists () {
        [[ ! -s \$1 ]] && \
        echo "
    ERROR - CIRCLE_FINDER - $prefix
    ===================================
    Error \$1 does not exist or is empty.
    Stopped circle_finder.
    No circular DNA was identified.
    ===================================
    " > ${prefix}.circle_finder_exit_log.txt && exit
    }
    
    awk '{print \$4}' ${split} | sort | uniq -c > ${prefix}.split.id-freq.txt
    #This file "${prefix}.split.id-freq.txt" will be used for collecting split id that have frequency equal to 4.
    awk '\$1=="2" {print \$2}' ${prefix}.split.id-freq.txt > ${prefix}.split.id-freq2.txt
    # awk '\$1=="4" {print \$2}' ${prefix}.split.id-freq.txt > ${prefix}.split.id-freq4.txt
    
    awk '{print \$4}' ${concordant} | sort | uniq -c > ${prefix}.concordant.id-freq.txt
    #The following command will chose (may not be always true) one concordant and 2 split read
    
    awk '\$1=="3" {print \$2}' ${prefix}.concordant.id-freq.txt > ${prefix}.concordant.id-freq3.txt
    awk '\$1>3 {print \$2}' ${prefix}.concordant.id-freq.txt > ${prefix}.concordant.id-freqGr3.txt
    
    # check if output files exists and are not empty
    file_exists ${prefix}.concordant.id-freq3.txt
    file_exists ${prefix}.concordant.id-freqGr3.txt
    
    grep -w -Ff ${prefix}.split.id-freq2.txt ${split} > ${prefix}.split_freq2.txt
    # grep -w -Ff ${prefix}.split.id-freq4.txt ${split} > ${prefix}.split_freq4.txt
    
    # check if output files exist and are not empty
    file_exists ${prefix}.split_freq2.txt
    # file_exists ${prefix}.split_freq4.txt
    
    #Selecting concordant pairs that were 1) mapped uniquely and 2) mapped on more than one loci (file "freqGr3.txt")
    grep -w -Ff ${prefix}.concordant.id-freq3.txt ${concordant} > ${prefix}.concordant_freq3.txt
    grep -w -Ff ${prefix}.concordant.id-freqGr3.txt ${concordant} > ${prefix}.concordant_freqGr3.txt
    
    # check if output files exist and are not empty
    file_exists ${prefix}.concordant_freq3.txt
    file_exists ${prefix}.concordant_freqGr3.txt
    
    #Step 7: Putting split read with same id in one line
    sed 'N;s/\\n/\\t/' ${prefix}.split_freq2.txt > ${prefix}.split_freq2.oneline.txt
    # sed 'N;s/\\n/\\t/' ${prefix}.split_freq4.txt > ${prefix}.split_freq4.oneline.txt
    
    # check if output files exist and are not empty
    file_exists ${prefix}.split_freq2.oneline.txt
    # file_exists ${prefix}.split_freq4.oneline.txt
    
    #Step 8: Split reads map on same chromosome and map on same strand. Finally extracting id (split read same chromosome, split read same strand), collecting all the split reads that had quality >0
    awk '\$1==\$10 && \$7==\$16 && \$6>0 && \$15>0 {print \$4} ' ${prefix}.split_freq2.oneline.txt > \
        ${prefix}.split_freq2.oneline.S-R-S-CHR-S-ST.ID.txt
    
    # check if output files exist and are not empty
    file_exists ${prefix}.split_freq2.oneline.S-R-S-CHR-S-ST.ID.txt
    
    #Step 9: Based on unique id I am extracting one continuously mapped reads and their partner mapped as split read (3 lines for each id)
    grep -w -Ff "${prefix}.split_freq2.oneline.S-R-S-CHR-S-ST.ID.txt" "${prefix}.concordant_freq3.txt" > \
        "${prefix}.concordant_freq3.2SPLIT-1M.txt"
    
    # check if output files exist and are not empty
    file_exists ${prefix}.concordant_freq3.2SPLIT-1M.txt
    
    #Step 10: Sorting based on read-id followed by length of mapped reads.
    awk 'BEGIN{FS=OFS="\\t"} {gsub("M", " M ", \$8)} 1' ${prefix}.concordant_freq3.2SPLIT-1M.txt | \
        awk 'BEGIN{FS=OFS="\\t"} {gsub("S", " S ", \$8)} 1' | \
        awk 'BEGIN{FS=OFS="\\t"} {gsub("H", " H ", \$8)} 1' | \
        awk 'BEGIN{FS=OFS=" "} {if ((\$9=="M" && \$NF=="H") || \
            (\$9=="M" && \$NF=="S"))  {printf ("%s\\tfirst\\n",\$0)} \
            else if ((\$9=="S" && \$NF=="M") || (\$9=="H" && \$NF=="M")) {printf ("%s\\tsecond\\n",\$0)} \
            else  {printf ("%s\\tconfusing\\n",\$0)}}' | \
        awk 'BEGIN{FS=OFS="\\t"} {gsub(" ", "", \$8)} 1' | \
        awk '{printf ("%s\\t%d\\n",\$0,(\$3-\$2)+1)}' | \
        sort -k4,4 -k10,10n | sed 'N;N;s/\\n/\\t/g' | \
        awk '{if (\$5==\$15) {print \$0}  \
            else if ((\$5=="1" && \$15=="2" && \$25=="1") || (\$5=="2" && \$15=="1" && \$25=="2")) \
                {printf ("%s\\t%d\\t%d\\t%s\\t%d\\t%d\\t%s\\t%s\\t%s\\t%d\\t%s\\t%d\\t%d\\t%s\\t%d\\t%d\\t%s\\t%s\\t%s\\t%d\\t%s\\t%d\\t%d\\t%s\\t%d\\t%d\\t%s\\t%s\\t%s\\t%d\\n", \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$21,\$22,\$23,\$24,\$25,\$26,\$27,\$28,\$29,\$30,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20)} \
            else if ((\$5=="1" && \$15=="2" && \$25=="2") || (\$5=="2" && \$15=="1" && \$25=="1")) \
            {printf ("%s\\t%d\\t%d\\t%s\\t%d\\t%d\\t%s\\t%s\\t%s\\t%d\\t%s\\t%d\\t%d\\t%s\\t%d\\t%d\\t%s\\t%s\\t%s\\t%d\\t%s\\t%d\\t%d\\t%s\\t%d\\t%d\\t%s\\t%s\\t%s\\t%d\\n", \$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$21,\$22,\$23,\$24,\$25,\$26,\$27,\$28,\$29,\$30,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10)} }' \
        > ${prefix}.concordant_freq3.2SPLIT-1M.inoneline.txt
    
    # check if output files exist and are not empty
    file_exists ${prefix}.concordant_freq3.2SPLIT-1M.inoneline.txt
    
    #Step 11: Unique number of microDNA with number of split reads
    awk '\$1==\$11 && \$1==\$21 && \$7==\$17 && length(\$8)<=12 && length(\$18)<=12 && length(\$28)<=12'  ${prefix}.concordant_freq3.2SPLIT-1M.inoneline.txt | \
        awk '(\$7=="+" && \$27=="-") || (\$7=="-" && \$27=="+")' | \
        awk '{if (\$17=="+" && \$19=="second" && \$12<\$2 && \$22>=\$12 && \$23<=\$3) {printf ("%s\\t%d\\t%d\\n",\$1,\$12,\$3)} \
            else if (\$7=="+" && \$9=="second" && \$2<\$12 && \$22>=\$2 && \$23<=\$13) {printf ("%s\\t%d\\t%d\\n",\$1,\$2,\$13)} \
            else if (\$17=="-" && \$19=="second" && \$12<\$2 && \$22>=\$12 && \$23<=\$3) {printf ("%s\\t%d\\t%d\\n",\$1,\$12,\$3)} \
        else if (\$7=="-" && \$9=="second" && \$2<\$12 && \$22>=\$2 && \$23<=\$13) {printf ("%s\\t%d\\t%d\\n",\$1,\$2,\$13)} }' | \
    sort | uniq -c | awk '{printf ("%s\\t%d\\t%d\\t%d\\n",\$2,\$3,\$4,\$1)}' > ${prefix}.microDNA-JT.txt  
    
    """
}
