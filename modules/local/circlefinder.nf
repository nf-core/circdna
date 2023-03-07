process CIRCLEFINDER {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(split), path(concordant)

    output:
    tuple val(meta), path("*.microDNA-JT.txt")              , optional: true, emit: circdna
    tuple val(meta), path("*.circle_finder_exit_log.txt")   , optional: true, emit: log
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env bash

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
    # awk '\$1>3 {print \$2}' ${prefix}.concordant.id-freq.txt > ${prefix}.concordant.id-freqGr3.txt

    file_exists ${prefix}.concordant.id-freq3.txt

    grep -w -Ff ${prefix}.split.id-freq2.txt ${split} > ${prefix}.split_freq2.txt
    # grep -w -Ff ${prefix}.split.id-freq4.txt ${split} > ${prefix}.split_freq4.txt

    file_exists ${prefix}.split_freq2.txt

    #Selecting concordant pairs that were 1) mapped uniquely and 2) mapped on more than one loci (file "freqGr3.txt")
    grep -w -Ff ${prefix}.concordant.id-freq3.txt ${concordant} > ${prefix}.concordant_freq3.txt
#    grep -w -Ff ${prefix}.concordant.id-freqGr3.txt ${concordant} > ${prefix}.concordant_freqGr3.txt

    file_exists ${prefix}.concordant_freq3.txt

    #Step 7: Putting split read with same id in one line
    sed 'N;s/\\n/\\t/' ${prefix}.split_freq2.txt > ${prefix}.split_freq2.oneline.txt

    file_exists ${prefix}.split_freq2.oneline.txt

    #Step 8: Split reads map on same chromosome and map on same strand. Finally extracting id (split read same chromosome, split read same strand), collecting all the split reads that had quality >0
    awk '\$1==\$10 && \$7==\$16 && \$6>0 && \$15>0 {print \$4} ' ${prefix}.split_freq2.oneline.txt > \
        ${prefix}.split_freq2.oneline.S-R-S-CHR-S-ST.ID.txt

    file_exists ${prefix}.split_freq2.oneline.S-R-S-CHR-S-ST.ID.txt

    #Step 9: Based on unique id I am extracting one continuously mapped reads and their partner mapped as split read (3 lines for each id)
    grep -w -Ff "${prefix}.split_freq2.oneline.S-R-S-CHR-S-ST.ID.txt" "${prefix}.concordant_freq3.txt" > \
        "${prefix}.concordant_freq3.2SPLIT-1M.txt"

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

    file_exists ${prefix}.concordant_freq3.2SPLIT-1M.inoneline.txt

    #Step 11: Unique number of microDNA with number of split reads
    awk '\$1==\$11 && \$1==\$21 && \$7==\$17 && length(\$8)<=12 && length(\$18)<=12 && length(\$28)<=12'  ${prefix}.concordant_freq3.2SPLIT-1M.inoneline.txt | \
        awk '(\$7=="+" && \$27=="-") || (\$7=="-" && \$27=="+")' | \
        awk '{if (\$17=="+" && \$19=="second" && \$12<\$2 && \$22>=\$12 && \$23<=\$3) {printf ("%s\\t%d\\t%d\\n",\$1,\$12,\$3)} \
            else if (\$7=="+" && \$9=="second" && \$2<\$12 && \$22>=\$2 && \$23<=\$13) {printf ("%s\\t%d\\t%d\\n",\$1,\$2,\$13)} \
            else if (\$17=="-" && \$19=="second" && \$12<\$2 && \$22>=\$12 && \$23<=\$3) {printf ("%s\\t%d\\t%d\\n",\$1,\$12,\$3)} \
        else if (\$7=="-" && \$9=="second" && \$2<\$12 && \$22>=\$2 && \$23<=\$13) {printf ("%s\\t%d\\t%d\\n",\$1,\$2,\$13)} }' | \
    sort | uniq -c | awk '{printf ("%s\\t%d\\t%d\\t%d\\n",\$2,\$3,\$4,\$1)}' > ${prefix}.microDNA-JT.txt

    awk <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk --version) | grep "GNU Awk" | sed 's/^GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
}
