#!/bin/bash

# directory structure:
# raw/batch/sample/fastq/sequencing_run/*.fastq.gz

## Trim data
find raw -mindepth 2 -maxdepth 2 -type d | parallel fastp --in1 {}/*/*/*R1*.fastq.gz --in2 {}/*/*/*R2*.fastq.gz --out1 clean/{/}_1.fq.gz --out2 clean/{/}_2.fq.gz --json clean/fastp_reports/{/}/fastp.json --html clean/fastp_reports/{/}/fastp.html -R {/}

## QC
find raw -mindepth 2 -maxdepth 2 -type d | parallel mkdir clean/fastp_reports/{/}
multiqc -d -f -o clean/fastp_reports clean/fastp_reports

## Quantify

for file in `find raw/ -maxdepth 2 -mindepth 2 -type d`
do
	mkdir quant quant/${file##*/}
done

find raw/ -maxdepth 2 -mindepth 2 -type d | parallel -j 4 \
	salmon quant -i genome_index/salmon/GRCm38.p6.gencode.vM22 \
	-l IU \
	-p 12 \
	-1 clean/{/}_1.fq.gz \
	-2 clean/{/}_2.fq.gz \
	--validateMappings \
	-o quant/{/}

multiqc quant -f -n mapping_stats -m salmon
