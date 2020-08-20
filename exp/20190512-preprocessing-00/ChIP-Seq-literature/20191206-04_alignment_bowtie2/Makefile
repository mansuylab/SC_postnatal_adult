############################
## Author: Deepak Tanwar  ##
## Date: 20191206	  ##
############################

# This makefile will run the Bowtie2 software to align the raw trimmed data.

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*_1_trimmed.fq.gz)
targets := $(patsubst $(source_dir)/%_1_trimmed.fq.gz, $(target_dir)/%.bam, $(files))

all: $(targets)

$(target_dir)/%.bam: $(source_dir)/%_1_trimmed.fq.gz
	bowtie2 --version 2> ./log/$(basename $(notdir $@)).log >> ./log/$(basename $(notdir $@)).log && bowtie2 -p 8 -x ./input/mm10/mm10 -q $< 2>> ./log/$(basename $(notdir $@)).log | perl -ne "print if((/XM:i:[0-3][^0-9]/) || (/^@/));" | samtools view -bS - > $(basename $@)_tmp.bam && samtools sort -@ 16 $(basename $@)_tmp.bam $(basename $@) && rm $(basename $@)_tmp.bam && samtools index $@
