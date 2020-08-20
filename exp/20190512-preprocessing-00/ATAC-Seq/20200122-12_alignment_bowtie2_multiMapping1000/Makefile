############################
## Author: Deepak Tanwar  ##
## Date: 20200122	  ##
############################

# This makefile will run the Bowtie2 software to align the raw trimmed data.

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*_R1_val_1.fq.gz)
targets := $(patsubst $(source_dir)/%_R1_val_1.fq.gz, $(target_dir)/%.bam, $(files))

all: $(targets)

$(target_dir)/%.bam: $(source_dir)/%_R1_val_1.fq.gz
	bowtie2 --version 2> ./log/$(basename $(notdir $@)).log >> ./log/$(basename $(notdir $@)).log && bowtie2 -X 2000 --end-to-end --no-mixed --no-discordant --score-min L,-0.4,-0.4 -p 8 -k 1000 -x ./input/mm10/mm10 -1 $< -2 input/$(basename $(notdir $@))_R2_val_2.fq.gz 2>> ./log/$(basename $(notdir $@)).log | samtools view -Sb -@ 8 - -o $(basename $@)_tmp.bam && samtools sort -@ 8 $(basename $@)_tmp.bam -o $@ && samtools index $@ && rm $(basename $@)_tmp.bam
