############################
## Author: Deepak Tanwar  ##
## Date: 20191218	  ##
############################

# This makefile will run `bismark` software to align trimmed FastQ files

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*_1_val_1.fq.gz)
targets := $(patsubst $(source_dir)/%_1_val_1.fq.gz, $(target_dir)/%_paired_bismark_bt2_pe.bam, $(files))

all: $(targets)

$(target_dir)/%_paired_bismark_bt2_pe.bam: $(source_dir)/%_1_val_1.fq.gz
	bismark --version > ./log/$(basename $(notdir $@)).log 2>> ./log/$(basename $(notdir $@)).log && bismark ./input/mm10 -1 $< -2 $(subst _1_val_1.fq.gz,_2_val_2.fq.gz, $<) --multicore 8 -o ./output/ --temp_dir ./output/ --bowtie2 2>> ./log/$(basename $(notdir $@)).log >> ./log/$(basename $(notdir $@)).log
#	samtools index $@.bam
