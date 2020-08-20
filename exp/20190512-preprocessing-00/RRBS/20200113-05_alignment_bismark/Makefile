############################
## Author: Deepak Tanwar  ##
## Date: 20200113	  ##
############################

# This makefile will run `bismark` software to align trimmed FastQ files

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*.fq.gz)
targets := $(patsubst $(source_dir)/%.fq.gz, $(target_dir)/%_bismark_bt2.bam, $(files))

all: $(targets)

$(target_dir)/%_bismark_bt2.bam: $(source_dir)/%.fq.gz
	bismark --version > ./log/$(basename $(notdir $@)).log 2>> ./log/$(basename $(notdir $@)).log && bismark ./input/mm10 $< --multicore 8 -o ./output/ --temp_dir ./output/ --bowtie2 2>> ./log/$(basename $(notdir $@)).log >> ./log/$(basename $(notdir $@)).log
