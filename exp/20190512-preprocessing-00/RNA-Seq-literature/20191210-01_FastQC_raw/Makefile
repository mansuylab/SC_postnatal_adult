############################
## Author: Deepak Tanwar  ##
## Date: 20191210	  ##
############################

# This makefile will run the FastQC software to check the quality of FastQ files

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*_1.fastq.gz)
targets := $(patsubst $(source_dir)/%_1.fastq.gz, $(target_dir)/%.zip, $(files))

all: $(targets)

$(target_dir)/%.zip: $(source_dir)/%_1.fastq.gz
	fastqc -v > ./log/$(basename $(notdir $@)).log && fastqc $< -o ./output/ --noextract -f fastq -t 4 2>>./log/$(basename $(notdir $@)).log
