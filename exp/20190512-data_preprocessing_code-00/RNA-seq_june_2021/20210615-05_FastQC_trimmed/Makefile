############################
## Author: Deepak Tanwar  ##
## Date: 20210616	  ##
############################

# This makefile will run the FastQC software to check the quality of FastQ files

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*.fq.gz)
targets := $(patsubst $(source_dir)/%.fq.gz, $(target_dir)/%.zip, $(files))

all: $(targets)

$(target_dir)/%.zip: $(source_dir)/%.fq.gz
	fastqc -v > ./log/$(basename $(notdir $@)).stdout && fastqc $< -o ./output/ --noextract -f fastq -t 4 2>./log/$(basename $(notdir $@)).stderr
