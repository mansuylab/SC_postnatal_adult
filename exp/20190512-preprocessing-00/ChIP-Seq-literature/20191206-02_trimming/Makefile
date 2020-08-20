############################
## Author: Deepak Tanwar  ##
## Date: 20191206	  ##
############################

# This makefile will run the TrimGalore software to improve the quality of FastQ files

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*_1.fastq.gz)
targets := $(patsubst $(source_dir)/%_1.fastq.gz, $(target_dir)/%_trimmed.fq.gz, $(files))

all: $(targets)

$(target_dir)/%_trimmed.fq.gz: $(source_dir)/%_1.fastq.gz
	trim_galore $< -q 30 --length 30 -o $(target_dir) > ./log/$(basename $(basename $(notdir $@))).log 2>>./log/$(basename $(basename $(notdir $@))).log
