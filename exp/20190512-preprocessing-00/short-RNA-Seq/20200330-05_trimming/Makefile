############################
## Author: Deepak Tanwar  ##
## Date: 20200330	  ##
############################

# This makefile will run the TrimGalore software to improve the quality of FastQ files

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*_trimmedPolyA.fq.gz)
targets := $(patsubst $(source_dir)/%_trimmedPolyA.fq.gz, $(target_dir)/%_trimmed.fq.gz, $(files))

all: $(targets)

$(target_dir)/%_trimmed.fq.gz: $(source_dir)/%_trimmedPolyA.fq.gz
	trim_galore $< -q 30 --trim-n --length 16 --small_rna -j 8 -o $(target_dir) > ./log/$(basename $(basename $(notdir $@))).log 2>>./log/$(basename $(basename $(notdir $@))).log && mv ./output/$(basename $(basename $(notdir $<)))_trimmed.fq.gz $@ && mv ./output/$(notdir $<)_trimming_report.txt $(basename $(basename $@))_report.txt
#	cutadapt -m 15 -q 30 -u 3 -a AAAAAAAAAA -j 8 --trim-n -o $@ $< > ./log/$(basename $(basename $(notdir $@))).log 2>>./log/$(basename $(basename $(notdir $@))).log
