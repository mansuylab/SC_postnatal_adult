##############################
## @Author : Deepak Tanwar  ##
## @Date :   20210624       ##
##############################

## MultiQC report of RNA-Seq data

all: task1
task1:
	multiqc --version > log.stdout && multiqc -d -f -i "MultiQC report for SC RNA-Seq data</b>" -b "<i>--by</i> <br> <b>Deepak Tanwar</b>" -c multiqc_config.yaml -o ./output -z -v ./input  2> log.stderr >> log.stdout
