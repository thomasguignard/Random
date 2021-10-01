task computePoorCoverage {
	String SampleID
	String OutDirSampleID = ""
	String OutDir
	String WorkflowType
	String GenomeVersion
	String BedToolsExe
	String AwkExe
	String SortExe
 	#task specific variables
	File IntervalBedFile
	Int BedtoolsLowCoverage
	Int BedToolsSmallInterval
	File BamFile
	#runtime attributes
	Int Cpu
	Int Memory
	String OutputDirSampleID = if OutDirSampleID == "" then SampleID else OutDirSampleID

	String poorCoverageDir
	File TsvCoverageFile

	command <<<
		${BedToolsExe} genomecov -ibam ${BamFile} -bga \
		| ${AwkExe} -v low_coverage="${BedtoolsLowCoverage}" '$4<low_coverage' \
		| ${BedToolsExe} intersect -a ${IntervalBedFile} -b - \
		| ${SortExe} -k1,1 -k2,2n -k3,3n \
		| ${BedToolsExe} merge -c 4,8,8 -o distinct,min,max -i - \
		| ${BedToolsExe} intersect -loj -a -  -b ${poorCoverageDir}/*tsv  \
		| ${BedToolsExe} intersect -a -  -b ${TsvCoverageFile} -wb  \
		| ${AwkExe} -v small_intervall="${BedToolsSmallInterval}" \
		'BEGIN {OFS="\t";print "#chr","start","end","region","size bp","type","Cov-min","Cov_max","Occurrence","Target_mean_coverage","gene","UCSC link"} {split($4,gene,":");a=($3-$2+1);if(a<small_intervall) {b="SMALL_INTERVAL"} else {b="OTHER"};url="http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db='${GenomeVersion}'&position="$1":"$2-10"-"$3+10"&highlight='${GenomeVersion}'."$1":"$2"-"$3;print $1, $2, $3, $4, a, b, $5, $6, $7, $8, $9, url}' \
		> "${OutDir}${OutputDirSampleID}/${WorkflowType}/coverage/${SampleID}_poor_coverage.tsv"
	>>>
	output {
		File poorCoverageFile = "${OutDir}${OutputDirSampleID}/${WorkflowType}/coverage/${SampleID}_poor_coverage.tsv"
	}
	runtime {
		cpu: "${Cpu}"
		requested_memory_mb_per_core: "${Memory}"
	}
}
