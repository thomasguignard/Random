#! /bin/bash


echo "USAGE: nohup $0  folderName  refgenome(hg19 or hg38) &"

log="`date +%d-%b~%H-%M-%S`"_pipeline.log

touch $log

echo -e "################## Log file for minion pipeline ##############\n\nRun Initiated: `date +%F,%R:%S`" >> $log

if [ -z $1 ] || [ -z $2 ]; then
	echo -e "Missing arguments \n USAGE: nohup ./pipelineNV.sh folderName  refgenome(hg19 or hg38) &"
	echo -e "Missing arguments \n USAGE: nohup ./pipelineNV.sh folderName  refgenome(hg19 or hg38) &" >> $log
	echo "exit 1" >> $log
	echo "exit 1"
	exit 1 
fi

echo -e "nohup ./pipelineNV.sh $1 $2 &" >> $log
echo -e "Analyzed data: $1 \nReference genome used: $2" >> $log

echo -e "\n##################### DIRECTORY EXISTENCE #######################\n">> $log
echo -e "\n##################### DIRECTORY EXISTENCE #######################\n"
ARG1=$PWD/$1

echo $ARG1

echo "Controling Directory existence" >> $log
echo "Controling Directory existence"
if [ ! -d "$ARG1" ]; then
  # Control will enter here if $DIRECTORY doesn't exist.
	echo "This folder doesn't exist, please check folder name.">>$log
	echo "This folder doesn't exist, please check folder name."
	echo "USAGE:  ./pipelineNV.sh folderName refgenome(hg19 or hg38)">>$log
	echo "exit 1" >> $logfile
	echo "exit 1"
	exit 1
fi

echo "Folder $1 exists." >> $log

echo -e "\n########################## SUBFOLDERS CREATION ##########################\n">> $log
echo -e "\n########################## SUBFOLDERS CREATION ##########################\n"
echo "Creating subfolders..."
echo "Creating subfolders..." >>  $log

#prepare-folder:

FOLDER=$(ls -d ${ARG1}*/*)

echo $FOLDER
mv $log $FOLDER/
logfile=$FOLDER/$log

#Create subfolders only if they do not exist

if [ ! -d ${FOLDER}/FASTQ/ ];then
	echo "Create ${FOLDER}/FASTQ/">>$logfile
	mkdir "${FOLDER}/FASTQ/" 2>> $logfile 
else
	echo "${FOLDER}/FASTQ/ already exists.">> $logfile
fi

if [ ! -d ${FOLDER}/FASTQ/ALL/ ];then
	echo "Create ${FOLDER}/FASTQ/ALL/">>$logfile
	mkdir "${FOLDER}/FASTQ/ALL/" 2>> $logfile
else
        echo "${FOLDER}/FASTQ/ALL already exists.">> $logfile
fi

if [ ! -d ${FOLDER}/BAM_$2/ ];then
	echo "Create ${FOLDER}/BAM_$2/">>$logfile
	mkdir "${FOLDER}/BAM_$2/" 2>> $logfile
else
	echo "${FOLDER}/BAM_$2/ already exists.">> $logfile
fi

if [ ! -d ${FOLDER}/NANOSV_$2/ ]; then
	echo "Create ${FOLDER}/NANOSV_$2/">>$logfile
	mkdir "${FOLDER}/NANOSV_$2/" 2>> $logfile
else
	echo "${FOLDER}/NANOSV_$2 already exists.">> $logfile
fi

if [ ! -d ${FOLDER}/NANOVAR_$2/ ]; then
	echo "Create ${FOLDER}/NANOVAR_$2/">>$logfile
	mkdir "${FOLDER}/NANOVAR_$2/" 2>> $logfile
else
	echo "${FOLDER}/NANOVAR_$2/ already exists.">> $logfile
fi

if [ ! -d ${FOLDER}/NANOVAR_$2/nanovar2annotSV/ ]; then
	echo "Create ${FOLDER}/NANOVAR_$2/nanovar2annotSV/">>$logfile
	mkdir "${FOLDER}/NANOVAR_$2/nanovar2annotSV/" 2>> $logfile
else
	echo "${FOLDER}/NANOVAR_$2/nanovar2annotSV/ already exists." >> $logfile
fi

if [ ! -d ${FOLDER}/NANOVAR_$2/hsblast4igv/ ]; then
        echo "Create ${FOLDER}/NANOVAR_$2/hsblast4igv/">>$logfile
        mkdir "${FOLDER}/NANOVAR_$2/hsblast4igv/" 2>> $logfile
else
        echo "${FOLDER}/NANOVAR_$2/hsblast4igv/ already exists." >> $logfile
fi

echo "done"
echo -e "\n############################ BASECALLING ##########################\n" >> $logfile
echo -e "\n############################ BASECALLING ##########################\n"
#first-call:
#guppy_basecaller --input_path $FOLDER"/fast5/"  --save_path $FOLDER"/FASTQ/" --flowcell FLO-MIN106 --kit SQK-LSK109 --cpu_threads_per_caller 1 --num_callers 11
if [ ! -d $FOLDER/fast5/ ]; then
	echo "WARNING: No fast5 files available. Please check the folder content." >> $logfile
	echo "WARNING: No fast5 files available. Please check the folder content."
	echo "exit 1"
	exit 1
elif [ ! -e $FOLDER/FASTQ/ALL/ALL.fastq ]; then
	echo "`date +%F,%R:%S`: launch basecalling..." >> $logfile
	echo "launch basecalling..."
	echo -e "\nguppy_basecaller --input_path $FOLDER/fast5/  --save_path $FOLDER/FASTQ/  --cpu_threads_per_caller 1 --num_callers 11 -m /opt/ont/guppy/data/template_r9.4.1_450bps_fast.jsn   -c /opt/ont/guppy/data/dna_r9.4.1_450bps_fast.cfg" >> $logfile
	guppy_basecaller --input_path $FOLDER"/fast5/"  --save_path $FOLDER"/FASTQ/"  --cpu_threads_per_caller 1 --num_callers 11 -m /opt/ont/guppy/data/template_r9.4.1_450bps_fast.jsn   -c /opt/ont/guppy/data/dna_r9.4.1_450bps_fast.cfg 2>> $logfile
	if [ $? -ne 0 ];then
        	#STOP the script if a problem happens in basecalling
		echo -e "\nWARNING: Error in basecalling. The script has been stopped." >> $logfile
        	echo  "WARNING: Error in basecalling. The script has been stopped."
		echo "exit 1" >> $logfile
		echo "exit 1"
		exit 1
	else
		echo -e "\n`date +%F,%R:%S`: Basecalling done" >> $logfile
		echo "Basecalling done"
	fi
else
	echo "\nALL.fastq file already exists." >> $logfile
fi

echo -e "\n######################## FASTQ CONCATENATION ######################\n" >> $logfile
echo -e "\n######################## FASTQ CONCATENATION ######################\n"

if [ ! -e $FOLDER/FASTQ/*_0_0.fastq ];then
	echo "WARNING: No fastq files available for fastq concatenation. Please check basecalling.">> $logfile
	echo "WARNING: No fastq files available for fastq concatenation. Please check basecalling."
elif [ ! -e $FOLDER/FASTQ/ALL/ALL.fastq ]; then
	echo -e "`date +%F,%R:%S`: Launch fastq concatenation in one" >> $logfile
	echo "Launch fastq concatenation in one"
	for i in $FOLDER"/FASTQ/*tq"; do cat $i >> $FOLDER"/FASTQ/ALL/ALL.fastq" ; done
	if [ $? -ne 0 ];then
		echo -e "\nWARNING: Error in FASTQ concatenation. Please check basecalling." >> $logfile
		echo "WARNING: Error in FASTQ concatenation. Please check basecalling."
	else
		echo -e "\n`date +%F,%R:%S`: Concatenation done." >> $logfile
		echo "Concatenation done."
	fi
else
	echo "ALL.fastq file already exists." >> $logfile
fi

echo -e "\n######################## NANOVAR CALLING ##########################\n" >> $logfile
echo -e "\n######################## NANOVAR CALLING ##########################\n"
#Nanovar calling with a confidence score threshold of 0.5 (i.e. -S 0.5)
	#Execute the script only when ALL.fastq exists and nanovar vcf files with same names do not exist.

if [ ! -e $FOLDER/FASTQ/ALL/ALL.fastq ]; then
	echo "WARNING: No ALL.fastq file available for NanoVar calling. Please check basecalling." >> $logfile
	echo "WARNING: No ALL.fastq file available for NanoVar calling. Please check basecalling."
elif [ -e $FOLDER/NANOVAR_$2/nanovar_results/ALL.ouput.total.vcf ]; then
	echo "Nanovar vcf files with same names already exist. Please check filenames." >> $logfile
	echo "Nanovar vcf files with same names already exist. Please check filenames."
else
	echo "Lauch Nanovar calling..."
	echo -e "`date +%F,%R:%S`: Lauch Nanovar calling..." >> $logfile
	echo "~/TOOLS/NanoVar/nanovar -l $FOLDER/FASTQ/ALL/ALL.fastq -o $FOLDER/NANOVAR_$2/ -t 11 -r ~/resources/$2/$2.fa -f $2 -S 0.5" >> $logfile
	~/TOOLS/NanoVar/nanovar -l $FOLDER"/FASTQ/ALL/ALL.fastq" -o $FOLDER"/NANOVAR_$2/" -t 11 -r ~/resources/$2/$2.fa -f $2 -S 0.5 2>> $logfile
	if [ $? -ne 0 ];then
		echo -e "\nWARNING: Error in Nanovar calling." >> $logfile
		echo -e "\nWARNING: Error in Nanovar calling."
	else
		echo "`date +%F,%R:%S`: Nanovar calling done." >>  $logfile
		echo "Nanovar calling done."
		grep "Seq_depth" $FOLDER/NANOVAR_$2/*NanoVar.log >> $logfile
		firefox $FOLDER"/NANOVAR_$2/nanovar_results/ALL.output.report.html" &
	fi
fi

echo -e "\n######################## NANOVAR VCF CONVERSION (ANNOTSV) #################\n" >> $logfile
echo -e "\n######################## NANOVAR VCF CONVERSION (ANNOTSV) #################\n"

#Creation of compliant vcf files for AnnotSV using the filtered nanovar vcf file ( confidence threshold 0.5)
	#Check wether the nanovar calling ended up successfully by looking at the existence or not of the nanovar vcf file.  

if [ ! -e $FOLDER"/NANOVAR_$2/nanovar_results/ALL.output.filtered-0.5.vcf" ];then
	echo "WARNING: No nanovar vcf available for conversion." >> $logfile
	echo "WARNING: No nanovar vcf available for conversion."
elif [ -e $FOLDER"/NANOVAR_$2/nanovar2annotSV/${1/\//}_nanovar_filtered-0.5_ANSV.vcf" ]; then
	echo "Nanovar converted vcf files (ANNOTSV) with same names already exist. Please check filenames." >> $logfile
	echo "Nanovar converted vcf files (ANNOTSV) with same names already exist. Please check filenames."
else
	echo "Lauch VCF conversion..."
	echo -e "`date +%F,%R:%S`: Launch VCF conversion...">> $logfile
	awk -F "\t|;" 'BEGIN{ OFS="\t"};{if($0 ~ /^#/) {print $0}  
					else { gsub ( "SVRATIO" , "AF" ) ;if ( $5 == "BND" || $5 == "INV" ) {
						if ($9 != "CHR2=.") { gsub ("CHR2=chr", "chr"); gsub ("END=","") ; print $1,$2,$3,"A","<"$5">",$6,$7,$8";CHR2="$9";POS2="$10";"$11";"$12";"$13";"$14";"$15";"$16;
                	                                print $9,$10,$3".end","A","<"$5">",$6,$7,$8";CHR2="$1";POS2="$2";"$11";"$12";"$13";"$14";"$15";"$16} 
						else { gsub ("END=","") ; print $1,$2,$3,"A","<"$5">",$6,$7,$8";"$9";POS2="$10";"$11";"$12";"$13";"$14";"$15";"$16}
                                	                }
                                        	else { print $1,$2,$3,"A","<"$5">",$6,$7,$8";"$9";"$10";"$11";"$12";"$13";"$14";"$15";"$16} 
                                        	}
                                	}' $FOLDER"/NANOVAR_$2/nanovar_results/ALL.output.filtered-0.5.vcf" > $FOLDER"/NANOVAR_$2/nanovar2annotSV/${1/\//}_nanovar_filtered-0.5_ANSV.vcf" 2>> $logfile
	awk -F "\t|;" 'BEGIN{ OFS="\t"};{if($0 ~ /^#/) {print $0}  
                                        else { gsub ( "SVRATIO" , "AF" ) ;if ( $5 == "BND" || $5 == "INV" ) {
                                                if ($9 != "CHR2=.") { gsub ("CHR2=chr", "chr"); gsub ("END=","") ; print $1,$2,$3,"A","<"$5">",$6,$7,$8";CHR2="$9";POS2="$10";"$11";"$12";"$13";"$14";"$15";"$16;
                                                        print $9,$10,$3".end","A","<"$5">",$6,$7,$8";CHR2="$1";POS2="$2";"$11";"$12";"$13";"$14";"$15";"$16} 
                                                else { gsub ("END=","") ; print $1,$2,$3,"A","<"$5">",$6,$7,$8";"$9";POS2="$10";"$11";"$12";"$13";"$14";"$15";"$16}
                                                        }
                                                else { print $1,$2,$3,"A","<"$5">",$6,$7,$8";"$9";"$10";"$11";"$12";"$13";"$14";"$15";"$16} 
                                                }
                                        }' $FOLDER"/NANOVAR_$2/nanovar_results/ALL.output.total.vcf" > $FOLDER"/NANOVAR_$2/nanovar2annotSV/${1/\//}_nanovar_total_ANSV.vcf" 2>> $logfile
	if [ $? -ne 0 ];then
		echo "\nWARNING: Error in nanovar vcf conversion." >> $logfile
		 echo "\nWARNING: Error in nanovar vcf conversion."
	else
		echo "VCF conversion done."
		echo -e "\n`date +%F,%R:%S`: VCF conversion done." >> $logfile
		echo "Now run annotSV to annotate converted vcf right here: https://lbgi.fr/AnnotSV/runjob"
		echo -e "\n`date +%F,%R:%S`: Now run annotSV to annotate converted vcf right here: https://lbgi.fr/AnnotSV/runjob" >> $logfile
		firefox https://lbgi.fr/AnnotSV/runjob &
	fi
fi

echo -e "\n######################## HSBLAST FILE CONVERSION FOR IGV #################\n"
echo -e "\n######################## HSBLAST FILE CONVERSION FOR IGV #################\n" >> $logfile

if [ ! -e $FOLDER/NANOVAR_$2/nanovar_run/hsblast_longreads/ALL.hsblast-$2.tsv ];then
	echo "WARNING: No nanovar hsblast file avaiblable for conversion." >> $logfile
	echo "WARNING: No nanovar hsblast file avaiblable for conversion."
elif [ -e $FOLDER/NANOVAR_$2/hsblast4igv/ALL.hsblast-$2_igv_sorted.bed ]; then
	echo "Hsblast file conversion for visualisation in igv has already been done." >> $logfile
	echo "Hsblast file conversion for visualisation in igv has already been done."
else
	echo "Launch Hsblast file conversion for igv..."
	echo -e "`date +%F,%R:%S`: Launch Hsblast file conversion for igv..."
	awk 'BEGIN{OFS="\t"}{print $1,$2,$2+$3,$5,$12,$8,$2,$2+$3}' $FOLDER"/NANOVAR_$2/nanovar_run/hsblast_longreads/ALL.hsblast-$2.tsv" |  awk -F '\t' -v OFS='\t' '{y=$1"\t"$2"\t"$4; a[y]=$0}END{for (y in a) print a[y]}' >  $FOLDER"/NANOVAR_$2/hsblast4igv/ALL.hsblast-$2.bed" 2>> $logfile
	awk -F '\t' -v OFS='\t' '{x=$4;
        	                y=$1"\t"$2"\t"$3"\t"$4;
                	        d[y]=$1"\t"$2"\t"$3;
                        	b[y]=$4;c[y]=$5"\t"$6"\t"$7"\t"$8;
                        	if(a[x] !~ $1":"){ g[x]=g[x]+1};
                        	a[x]=a[x]"~"$1":"$2"-"$3;
	                        f[x]=f[x]+1}
        	                END{for (y in b ){ if (f[b[y]] == 1 ) {print d[y],b[y]""a[b[y]],c[y],"100,100,100"}
                	                         else { if (f[b[y]] == 2 ) { if (g[b[y]] == 1) { print d[y],b[y]""a[b[y]],c[y],"0,0,255"}
                        	                                                else { print d[y],b[y]""a[b[y]],c[y],"255,0,0" } }
                                	                else {if (f[b[y]] >= 3 ) { if (g[b[y]] == 1) { print d[y],b[y]""a[b[y]],c[y],"0,255,0"}
                                        	                                 else {if(g[b[y]] < f[b[y]] && g[b[y]]==2) {print d[y],b[y]""a[b[y]],c[y],"255,0,0" }
                                                	                                else{print d[y],b[y]""a[b[y]],c[y],"255,127,0"} } 
                                                	}}}}}' $FOLDER"/NANOVAR_$2/hsblast4igv/ALL.hsblast-$2.bed" > $FOLDER"/NANOVAR_$2/hsblast4igv/ALL.hsblast-$2_igv.bed" 2>>$logfile
	bedtools sort -i $FOLDER"/NANOVAR_$2/hsblast4igv/ALL.hsblast-$2_igv.bed" > $FOLDER"/NANOVAR_$2/hsblast4igv/ALL.hsblast-$2_igv_sorted.bed" 2>> $logfile
	if [ $? -ne 0 ];then
                echo "\nWARNING: Error in hsblast file conversion." >> $logfile
		echo "\nWARNING: Error in hsblast file conversion."
        else
                echo "Hsblast file conversion done."
                echo -e "\n`date +%F,%R:%S`: hsblast file conversion done." >> $logfile
        fi

fi

# f = nombre de positions de blast differentes
# g = nombre de chromosome différents trouvé par blast

#nbrChr\nbrBlast	1       2       >=3     
#NORMAL			GRIS    /       /
#INTRA			/       BLEU    VERT
#INTER			/       ROUGE   ORANGE

echo -e "\n######################### FIRST ALIGNMENT #########################\n" >> $logfile
echo -e "\n######################### FIRST ALIGNMENT #########################\n"

#first-alignement:
	#Execute the script only when FASTQ folder is not empty (doesn't only contain ALL folder) and ALL.sam file doesn't already exist.
if [ ! -e $FOLDER/FASTQ/ALL/ALL.fastq ]; then
	echo "WARNING: No ALL.fastq available for alignement. Please check basecalling" >> $logfile
	echo "WARNING: No ALL.fastq available for alignement. Please check basecalling"
elif [ -e $FOLDER"/BAM_$2/ALL.sam" ]; then
	echo "Alignment has already been done." >> $logfile
	echo "Alignment has already been done."
else
	echo "Launch aligner..."
	echo -e "`date +%F,%R:%S`: Lauch aligner..." >> $logfile
	echo "guppy_aligner -i $FOLDER/FASTQ/ALL/ -s $FOLDER/BAM_$2/  --align_ref ~/resources/$2/$2.fa -t 11" >> $logfile
	guppy_aligner -i $FOLDER"/FASTQ/ALL/" -s $FOLDER"/BAM_$2/"  --align_ref ~/resources/$2/$2.fa -t 11 2>> $logfile
	if [ $? -ne 0 ];then
                echo "\nWARNING: Error in Alignment." >> $logfile
		echo "\nWARNING: Error in Alignment."
        else
		echo -e "`date +%F,%R:%S`: Alignment done." >> $logfile
		echo "Alignment done"
	fi
fi

echo -e "\n######################## BAM CONVERSION ###########################\n" >> $logfile
echo -e "\n######################## BAM CONVERSION ###########################\n"
#sam-to-bam:
	#Execute the script only when the ALL.bam file doesn't exist
if [ ! -e $FOLDER/BAM_$2/ALL.sam ];then
	echo "WARNING: No ALL.sam file available. Please check alignment.">> $logfile
	echo "WARNING: No ALL.sam file available. Please check alignment."
elif [ ! -e $FOLDER"/BAM_$2/ALL.bam" ]; then
	echo "Launch BAM conversion..."
	echo -e "`date +%F,%R:%S`: Launch BAM conversion..." >> $logfile
	#~/TOOLS/sambamba_v0.6.6 view -h $FOLDER"/BAM_$2/*sam" -t 11 -f bam
     	# we use samtools view to convert, bgzf error with sambamba
	samtools view -bSh -@ 10 $FOLDER"/BAM_$2/ALL.sam" |samtools sort -@ 10  - > $FOLDER"/BAM_$2/ALL.bam"
	samtools index -@ 10 $FOLDER"/BAM_$2/ALL.bam" 2>> $logfile
	if [ $? -ne 0 ];then
                echo "\nWARNING: Error in BAM conversion." >> $logfile
		echo "\nWARNING: Error in BAM conversion."
        else
		echo -e "`date +%F,%R:%S`: Bam conversion done." >> $logfile
		echo "Bam conversion done."
	fi
else
	echo "\nBAM conversion has already been done." >> $logfile
fi

echo -e "\n######################### NANOSV CALLING ##########################\n" >> $logfile
echo -e "\n######################### NANOSV CALLING ##########################\n"
#NanoSV calling: 
	#Execute the script only when ALL.bam file exists and ALL.vcf doesn't

if [ ! -e  $FOLDER/BAM_$2/ALL.bam ]; then
	echo "WARNING: ALL.bam file doesn't exist. NanoSV couldn't be lauchned. Please check alignment." >> $logfile
	echo "WARNING: ALL.bam file doesn't exist. NanoSV couldn't be lauchned. Please check alignment."
elif [ -e $FOLDER/NANOSV_$2/ALL.vcf ];then
	echo "NanoSV calling has already been done." >> $logfile
	echo "NanoSV calling has already been done."
else
	echo "Launch NanoSV calling..."
	echo "`date +%F,%R:%S`: Launch NanoSV calling..." >> $logfile
	echo "python3 ~/TOOLS/nanosv/nanosv/NanoSV.py -t 10 -s ~/TOOLS/sambamba-0.6.8-linux-static -c ~/TOOLS/nanosv/nanosv/config.ini -b ~/resources/$2/nanosv_human_$2.bed -o $FOLDER/NANOSV_$2/ALL.vcf   $FOLDER/BAM_$2/ALL.bam" >> $logfile
	python3 ~/TOOLS/nanosv/nanosv/NanoSV.py -t 10 -s ~/TOOLS/sambamba-0.6.8-linux-static -c ~/TOOLS/nanosv/nanosv/config.ini -b ~/resources/$2/nanosv_human_$2.bed -o $FOLDER"/NANOSV_$2/ALL.vcf"    $FOLDER"/BAM_$2/ALL.bam"  2>> $logfile
	if [ $? -ne 0 ];then
                echo "\nWARNING: Error in NanovaSV calling." >> $logfile
		echo "\nWARNING: Error in NanovaSV calling."
        else
		echo "NanoSV calling done."
		echo -e "`date +%F,%R:%S`: NanoSV calling done." >> $logfile
	fi
fi

echo -e "\nRun ended: `date +%F,%R:%S`" >> $logfile
echo "Run done"
exit 0
