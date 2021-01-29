#!/usr/bin/perl

##############################################################################################
# knotAnnotSV 1.0                                                                            #	
#                                                                                            #
# knotAnnotSV: Creation of a customizable html file to visualize, filter                     # 
#                   and analyze an AnnotSV output                                            #
#                                                                                            #
# Author: Thomas Guignard 2020-2021                                                          #
#                                                                                            #
# Copyright (C) 2020-2021 Thomas Guignard (t-guignard@chu-montpellier.fr)                    #
#                                                                                            #
# This is part of knotAnnotSV source code.                                                   #
#                                                                                            #
# This program is free software; you can redistribute it and/or                              #
# modify it under the terms of the GNU General Public License                                #
# as published by the Free Software Foundation; either version 3                             #
# of the License, or (at your option) any later version.                                     #
#                                                                                            #
# This program is distributed in the hope that it will be useful,                            #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                             #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                              #
# GNU General Public License for more details.                                               #
#                                                                                            #
# You should have received a copy of the GNU General Public License                          #
# along with this program; If not, see <http://www.gnu.org/licenses/>.                       #
##############################################################################################


use strict; 
use warnings;
use Getopt::Long; 
use YAML::XS 'LoadFile';
use Sort::Key::Natural qw(rnatkeysort);
use File::Basename;

#use Switch;
#use Data::Dumper;

#parameters
my $man = "USAGE : \nperl knotAnnotSV.pl 
\n--configFile <YAML config file for customizing output>
\n--annotSVfile <AnnotSV annotated file> 
\n--outDir <output directory (default = current dir)> 
\n--outPrefix <output file prefix (default = \"\")> 
\n--genomeBuild <Genome Assembly reference (default = hg19)> 
\n--LOEUFcolorRange <Number to define which color to use for LOEUF bin: 1 (red-to-green), 2 (red-shades-only). (default = 1)> 
\n--datatableDir <Local Path to DataTables directory containing css and js files (default = \"\", requires web connection)>"; 


my $configFile = ".";

my $help;
my $current_line;
my $incfile = "";
my $outDir = ".";
my $outPrefix = "";
my $annotSVranking = "";
my $datatableDir = "";
my $LOEUFcolorRange = "";

#vcf parsing and output
my @line;
my $variantID; 
my $count = 0;

my @finalSortData;
my %hashFinalSortData;
	
my $config;
my %columnHash;
my %InColHash;
my %OutColHash;
my %NameColHash;
my %dataHash;
my %dataCommentHash;
my %dataColorHash;

my $scorePenalty;
my $fullSplitScore;
my %SV_ID;
my $SV_type;
my $fullRowColor;

my $url2UCSC="";
my $url2OMIM="";
my $url2DECIPHER="";
my $url2ensembl="";
						
my @criteria;	
my @OMIM_ID_array;
my @OMIM_phen_array;
my @GeneName_array;
my @RE_gene_array;

my %debugHash;

my $genomeBuild="";
#style alignment for tooltiptext
my $alignTooltiptext="";

GetOptions( "annotSVfile=s"		=> \$incfile,
			"configFile=s"		=> \$config,
			"outDir=s"			=> \$outDir,
			"outPrefix:s"		=> \$outPrefix,
			"datatableDir=s"	=> \$datatableDir,
			"genomeBuild=s"		=> \$genomeBuild,
			"LOEUFcolorRange=s"	=> \$LOEUFcolorRange,
			"help|h"			=> \$help);
				
				


#check mandatory arguments
if(defined $help || $incfile eq ""){
	die($man."\n");
}

#add underscore to output prefix
if($outPrefix ne ""){
	$outPrefix .= "_";
}

#check datatable pathname
if ($datatableDir ne "" && ! -d $datatableDir){
    print "\nError with datatableDir argument: ".$datatableDir." is not an existing directory.\n\n";
    exit 1;    
}
			
#check if genomeBuild exists (hg19 grch37 hg38, mm etc..) and put hg19 as default build 
if($genomeBuild eq ""){
	$genomeBuild = "hg19";
}


#open( CONFIG, "<$config" ) or die ("Cannot open vcf file $config") ;

#CONFIG FILE PARSING
print STDERR "Parsing config file....\n";
my %configHash;
my $configHash;
my $OutColCounter = 0;
my $outPOSITION = 0; 
$configHash = LoadFile($config);
#print Dumper($configHash);


foreach my $item (keys %{$configHash}){
		#print $item."\n";
		my $field = $item;
		#print $configHash->{$field}->{POSITION}."\n";

        if ($configHash->{$field}->{POSITION} != 0 && defined $OutColHash{$configHash->{$field}->{POSITION}}){
            print "\nError in config file: 'POSITION : ".$configHash->{$field}->{POSITION}."' is assigned twice.\nPlease correct config file to avoid 2 identical positions.\n\n";
            exit 1;
        }else {
		    if ($configHash->{$field}->{POSITION} != 0){
				#$OutColHash{$configHash->{$field}->{POSITION}} = $field;
				$outPOSITION = $configHash->{$field}->{POSITION};
				$OutColHash{$outPOSITION}{'field'} = $field;
				#$OutColHash{$configHash->{$field}->{POSITION}}{'field'} = $field;
				$OutColCounter ++;

            }
		    $NameColHash{$field} = $configHash->{$field}->{POSITION};
			#TODO $NameColHash{$field}{'POSITION'} = $configHash->{$field}->{POSITION};
			
			if (defined $configHash->{$field}->{COMMENTLIST}){
				$dataCommentHash{$field}{'commentFieldList'}=$configHash->{$field}->{COMMENTLIST} ;
			}
			#get custom name for column
			if (defined $configHash->{$field}->{RENAME}){
			#$NameColHash{$field}{'RENAME'} = $configHash->{$field}->{RENAME};
				$OutColHash{$configHash->{$field}->{POSITION}}{'RENAME'} = $configHash->{$field}->{RENAME};
			}else{
			#$NameColHash{$field}{'RENAME'} = $field;
				$OutColHash{$configHash->{$field}->{POSITION}}{'RENAME'} = $field;
			}
			#get custom string for column header tooltip
			if (defined $configHash->{$field}->{HEADERTIPS}){
				#$NameColHash{$field}{'HEADERTIPS'} = $configHash->{$field}->{HEADERTIPS};
				$OutColHash{$configHash->{$field}->{POSITION}}{'HEADERTIPS'} = $configHash->{$field}->{HEADERTIPS};
			}
        }
		
		#print Dump($field)."\n";
}


#ANNOTSV RANKING PARSING to get ranking decision
my @SVrankLine;
my %SVrankHash;

if ($annotSVranking ne ""){
	open( SVRANK , "<$annotSVranking" )or die("Cannot open ranking file ".$annotSVranking."\n") ;
	
	
	while( <SVRANK> ){

		chomp $_;
		@SVrankLine = split( /\t/, $_ );	

		if (defined $SVrankHash{$SVrankLine[0]."_".$SVrankLine[2]}){
			print "Doublon detected in ranking file : ".$_."\n\n";
		}else{
			$SVrankHash{$SVrankLine[0]."_".$SVrankLine[2]} = $SVrankLine[4];
		}	
	}

	#add comment with rank to AnnotSV ranking column
	if (! defined $dataCommentHash{'ACMG_class'}{'commentFieldList'}){
		$dataCommentHash{'ACMG_class'}{'SVrank'} = "OK";
	}

}





#Hash of ACMG incidentalome genes  => grey color #808080
my %ACMGgene = ("ACTA2" =>1,"ACTC1" =>1,"APC" =>1,"APOB" =>1,"ATP7B" =>1,"BMPR1A" =>1,"BRCA1" =>1,"BRCA2" =>1,"CACNA1S" =>1,"COL3A1" =>1,"DSC2" =>1,"DSG2" =>1,"DSP" =>1,"FBN1" =>1,"GLA" =>1,"KCNH2" =>1,"KCNQ1" =>1,"LDLR" =>1,"LMNA" =>1,"MEN1" =>1,"MLH1" =>1,"MSH2" =>1,"MSH6" =>1,"MUTYH" =>1,"MYBPC3" =>1,"MYH11" =>1,"MYH7" =>1,"MYL2" =>1,"MYL3" =>1,"NF2" =>1,"OTC" =>1,"PCSK9" =>1,"PKP2" =>1,"PMS2" =>1,"PRKAG2" =>1,"PTEN" =>1,"RB1" =>1,"RET" =>1,"RYR1" =>1,"RYR2" =>1,"SCN5A" =>1,"SDHAF2" =>1,"SDHB" =>1,"SDHC" =>1,"SDHD" =>1,"SMAD3" =>1,"SMAD4" =>1,"STK11" =>1,"TGFBR1" =>1,"TGFBR2" =>1,"TMEM43" =>1,"TNNI3" =>1,"TNNT2" =>1,"TP53" =>1,"TPM1" =>1,"TSC1" =>1,"TSC2" =>1,"VHL" =>1,"WT1"=>1);


#initialize gene colore according to pLI/LOEUF
#my %pLI_ColorHash;
#my $pLI_Color;

#$pLI_ColorHash{'0.9'} = '#FF0000';
#$pLI_ColorHash{'0.8'} = '#FF3300';
#$pLI_ColorHash{'0.7'} = '#FF6600';
#$pLI_ColorHash{'0.6'} = '#FF9900';
#$pLI_ColorHash{'0.5'} = '#FFCC00';
#$pLI_ColorHash{'0.4'} = '#FFFF00';
#$pLI_ColorHash{'0.3'} = '#BFFF00';
#$pLI_ColorHash{'0.2'} = '#7FFF00';
#$pLI_ColorHash{'0.1'} = '#3FFF00';
#$pLI_ColorHash{'0.0'} = '#00FF00';


#Select and initialize color Range for LOEUF bin
my %LOEUF_ColorHash;
my $LOEUF_Color;
if($LOEUFcolorRange eq "" ||$LOEUFcolorRange == 1 ){
	$LOEUF_ColorHash{'0.0'} = '#FF0000';
	$LOEUF_ColorHash{'1.0'} = '#FF3300';
	$LOEUF_ColorHash{'2.0'} = '#FF6600';
	$LOEUF_ColorHash{'3.0'} = '#FF9900';
	$LOEUF_ColorHash{'4.0'} = '#FFCC00';
	$LOEUF_ColorHash{'5.0'} = '#FFFF00';
	$LOEUF_ColorHash{'6.0'} = '#BFFF00';
	$LOEUF_ColorHash{'7.0'} = '#7FFF00';
	$LOEUF_ColorHash{'8.0'} = '#3FFF00';
	$LOEUF_ColorHash{'9.0'} = '#00FF00';

}elsif($LOEUFcolorRange == 2){
	$LOEUF_ColorHash{'0.0'} = '#FF0606';
	$LOEUF_ColorHash{'1.0'} = '#FF3737';
	$LOEUF_ColorHash{'2.0'} = '#FF5050';
	$LOEUF_ColorHash{'3.0'} = '#FF6363';
	$LOEUF_ColorHash{'4.0'} = '#FF6F6F';
	$LOEUF_ColorHash{'5.0'} = '#FF8888';
	$LOEUF_ColorHash{'6.0'} = '#FF9B9B';
	$LOEUF_ColorHash{'7.0'} = '#FFB4B4';
	$LOEUF_ColorHash{'8.0'} = '#FFD9D9';
	$LOEUF_ColorHash{'9.0'} = '#FFF2F2';

}else{
	die("LOEUFcolorRange argument is unacceptable: ".$LOEUFcolorRange."\n") ;
}


### Initialize column Name Mode Hash

my %colNameMode = (
"AnnotSV_ID"=>"fullsplit",
"SV_chrom"=>"fullsplit",
"SV_start"=>"fullsplit",
"SV_end"=>"fullsplit",
"SV_length"=>"fullsplit",
"SV_type"=>"fullsplit",
"Samples_ID"=>"fullsplit",
"REF"=>"fullsplit",
"ALT"=>"fullsplit",
"FORMAT"=>"fullsplit",
"Annotation_mode"=>"fullsplit",
"Gene_name"=>"fullsplit",
"Gene_count"=>"full",
"Tx"=>"split",
"Tx_start"=>"split",
"Tx_end"=>"split",
"Overlapped_tx_length"=>"split",
"Overlapped_CDS_length"=>"split",
"Overlapped_CDS_percent"=>"split",
"Frameshift"=>"split",
"Exon_count"=>"split",
"Location"=>"split",
"Location2"=>"split",
"Dist_nearest_SS"=>"split",
"Nearest_SS_type"=>"split",
"Intersect_start"=>"split",
"Intersect_end"=>"split",
"RE_gene"=>"full",
"B_gain_source"=>"fullsplit",
"B_gain_coord"=>"fullsplit",
"B_loss_source"=>"fullsplit",
"B_loss_coord"=>"fullsplit",
"B_ins_source"=>"fullsplit",
"B_ins_coord"=>"fullsplit",
"B_inv_source"=>"fullsplit",
"B_inv_coord"=>"fullsplit",
"P_gain_phen"=>"fullsplit",
"P_gain_hpo"=>"fullsplit",
"P_gain_source"=>"fullsplit",
"P_gain_coord"=>"fullsplit",
"P_loss_phen"=>"fullsplit",
"P_loss_hpo"=>"fullsplit",
"P_loss_source"=>"fullsplit",
"P_loss_coord"=>"fullsplit",
"P_ins_phen"=>"fullsplit",
"P_ins_hpo"=>"fullsplit",
"P_ins_source"=>"fullsplit",
"P_ins_coord"=>"fullsplit",
"P_snvindel_nb"=>"fullsplit",
"P_snvindel_phen"=>"fullsplit",
"Cosmic_ID"=>"fullsplit",
"Cosmic_mut_typ"=>"fullsplit",
"TAD_coordinate"=>"full",
"ENCODE_experiment"=>"full",
"GC_content_left"=>"full",
"GC_content_right"=>"full",
"Repeat_coord_left"=>"full",
"Repeat_type_left"=>"full",
"Repeat_coord_right"=>"full",
"Repeat_type_right"=>"full",
"Gap_left"=>"full",
"Gap_right"=>"full",
"ENCODE_blacklist_left"=>"full",
"ENCODE_blacklist_characteristics_left"=>"full",
"ENCODE_blacklist_right"=>"full",
"ENCODE_blacklist_characteristics_right"=>"full",
"SegDup_left"=>"full",
"SegDup_right"=>"full",
"ACMG"=>"split",
"HI"=>"fullsplit",
"TS"=>"fullsplit",
"DDD_mode"=>"split",
"DDD_consequence"=>"split",
"DDD_disease"=>"split",
"DDD_pmid"=>"split",
"DDD_HI_percent"=>"fullsplit",
"ExAC_synZ"=>"fullsplit",
"ExAC_misZ"=>"fullsplit",
"ExAC_delZ"=>"fullsplit",
"ExAC_dupZ"=>"fullsplit",
"ExAC_cnvZ"=>"fullsplit",
"LOEUF_bin"=>"fullsplit",
"GnomAD_pLI"=>"fullsplit",
"ExAC_pLI"=>"fullsplit",
"OMIM_morbid"=>"fullsplit",
"OMIM_morbid_candidate"=>"fullsplit",
"OMIM_ID"=>"fullsplit",
"OMIM_phenotype"=>"split",
"OMIM_inheritance"=>"split",
"Exomiser_gene_pheno_score"=>"fullsplit",
"Human_pheno_evidence"=>"split",
"Mouse_pheno_evidence"=>"split",
"Fish_pheno_evidence"=>"split",
"Compound_htz(sample)"=>"fullsplit",
"Count_hom(sample)"=>"fullsplit",
"Count_htz(sample)"=>"fullsplit",
"Count_htz/allHom(sample)"=>"fullsplit",
"Count_htz/total(cohort) "=>"fullsplit",
"Count_total(cohort)"=>"fullsplit",
"AnnotSV_ranking_score"=>"full",
"AnnotSV_ranking_criteria"=>"full",
"ACMG_class"=>"full");

# Initialisation of ACMG criteria
my %gainRankCriteria = (
"1A"=>"Contains protein-coding or other known functionally important elements. ",
"1B"=>"Does NOT contain protein-coding or any known functionally important elements.",
"2A"=>"Complete overlap; the known pathogenic gain SV is fully contained within the observed copy-number gain.",
"2B"=>"Partial overlap of a known pathogenic gain SV. The observed CNV does NOT contain the TS gene or critical region for this known pathogenic Gain SV OR Unclear if the known causative gene or critical region is affected OR No specific causative gene or critical region has been established for this known pathogenic SV.",
"2C"=>"Identical in gene content to the established benign copy-number gain.",
"2D"=>"Smaller than established benign copy-number gain, breakpoint(s) does not interrupt protein-coding genes.",
"2E"=>"Smaller than established benign copy-number gain, breakpoint(s) potentially interrupts protein-coding gene.",
"2F"=>"Larger than known benign copy-number gain, does not include additional protein-coding genes.",
"2G"=>"Overlaps a benign copy-number gain but includes additional genomic material.",
"2H-1"=>"HI gene / morbid gene fully contained within observed copy-number gain and patient's phenotype is highly specific and consistent with what is expected for LOF of that gene (Exomiser_gene_pheno_score >= 0.7).",
"2H-2"=>"HI gene / morbid gene fully contained within observed copy-number gain and patient's phenotype is nonspecific with what is expected for LOF of that gene (Exomiser_gene_pheno_score < 0.7).",
"2I"=>"Both breakpoints are within the same HI gene / morbid gene (gene-level sequence variant, possibly resulting in loss of function [LOF]) and disrupts the reading frame.",
"2J"=>"One breakpoint is within an established HI gene / morbid gene, patient’s phenotype is either inconsistent with what is expected for LOF of that gene OR unknown.",
"2K"=>"One breakpoint is within an established HI gene / morbid gene, patient’s phenotype is highly specific and consistent with what is expected for LOF of that gene (Exomiser_gene_pheno_score > 0.7).",
"2L"=>"One or both breakpoints are within gene(s) of no established clinical significance.",
"3A"=>"0–34 genes wholly or partially included",
"3B"=>"35–49 genes wholly or partially included",
"3C"=>"50 or more genes wholly or partially included",
"5F"=>"Inheritance information is unavailable or uninformative.",
"5G"=>"The patient phenotype is nonspecific, but is consistent with what has been described in similar cases (EXOMISER_GENE_PHENO_SCORE > 0.5).",
"5H"=>"The patient phenotype is highly specific and consistent with what has been described in similar cases (EXOMISER_GENE_PHENO_SCORE > 0.7).");


my %lossRankCriteria = (
"1A"=>"Contains protein-coding or other known functionally important elements. ",
"1B"=>"Does NOT contain protein-coding or any known functionally important elements.",
"2A"=>"Complete overlap of a known pathogenic Loss SV.",
"2B"=>"Partial overlap of a known pathogenic Loss SV. The observed CNV does NOT contain a HI gene OR Unclear if known causative gene or critical region is affected OR No specific causative gene or critical region has been established for this known pathogenic Loss SV",
"2C-1"=>"Partial overlap with the 5' end of an established HI gene / morbid gene (3' end of the gene not involved) and coding sequence is involved",
"2C-2"=>"Partial overlap with the 5' end of an established HI gene / morbid gene (3' end of the gene not involved) and only the 5' UTR is involved",
"2D"=>"Partial overlap with the 3' end of an established HI gene / morbid gene (5' end of the gene not involved)...",
"2D-1"=>"Partial overlap with the 3' end of an established HI gene / morbid gene (5' end of the gene not involved) and only the 3' untranslated region is involved.",
"2D-2"=>"Partial overlap with the 3' end of an established HI gene / morbid gene (5' end of the gene not involved) and only the last exon is involved. Other established pathogenic snv/indel have been reported in the observed CNV",
"2D-3"=>"Partial overlap with the 3' end of an established HI gene / morbid gene (5' end of the gene not involved) and only the last exon is involved. No other established pathogenic variants have been reported in the observed CNV.",
"2D-4"=>"Partial overlap with the 3' end of an established HI gene / morbid gene (5' end of the gene not involved) and it includes other exons in addition to the last exon. Nonsense-mediated decay is expected to occur.",
"2E"=>"Both breakpoints are within the same HI gene / morbid gene (intragenic CNV; gene-level sequence variant)...",
"2E-1"=>"Both breakpoints are within the same HI gene / morbid gene (intragenic CNV; gene-level sequence variant) and disrupts the reading frame",
"2E-2"=>"Both breakpoints are within the same HI gene / morbid gene (intragenic CNV; gene-level sequence variant) and >=1 exon deleted AND other established pathogenic snv/indel have been reported in the observed CNV AND variant removes >= 10% of protein",
"2E-3"=>"Both breakpoints are within the same HI gene / morbid gene (intragenic CNV; gene-level sequence variant) and >=1 exon deleted AND other established pathogenic snv/indel have been reported in the observed CNV AND variant removes < 10% of protein",
"2E-4"=>"Both breakpoints are within the same HI gene / morbid gene (intragenic CNV; gene-level sequence variant) and >=1 exon deleted AND no established pathogenic snv/indel have been reported in the observed CNV AND variant removes > 10% of protein",
"2F"=>"Completely contained within an established benign CNV region.",
"2G"=>"Overlaps an established benign CNV, but includes additional genomic material.",
"2H"=>"Two or more HI predictors suggest that AT LEAST ONE gene in the interval is HI (gnomAD pLI >=0.9 and DECIPHER HI index <=10%)",
"3A"=>"0–24 genes wholly or partially included",
"3B"=>"25–34 genes wholly or partially included",
"3C"=>"35+ genes wholly or partially included",
"5F"=>"Inheritance information is unavailable or uninformative.",
"5G"=>"The patient phenotype is nonspecific, but is consistent with what has been described in similar cases (EXOMISER_GENE_PHENO_SCORE > 0.5).",
"5H"=>"The patient phenotype is highly specific and consistent with what has been described in similar cases (EXOMISER_GENE_PHENO_SCORE > 0.7).");



####################################################################
#############################################
##################   Start parsing VCF

open( VCF , "<$incfile" )or die("Cannot open annotSV file ".$incfile."\n") ;

my ($outBasename,$dir,$ext) = fileparse($incfile,'\.tsv$');


while( <VCF> ){
  	$current_line = $_;
	$count++;


	chomp $current_line;
	@line = split( /\t/, $current_line );	
	

#############################################
##############   Treatment for First line to create header of the output

	if ( $line[0] eq "AnnotSV_ID" )   {
		
		#TODO check if header contains required INFO
		print STDERR "Parsing AnnotSV header in order to get column names ... \n";

		#initialize column order
		for( my $fieldNbr = 0 ; $fieldNbr < scalar @line; $fieldNbr++){
			
			#print STDERR "Field found:  ".$line[$fieldNbr]."\n";
			$InColHash{$fieldNbr} = $line[$fieldNbr];
			
			$dataHash{$line[$fieldNbr]} = ".";
			
		}
		
		foreach my $field (keys %NameColHash){
			
			if (! defined $dataHash{$field}){
				$dataHash{$field} = "NA";
			}
		}

		next;
		
#############################################
##############################
##########  start to compute variant lines	

	}else {

		#initialise final printable string
		@finalSortData = ("");
		#$pLI_Color=".";
		$LOEUF_Color=".";

        #fill nbr of SV_ID;
        if (defined $SV_ID{$line[0]} ){
            $SV_ID{$line[0]}{'nbr'}++;
        }else{
            $SV_ID{$line[0]}{'nbr'} = 1;
        }
        
		#reinitialize Comment Values
		foreach my $field (keys %dataCommentHash){
			delete $dataCommentHash{$field}{'values'};
		}
		
		#reinitialize dataHash Values
		foreach my $field (keys %dataHash){
			unless($dataHash{$field} eq "NA" ){$dataHash{$field} = ".";}
		}


		#fill printable string
		for( my $fieldNbr = 0 ; $fieldNbr < scalar @line; $fieldNbr++){
			
			if($line[$fieldNbr] ne ""){
				$line[$fieldNbr] =~ s/;/; /g ;  
				#$line[$fieldNbr] =~ s/\//\/ /g ;  	
				$line[$fieldNbr] =~ s/</&lt;/g ;  
				$line[$fieldNbr] =~ s/>/&gt;/g ;  

				$dataHash{$InColHash{$fieldNbr}} = $line[$fieldNbr] ; 
				
			}else{
				$dataHash{$InColHash{$fieldNbr}} = ".";
			}
		}

		#remove "-" sign for SV length
		if (defined $dataHash{'SV_length'}){
			$dataHash{'SV_length'} =~ s/^\-//;
		}



		##### attribute URL to genes and OMIM

		my $GeneName_link_string="";
		my $OMIM_ID_link_string="";
		my $OMIM_phen_link_string="";
		@OMIM_ID_array = "";
		@OMIM_phen_array = "";
		@GeneName_array = "";
		@RE_gene_array = "";
		my %GeneName_hash ;
		my %RE_gene_hash ;
		my $tempString;
		#add url to OMIM_ID
			if ($dataHash{"OMIM_ID"} ne "."){
				$tempString = "";
				@OMIM_ID_array = split(/; /, $dataHash{"OMIM_ID"} );
				for( my $ID = 0 ; $ID < scalar @OMIM_ID_array; $ID++){
					$tempString .=    "<a href=\"https://www.omim.org/entry/".$OMIM_ID_array[$ID]."\" target=\"_blank\" rel=\"noopener noreferrer\" >".$OMIM_ID_array[$ID]."; </a>";
					if($ID < 5){
						$OMIM_ID_link_string .= "<a href=\"https://www.omim.org/entry/".$OMIM_ID_array[$ID]."\" target=\"_blank\" rel=\"noopener noreferrer\" >".$OMIM_ID_array[$ID]."; </a>";
					}
				}				
				$dataHash{"OMIM_ID"} = $tempString;
			}else{
				$OMIM_ID_link_string = ".";
			}

			if ($dataHash{"OMIM_phenotype"} ne "."){
				$tempString = "";
				@OMIM_phen_array = split(/; /, $dataHash{"OMIM_phenotype"} );
				if(@OMIM_phen_array){
					for ( my $ID = 0 ; $ID < scalar @OMIM_phen_array ; $ID++){
						if( $OMIM_phen_array[$ID] =~ m/^(.+?)(\d{6})(.+?)$/){
							$tempString .=   $1."<a href=\"https://www.omim.org/entry/".$2."\" target=\"_blank\" rel=\"noopener noreferrer\" >".$2."</a>".$3.";<br>";   	
						
							#DEBUG TODO not sure
							if ( $ID == 0 ){
								if ( scalar @OMIM_phen_array > 1){
									$OMIM_phen_link_string = $1."<a href=\"https://www.omim.org/entry/".$2."\" target=\"_blank\" rel=\"noopener noreferrer\" >".$2."</a>".$3."; [...".scalar @OMIM_phen_array." pheno]";
								}else{
									$OMIM_phen_link_string = $1."<a href=\"https://www.omim.org/entry/".$2."\" target=\"_blank\" rel=\"noopener noreferrer\" >".$2."</a>".$3;
								}
							}
						}else{
							$tempString .= $OMIM_phen_array[$ID].";<br>";
						}
					}
					$tempString =~ s/<br>$//;
					$dataHash{"OMIM_phenotype"} = $tempString ;
					if ($OMIM_phen_link_string eq ""){
						$OMIM_phen_link_string = $dataHash{"OMIM_phenotype"};
					}

				}else{
					$dataHash{"OMIM_phenotype"} =~ s/; /;<br>/g ;
					$OMIM_phen_link_string = ".";
				}
			}else{
				if ($dataHash{"OMIM_morbid"} eq "yes"){
					$OMIM_phen_link_string = "yes";
				}else{
					$OMIM_phen_link_string = ".";
				}

			}
	

			#URL for gene name
			if ($dataHash{"Gene_name"} ne "."){
				$tempString = "";
				@GeneName_array = split(/; /, $dataHash{"Gene_name"} );
				for( my $ID = 0 ; $ID < scalar @GeneName_array; $ID++){
					$GeneName_hash{$GeneName_array[$ID]} = 1;
					$tempString .=    "<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=".$GeneName_array[$ID]."\" target=\"_blank\" rel=\"noopener noreferrer\">".$GeneName_array[$ID]."; </a>";
					if($ID < 5){
						$GeneName_link_string .= "<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=".$GeneName_array[$ID]."\" target=\"_blank\" rel=\"noopener noreferrer\" >".$GeneName_array[$ID]."</a>; ";
					}
				}			

				$dataHash{"Gene_name"} = $tempString  ;
				if (scalar @GeneName_array > 1){
					$GeneName_link_string .= " [...".$dataHash{"Gene_count"}."genes]";
				}else{
					$GeneName_link_string =~ s/; $//;
				}
			}else{
				$GeneName_link_string = ".";
			}
			

			#Classify RE_Gene elements according to Gene name list to highlight missing elements
			if ($dataHash{"RE_gene"} ne "."){
				$tempString = "";
				my $splicebool = 0;
				$dataHash{"RE_gene"} =~ s/; morbid/\/morbid/g;
				@RE_gene_array = split(/; /, $dataHash{"RE_gene"} );
				for ( my $ID = 0 ; $ID < scalar @RE_gene_array ; $ID++){
					if( $RE_gene_array[$ID] =~ m/^(\S+)\s?/ ){
						#print $RE_gene_array[$ID]."__ge__".$1."\n";
						if (defined $GeneName_hash{$1}){
							$RE_gene_hash{$RE_gene_array[$ID]} = 0;
						}else{
							$RE_gene_hash{$RE_gene_array[$ID]} = 1000;
							if ($RE_gene_array[$ID] =~ /morbid/){
								$RE_gene_hash{$RE_gene_array[$ID]} += 10;
							}
							if ($RE_gene_array[$ID] =~ m/EX=(\d\.\d{4})/){
								$RE_gene_hash{$RE_gene_array[$ID]} += 20 * $1;
								#print $RE_gene_array[$ID]."__ex__".$1."\n";
							}

						}
					}
				}	
                foreach my $ReGene (sort {$RE_gene_hash{$b} <=> $RE_gene_hash{$a}} keys %RE_gene_hash){
					if ($splicebool == 0 && $RE_gene_hash{$ReGene} == 0 && $tempString ne ""){
						$splicebool =1;
						$tempString .= "<br>------<br>";
					}
					$tempString .= $ReGene."; ";
				}
				$dataHash{"RE_gene"} = $tempString;
			}
			



		#get SVtype
		$SV_type = $dataHash{'SV_type'};

		#check what type among different writting allowed
		if ($SV_type =~ /DUP|GAIN|CN[3-9]/i){
			$SV_type = "DUP";
			$fullRowColor = "#DFE9FF";
		}elsif ($SV_type =~ /DEL|LOSS|CN[0-1]/i){
			$SV_type = "DEL";
			$fullRowColor = "#F7D8DD";
		}elsif ($SV_type =~ /INV/i){
			$SV_type = "INV";
			$fullRowColor = "#FEE7CD";
		}elsif ($SV_type =~ /INS/i){
			$SV_type = "INS";
			$fullRowColor = "#FFCCEF";
		}elsif ($SV_type =~ /CPX/i){
			$SV_type = "CPX";
			$fullRowColor = "#D4F7DC";
		}elsif ($SV_type =~ /BND/i){
			$SV_type = "BND";
			$fullRowColor = "#FFFF99";
		}else{
			$fullRowColor = "#E0CCFF";
		}

		

		#hash to avoid duplicated comments for P_ and B_ loss and gain
		my %commentDuplicate;
		my $correctFieldCom;

		#get all comments values
		foreach my $field (keys %dataCommentHash){
            #print $field."\n";
			#if (defined $dataCommentHash{$field})
			#	foreach my $fieldCom (@{$dataCommentHash{$field}{'commentFieldList'}})
			if (defined $dataCommentHash{$field}){
                if (! defined $dataCommentHash{$field}{'values'}){
					#special treatment for field , adding to comment the switched name
					if ( $field =~ /^[BP]_/){
						if ($SV_type eq "DEL" ){ 
							$correctFieldCom = $field =~ s/gain|ins|inv/loss/r;
						}elsif ($SV_type eq "DUP" ){ 
							$correctFieldCom = $field =~ s/loss|ins|inv/gain/r;
						}else{
							$correctFieldCom = $field;
						}
						$commentDuplicate{$correctFieldCom} += 1;	
						#add in first line of comment
               		 	$dataCommentHash{$field}{'values'} .= "<span class=\"commentTitle\">".$correctFieldCom . " :</span> ".$dataHash{$correctFieldCom};
					}


				    foreach my $fieldCom (@{$dataCommentHash{$field}{'commentFieldList'}}){
					    
						#skip this field if it doesn't belong to annotation mode
						if (defined $colNameMode{$fieldCom}){ 
							if ($dataHash{"Annotation_mode"} eq "full"){
								if ($colNameMode{$fieldCom} =~ /^split/){
									next;
								}
							}else{
								if ($colNameMode{$fieldCom} =~ /full$/){
									next;
								}
							}
						}


						if (defined $dataHash{$fieldCom}){
							if ( $fieldCom =~ /^[BP]_/){
								if ($SV_type eq "DEL" ){ 
									$correctFieldCom = $fieldCom =~ s/gain|ins|inv/loss/r;
								}elsif ($SV_type eq "DUP" ){ 
									$correctFieldCom = $fieldCom =~ s/loss|ins|inv/gain/r;
								}else{
									$correctFieldCom = $fieldCom;
								}
								$commentDuplicate{$correctFieldCom} += 1;	
								if (defined $commentDuplicate{$correctFieldCom} && $commentDuplicate{$correctFieldCom} > 1){
									next;
								}
								#print $field."_".$correctFieldCom."\n";
                            	$dataCommentHash{$field}{'values'} .= "<br><span class=\"commentTitle\">".$correctFieldCom . " :</span> ".$dataHash{$correctFieldCom};
							}else{
								
								
								if($fieldCom eq "AnnotSV_ranking_criteria" && $dataHash{$fieldCom} ne "." ){
                            		$dataCommentHash{$field}{'values'} .= "<br><span class=\"commentTitle\">".$fieldCom . " :</span><br>";
									
									@criteria = split( /; /, $dataHash{$fieldCom} );
									
									foreach my $crit (@criteria){
										$crit =~ /^([1-5]\w\-?\d?\d?)\s?\(?/;
								
										#print "DEBUG2  :  ".$crit."_____".$1."\n";

										if ($SV_type eq "DEL" ){ 
											if(defined $lossRankCriteria{$1}){
												$dataCommentHash{$field}{'values'} .= $crit.": ".$lossRankCriteria{$1}."<br>";  	
											}else{
												$dataCommentHash{$field}{'values'} .= $crit.": NA<br>";  	
											}
										}elsif ($SV_type eq "DUP" ){ 
											if(defined $gainRankCriteria{$1}){
												$dataCommentHash{$field}{'values'} .= $crit.": ".$gainRankCriteria{$1}."<br>";   	
											}else{
												$dataCommentHash{$field}{'values'} .= $crit.": NA<br>";   	
											}
										}else{
											$dataCommentHash{$field}{'values'} .= $dataHash{$fieldCom};
										}
									}
								

								}else{
                            		$dataCommentHash{$field}{'values'} .= "<br><span class=\"commentTitle\">".$fieldCom . " :</span> ".$dataHash{$fieldCom};
								}
							}

                        }else{
                            
				if (defined $debugHash{$fieldCom."_".$field}){
					#Do nothing
				}else{
			    		$debugHash{$fieldCom."_".$field} = 1;
			        	print "Undefined field (typo error in config file or absent in annotation file (ex: exomiser):  ".$fieldCom." in ".$field."\n";
				}
				
                            $dataCommentHash{$field}{'values'} .= "<br><span class=\"commentTitle\">".$fieldCom . " :</span> NA";
                        }
                        #print $field.":\t".$fieldCom.":\t".$dataCommentHash{'Gene name'}{'values'}."\n";
				    }
			    }
				
				# add ranking decision as comment of AnnotSV ranking field
				#if (defined $dataCommentHash{$field}{'SVrank'}){
	  			#	if (defined $SVrankHash{$dataHash{'AnnotSV ID'}."_".$dataHash{'Gene name'}}){
                #		if (! defined $dataCommentHash{$field}{'values'}){
				#			$dataCommentHash{$field}{'values'}= "<br><b>Ranking :</b> " . $SVrankHash{$dataHash{'AnnotSV ID'}."_".$dataHash{'Gene name'}} ;
				#		}else{
				#			$dataCommentHash{$field}{'values'} = "<br><b>Ranking :</b> " .$SVrankHash{$dataHash{'AnnotSV ID'}."_".$dataHash{'Gene name'}} . $dataCommentHash{$field}{'values'};  	
				#		}
				#	}
				#}
			}
		}


		my $correctField;

		#fill finalSortData array   (try to invert foreach with NameColHash for external name)
		foreach my $field (keys %dataHash){
				
			if (defined $NameColHash{$field} && $NameColHash{$field} != 0){
				
				if ( $field =~ /^[PB]_/){

					if ($SV_type eq "DEL" ){ 
						$correctField = $field =~ s/gain|ins|inv/loss/r;
						$finalSortData[$NameColHash{$field} - 1] = $dataHash{$correctField};
					}elsif ($SV_type eq "DUP" ){
						$correctField = $field =~ s/loss|ins|inv/gain/r;
						$finalSortData[$NameColHash{$field} - 1] = $dataHash{$correctField};
					}else{
						$finalSortData[$NameColHash{$field} - 1] = $dataHash{$field};
					}
				}else{
					$finalSortData[$NameColHash{$field} - 1] = $dataHash{$field};
            	}
			}
		}

		#check missing fields (typo or error) and fill finalSortData array
		foreach my $field (keys %NameColHash){
			
			if (! defined $dataHash{$field}){
				$finalSortData[$NameColHash{$field} - 1] = "NA";
			}
		}



		#fill color LOEUF bin (old=pLI) 
		foreach my $field (keys %dataHash){
			#if (defined $dataHash{'pLI_ExAC'} ){
			if (defined $dataHash{'LOEUF_bin'} ){
                foreach my $bin (sort {$b <=> $a} keys %LOEUF_ColorHash){
                    if ( $dataHash{'LOEUF_bin'} eq "."){
                        $LOEUF_Color = '#FFFFFF';
                        last;
                    }
                    if ( $dataHash{'LOEUF_bin'}  >= $bin){
                        $LOEUF_Color = $LOEUF_ColorHash{$bin};
                        last;
                    }
                }
			}else{
                $LOEUF_Color = '#FFFFFF';            
            }
		}
        
		#change LOEUF color if ACMG gene
		if (defined $ACMGgene{$dataHash{'Gene_name'}}){
			$LOEUF_Color = '#808080';
		}



		# check if column names and data have the same size
		#if (scalar @finalSortData != keys %OutColHash){
		if (scalar @finalSortData != $OutColCounter){
			print "\nError in config file: it seems that some POSITION are missing.\nPlease correct config file with continuous 'POSITION' number.\n\n"; 
			exit 1;
		}

		
		#	print Dumper(\@finalSortData);

#############################################################################
######### FILL HASH STRUCTURE WITH COMMENTS 
		
		    
            $variantID = $line[0];

			#key for good sorting should be rank concatenated with annotSV_ID
			#
			#grant full line to be first
			if ($dataHash{"Annotation_mode"} eq "full"){
				#TODO compute CNV penalty
				$scorePenalty = 0;	
				if(defined $dataHash{"ACMG_class"}){
					
						if($dataHash{"ACMG_class"} ne "." && $dataHash{"ACMG_class"} ne "NA" ){
							$scorePenalty = $dataHash{"ACMG_class"}."_" ;
						}else{
							$scorePenalty .= "_" ;
						}

						if($SV_type eq "DEL"){
							$scorePenalty .= "3_";
						}elsif($SV_type eq "DUP"){
							$scorePenalty .= "2_";
						}else{
							$scorePenalty .= "1_";
						}

						if($dataHash{"Gene_count"} ne "."){
							$scorePenalty .= $dataHash{"Gene_count"}."_";
						}else{
							$scorePenalty .= "0_";
						}	
						if (defined $dataHash{"Exomiser_gene_pheno_score"} && $dataHash{"Exomiser_gene_pheno_score"} ne "-1.0" && $dataHash{"Exomiser_gene_pheno_score"} ne "NA"){
							$scorePenalty .= $dataHash{"Exomiser_gene_pheno_score"}."_";
						}else {
							$scorePenalty .= "0.0000_";
						}
						if(defined $dataHash{"OMIM_morbid"} && $dataHash{"OMIM_morbid"} eq "yes"){
							$scorePenalty .= "1_" ;
						}else{
							$scorePenalty .= "0_" ;
						}
					
				}
				$SV_ID{$dataHash{'AnnotSV_ID'}}{'finalScore'} = $scorePenalty;
				$fullSplitScore = 1000;
			}else{
				#TODO compute gene penalty exomiser > OMIM_morbid > LOEUF_bin
				$fullSplitScore = 10;
				if (defined $dataHash{"Exomiser_gene_pheno_score"} && $dataHash{"Exomiser_gene_pheno_score"} ne "-1" && $dataHash{"Exomiser_gene_pheno_score"} ne "NA"){
					$fullSplitScore += (20 * $dataHash{"Exomiser_gene_pheno_score"});	
				}	
				if(defined $dataHash{"OMIM_morbid"} && $dataHash{"OMIM_morbid"} eq "yes"){
					$fullSplitScore += 10;
				}
				if(defined $dataHash{"LOEUF_bin"} && $dataHash{"LOEUF_bin"} ne "."){
					$fullSplitScore +=  (10 - $dataHash{"LOEUF_bin"}) ; 		
				}

			}

			#finalsortData assigment according to full or split
			$hashFinalSortData{$SV_ID{$dataHash{'AnnotSV_ID'}}{'finalScore'}."_".$variantID}{$dataHash{'AnnotSV_ID'}}{$fullSplitScore}{$count}{'finalArray'} = [@finalSortData] ; 


			
			#url to UCSC for SV , highlight in blue and zoomout x1.5
			$hashFinalSortData{$SV_ID{$dataHash{'AnnotSV_ID'}}{'finalScore'}."_".$variantID}{$dataHash{'AnnotSV_ID'}}{$fullSplitScore}{$count}{'url2UCSC'} = "\"http://genome.ucsc.edu/cgi-bin/hgTracks?db=".$genomeBuild."&position=chr".$dataHash{"SV_chrom"}.":".$dataHash{"SV_start"}."-".$dataHash{"SV_end"}."&hgt.out1=submit&highlight=".$genomeBuild.".chr".$dataHash{"SV_chrom"}.":".$dataHash{"SV_start"}."-".$dataHash{"SV_end"}."#aaedff\" target=\"_blank\" rel=\"noopener noreferrer\"" ; 

			#URL to HUGO HGNC for split line only
			$hashFinalSortData{$SV_ID{$dataHash{'AnnotSV_ID'}}{'finalScore'}."_".$variantID}{$dataHash{'AnnotSV_ID'}}{$fullSplitScore}{$count}{'GeneName_short'} = $GeneName_link_string ; 

			#url to OMIM
			#$hashFinalSortData{$SV_ID{$dataHash{'AnnotSV_ID'}}{'finalScore'}."_".$variantID}{$dataHash{'AnnotSV_ID'}}{$fullSplitScore}{$count}{'url2OMIM'} = "https://www.omim.org/entry/".$dataHash{"OMIM_ID"}  ;

			$hashFinalSortData{$SV_ID{$dataHash{'AnnotSV_ID'}}{'finalScore'}."_".$variantID}{$dataHash{'AnnotSV_ID'}}{$fullSplitScore}{$count}{'OMIM_ID_short'} = $OMIM_ID_link_string ; 
			$hashFinalSortData{$SV_ID{$dataHash{'AnnotSV_ID'}}{'finalScore'}."_".$variantID}{$dataHash{'AnnotSV_ID'}}{$fullSplitScore}{$count}{'OMIM_phen_short'} = $OMIM_phen_link_string ; 


			#comment assigment
			foreach my $fieldCom (keys %dataCommentHash){
				$hashFinalSortData{$SV_ID{$dataHash{'AnnotSV_ID'}}{'finalScore'}."_".$variantID}{$dataHash{'AnnotSV_ID'}}{$fullSplitScore}{$count}{'hashComments'}{$NameColHash{$fieldCom} -1} = $dataCommentHash{$fieldCom}{'values'} ; 
				#print $dataCommentHash{'Gene name'}{'values'}."\n\n";
				#TODO check color for gene according to output position of gene field 
			}


            #assign color to Gene Name
			$hashFinalSortData{$SV_ID{$dataHash{'AnnotSV_ID'}}{'finalScore'}."_".$variantID}{$dataHash{'AnnotSV_ID'}}{$fullSplitScore}{$count}{'hashColor'}{$NameColHash{'Gene_name'}-1} = $LOEUF_Color ;



            #assign color to Gene Name
			$hashFinalSortData{$SV_ID{$dataHash{'AnnotSV_ID'}}{'finalScore'}."_".$variantID}{$dataHash{'AnnotSV_ID'}}{$fullSplitScore}{$count}{'fullRowColor'} = $fullRowColor ;




	} #END of IF-ELSE(#CHROM)	

}#END WHILE VCF


#prepare path to lib css and js


#original url used
#'https://code.jquery.com/jquery-3.5.1.js'
#'https://cdn.datatables.net/1.10.22/js/jquery.dataTables.min.js'
#'https://cdn.datatables.net/fixedheader/3.1.7/js/dataTables.fixedHeader.min.js'
#'https://cdn.datatables.net/fixedcolumns/3.3.2/js/dataTables.fixedColumns.min.js'

#'https://cdn.datatables.net/1.10.22/css/jquery.dataTables.min.css'
#'https://cdn.datatables.net/fixedheader/3.1.7/css/fixedHeader.dataTables.min.css'
#'https://cdn.datatables.net/fixedcolumns/3.3.2/css/fixedColumns.dataTables.min.css'


my $path2jquery = "https://code.jquery.com/";
my $path2jqueryDT = "https://cdn.datatables.net/1.10.22/js/";
my $path2jsFCDT = "https://cdn.datatables.net/fixedcolumns/3.3.2/js/";

my $path2css = "https://cdn.datatables.net/1.10.22/css/"; 
my $path2FCcss = "https://cdn.datatables.net/fixedcolumns/3.3.2/css/";

if ($datatableDir ne ""){

	#remove trailing slash from path
	$datatableDir =~ s/\/$//;

	$path2jquery = $datatableDir."/js/";
	$path2jqueryDT = $datatableDir."/js/";

	$path2css = $datatableDir."/css/"; 
	$path2FCcss = $datatableDir."/css/";
}





############ THE TIME HAS COME TO FILL THE HTML OUTPUT FILE
####################################################################
######################### HTML  Initialisation ######################
#
#
my $cleanOutBasename = $outBasename =~ s/\.//rg;

#file header
my $htmlStart = "<!DOCTYPE html>\n<html>
\n<head>
\n<meta charset=\"utf-8\">
\n<title>".$outPrefix.$outBasename."</title>\n
\n<link rel=\"shortcut icon\" href=\"https://github.com/mobidic/knotAnnotSV/raw/master/images/knot_favicon.png\" type=\"image/png\"/>

\n<link rel=\"stylesheet\" href=\"/home/tg/TOOLS/w2ui-1.5.rc1.css\" type=\"text/css\"/>

\n<style>
\n.loeufbin0{
\n\tcolor: ".$LOEUF_ColorHash{'0.0'}."	
\n}
\n.loeufbin1{
\n\tcolor: ".$LOEUF_ColorHash{'1.0'}."	
\n}
\n.loeufbin2{
\n\tcolor: ".$LOEUF_ColorHash{'2.0'}."	
\n}
\n.loeufbin3{
\n\tcolor: ".$LOEUF_ColorHash{'3.0'}."	
\n}
\n.loeufbin4{
\n\tcolor: ".$LOEUF_ColorHash{'4.0'}."	
\n}
\n.loeufbin5{
\n\tcolor: ".$LOEUF_ColorHash{'5.0'}."	
\n}
\n.loeufbin6{
\n\tcolor: ".$LOEUF_ColorHash{'6.0'}."	
\n}
\n.loeufbin7{
\n\tcolor: ".$LOEUF_ColorHash{'7.0'}."	
\n}
\n.loeufbin8{
\n\tcolor: ".$LOEUF_ColorHash{'8.0'}."	
\n}
\n.loeufbin9{
\n\tbackground-color: ".$LOEUF_ColorHash{'9.0'}.";	
\n}
\n</style>
\n<script src=\"http://ajax.googleapis.com/ajax/libs/jquery/2.1.1/jquery.min.js\"></script>
\n<script type=\"text/javascript\" src=\"/home/tg/TOOLS/w2ui-1.5.rc1.min.js\"></script>


\n</head>
\n<body>
\n<div id=\"grid\" style=\"width: 100%; height: 400px;\"></div>
\n</body>

\n<script type=\"text/javascript\">
\$(function () {
	\$('#grid').w2grid({ 
		name: 'grid', 
		columns: [";

                

foreach my $col (sort {$a <=> $b} keys %OutColHash){
			#print HTML "\t<th style=\"word-wrap: break-word\"   >";
	print "colNbr:   ".$col."\t".$OutColHash{$col}{'field'}."\n";  
	if (defined $OutColHash{$col}{'field'} && $col != 0){
		if (defined $OutColHash{$col}{'field'} && $OutColHash{$col}{'field'} eq "Gene_name"){
			$htmlStart .= "\t\t\t{field: \"".$OutColHash{$col}{'field'}."\", caption: \"".$OutColHash{$col}{'RENAME'}."\", size: '30%'},\n";
		}elsif (defined $OutColHash{$col}{'field'} && $OutColHash{$col}{'field'} eq "AnnotSV_ID"){
			$htmlStart .= "\t\t\t{field: \"".$OutColHash{$col}{'field'}."\", caption: \"TOTO\", size: '30%'},\n";
			#$htmlStart .= "{field: \"".$OutColHash{$col}{'field'}."\", caption: \"".$OutColHash{$col}{'RENAME'}."\", render: function (record, index, column_index) { 
			#  var html =  record.AnnotSV_ID;
			#  return html;
			#}  },\n";
		   
		}else{
			$htmlStart .= "\t\t\t{field: \"".$OutColHash{$col}{'field'}."\", caption: \"".$OutColHash{$col}{'RENAME'}."\", size: '30%'},\n";	
		}
	}
}
$htmlStart .= "],\n\n";
        


$htmlStart .= "\t\trecords: [\n";

my $kindRank=0;

# check if last line was "FULL" to finish the raw with </tr>


foreach my $rank (rnatkeysort { "$_-$hashFinalSortData{$_}" } keys %hashFinalSortData){
	print $rank."\n";
	
	foreach my $ID ( keys %{$hashFinalSortData{$rank}}){
		
		#foreach my $rankSplit (sort {$hashFinalSortData{$rank}{$ID}{$b} <=> $hashFinalSortData{$rank}{$ID}{$a} } ( keys %{$hashFinalSortData{$rank}{$ID}})){
		foreach my $rankSplit (rnatkeysort { "$_-$hashFinalSortData{$rank}{$ID}{$_}" }  keys %{$hashFinalSortData{$rank}{$ID}}){
	
			#print "DEBUG ranksplit:  ".$rankSplit."\n";

			foreach my $variant ( keys %{$hashFinalSortData{$rank}{$ID}{$rankSplit}}){
			
				#increse rank number then change final array
				$kindRank++;

				#FILL tab 'ALL';
				#if (    $hashFinalSortData{$rank}{$ID}{$rankSplit}{$variant}{'finalArray'}[$NameColHash{'Annotation_mode'} - 1] eq "full") {
			
				$htmlStart .= "\t\t\t{recid: ".$kindRank.", ";

				#Once for "ALL"  = FULL+SPLIT
				for( my $fieldNbr = 0 ; $fieldNbr < scalar @{$hashFinalSortData{$rank}{$ID}{$rankSplit}{$variant}{'finalArray'}} ; $fieldNbr++){
					#print $fieldNbr."\n";
		

					if (defined $hashFinalSortData{$rank}{$ID}{$rankSplit}{$variant}{'finalArray'}[$fieldNbr]){
						if ($fieldNbr eq $NameColHash{'AnnotSV_ID'} - 1){
							$htmlStart .= "".$OutColHash{$fieldNbr + 1 }{'field'}.": '"."<a href=".$hashFinalSortData{$rank}{$ID}{$rankSplit}{$variant}{'url2UCSC'}.">".$hashFinalSortData{$rank}{$ID}{$rankSplit}{$variant}{'finalArray'}[$fieldNbr]."</a>',";
						}else{
							$htmlStart .= "".$OutColHash{$fieldNbr + 1 }{'field'}.": '".$hashFinalSortData{$rank}{$ID}{$rankSplit}{$variant}{'finalArray'}[$fieldNbr]."',";
						}
					}else{
							$htmlStart .= "".$OutColHash{$fieldNbr + 1 }{'field'}.": '.',";
					}

					
				} #END FOR FULL+SPLIT tab
				
				$htmlStart .= "LOEUF_color: '".$hashFinalSortData{$rank}{$ID}{$rankSplit}{$variant}{'hashColor'}{$NameColHash{'Gene_name'} -1}."',";
				if (    $hashFinalSortData{$rank}{$ID}{$rankSplit}{$variant}{'finalArray'}[$NameColHash{'Annotation_mode'} - 1] eq "full") {
					$htmlStart .=  "w2ui: { style: \"background-color:".$hashFinalSortData{$rank}{$ID}{$rankSplit}{$variant}{'fullRowColor'}  .";\"}},\n"; 
				}else{
					$htmlStart .=  "w2ui: { style: {".($NameColHash{'Gene_name'}-1).": 'background-color: ".$hashFinalSortData{$rank}{$ID}{$rankSplit}{$variant}{'hashColor'}{$NameColHash{'Gene_name'} -1}.";'}}},\n"; 

				} 
				#$htmlStart .= "};\n";

			}   # END FOREACH VARIANT


		}   # END FOREACH RANKSPLIT
	}   # END FOREACH ID
}	#END FOREACH RANK



        

$htmlStart .= "],
    });    
});
\n</script>
</html>";


close(VCF);

open(HTML, '>', $outDir."/".$outPrefix.$outBasename.".html") or die $!;

#Write in HTML file
print HTML $htmlStart;

close(HTML);

print STDERR "\n\n\nDone!\n\n\n";


exit 0; 


