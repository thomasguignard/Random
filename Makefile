get-phenotype-data:

	#get gene_to_phenotype from http://compbio.charite.de/jenkins/job/hpo.annotations/lastBuild/
	#wget http://compbio.charite.de/jenkins/job/hpo.annotations/lastBuild/artifact/util/annotation/genes_to_phenotype.txt
	wget http://compbio.charite.de/jenkins/job/hpo.annotations/lastSuccessfulBuild/artifact/util/annotation/genes_to_phenotype.txt

	#extract data to get right format
	awk -F "\t" '{print $2"\t"$9}' genes_to_phenotype.txt | sort | uniq > gene_to_disease.build1273.txt
	
	#get hpoa from  
	wget http://compbio.charite.de/jenkins/job/hpo.annotations.current/lastSuccessfulBuild/artifact/current/phenotype.hpoa
	# extract and format
	awk -F "\t" '$1!~/#/{print $4"\t"$1}'  phenotype.hpoa |sort -u > disease_to_pheno.hpoa.txt

	#merge to get final 
	cat gene_to_disease.build1273.txt > gene_to_disease.hpoa_build1273.txt
	cat disease_to_pheno.hpoa.txt >> gene_to_disease.hpoa_build1273.txt

	sort -u gene_to_disease.hpoa_build1273.txt > gene_to_disease.hpoa_build1273_uniq.txt



get-hp.obo-to-create-hpo-dag:

	#download
	wget	https://github.com/obophenotype/human-phenotype-ontology/blob/master/hp.obo

	# create hpo dag
	awk -F ": | !" 'BEGIN{ID="";ISA=""}{if($1=="id"){ID=$2}else if($1=="is_a"){print ID"\t"$2} }' hp.obo > hpodag_2020-08-11.txt


get-omim-genelist:

	awk -F "\t" '$9!="" && $1!~/#/{print $9}'  ~/DATABASES/OMIM/2020/genemap2.txt |sort -u > genemap2_genelist.txt
