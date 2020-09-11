#!/bin/python3

from phrank import Phrank
from phrank import utils as phrank_utils


DAG="hpodag_current.txt"
DISEASE_TO_PHENO="disease_to_pheno.current.txt"
DISEASE_TO_GENE="gene_to_disease.current.txt"
GENE_TO_PHENO="data/gene_to_pheno.amelie.txt"
p_hpo = Phrank(DAG, diseaseannotationsfile=DISEASE_TO_PHENO, diseasegenefile=DISEASE_TO_GENE)

#Fill gene list (from genemap2)
with open('genemap2_genelist.txt') as f:
	genes= f.read().splitlines()

patient_genes = set(genes)

# defining the phenotype sets

with open('HPO.txt') as h:
	phenotypeset1 = h.read().splitlines()



#phenotypeset1 = ['HP:0000077','HP:0030765','HP:0012115','HP:0002088','HP:0002099','HP:0001945','HP:0000719']
#phenotypeset2 = ['HP:0000975','HP:0002018','HP:0000421','HP:0012393','HP:0004406','HP:0002321']

# computing the similarity between two sets of phenotypes
#matchscore = p_hpo.compute_phenotype_match(phenotypeset1, phenotypeset2)
#print ("the phenotype similarity score is %.2f"%matchscore)


# defining patient genes and phenotypes
#patient_genes = set(['ENSG00000000419','ENSG00000000971','ENSG00000000971','ENSG00000001626','ENSG00000001626','ENSG00000001631','ENSG00000002822','ENSG00000003137'])
patient_phenotypes = phenotypeset1


#open output disease file
fd = open("output_disease_Phrank.txt", "w")
#open output gene file
fg = open("output_gene_Phrank.txt", "w")




# sorting the disease by best match
disease_ranking = p_hpo.rank_diseases(patient_genes, patient_phenotypes)
print ("Disease ranking...")
for disease_info in disease_ranking:
	#print ("disease id: %s\tsimilarity score: %.2f"%(disease_info[1],disease_info[0]))
	fd.write("disease id: %s\tsimilarity score: %.2f\n"%(disease_info[1],disease_info[0]))




# sorting the genes by best match
gene_ranking = p_hpo.rank_genes(patient_genes, patient_phenotypes)
print ("\nGene ranking...")
for gene_info in gene_ranking:
	#print ("Gene symbol: %s\tsimilarity score: %.2f"%(gene_info[1],gene_info[0]))
	fg.write("Gene symbol: %s\tsimilarity score: %.2f\n"%(gene_info[1],gene_info[0]))





fd.close()
fg.close()



print ("Done!")


