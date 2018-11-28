__version__ = "0.5.1"
__author__ = "Forest Dussault"
__email__ = "forest.dussault@canada.ca"


def get_phenotype(gene: str):
    return gene


def pipeline(gene_list: list) -> list:
    phenotype_list = []
    for gene in gene_list:
        phenotype = get_phenotype(gene)
        phenotype_list.append(phenotype)
    return phenotype_list
