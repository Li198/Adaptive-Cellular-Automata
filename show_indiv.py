import numpy as np
import matplotlib.pyplot as plt
from individual import individual 



if __name__ == '__main__':
    gene_pool_initial = np.loadtxt('gene_pool4')
    gene_pool = gene_pool_initial.astype(int)

    indiv_index = 36

    gene = gene_pool[indiv_index]

    indiv = individual(gene)
    indiv.show_cells()

        



