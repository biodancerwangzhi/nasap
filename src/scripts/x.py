import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from parse_gtf import get_gene_df

def scatterplot( x_list, y_list, title, x_label, y_label, text, output_file ):
  fig = plt.figure(figsize=(12, 7))
  ax = plt.gca()
  plt.scatter(x=x_list, y=y_list)
  plt.text(0.5, 0.75, text, transform = ax.transAxes)
  plt.title( title )
  plt.xlabel(x_label)
  plt.ylabel(y_label)
  fig.savefig(output_file)
  plt.close(fig)

output_root = '../../nasap_GSM3309956/'
gtf = '../../../data/Homo_sapiens.GRCh38.93.gtf'


print( protein_data.head() )
filter_genes = protein_data[pd.to_numeric( protein_data.pi, errors='coerce' ) > 2]


gene_df = get_gene_df( gtf )
print( gene_df.head )
filter_gene_df = pd.DataFrame({})
forward_bw_signal = pyBigWig.open(  forward_bw )
reverse_bw_signal = pyBigWig.open(  reverse_bw )

