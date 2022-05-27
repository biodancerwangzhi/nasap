import pandas as pd
from collections import defaultdict
import itertools

def get_gene_df( gtf ):
  feature_list = []
  for ln in open(gtf):
    if ln.startswith('#'): continue
    fields = ln.split('\t')
    if fields[2] != 'gene': continue
    chrom, start, end, strand = fields[0], int( fields[3]), int( fields[4]), fields[6]
    attr_field = fields[-1]
    if not chrom.startswith('chr'): chrom = 'chr' + chrom

    # if 'gene_name "' not in attr_field: continue
    # if 'gene_biotype "' not in attr_field: continue
    for keywords in [ ['gene_type "', 'gene_name "'],  ['gene_biotype "', 'gene_name "'],  ['gene_biotype "', 'gene "']]:
      try:
        gene_name, genetype = attr_field.split(keywords[1])[1].split('"')[0], attr_field.split(keywords[0])[1].split('"')[0]
        feature_list.append( [chrom, start, end, strand, gene_name, genetype] )
        break
      except:
        pass

  df = pd.DataFrame(feature_list, columns =["chrom", "start", "end", "strand", "gene_name", "genetype"] )
  return df

# df  = get_gene_df( r'H:\tmpDownload\gencode.v38.annotation.gtf\gencode.v38.annotation.gtf' )

# from collections import Counter
# print( df.shape )
# print( Counter( df['biotype'] ) )


def to_ranges(iterable):
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable),
                                        lambda t: t[1] - t[0]):
        group = list(group)
        yield group[0][1], group[-1][1]

def get_gene_merge_exon_dic(gtf):
  dd = defaultdict( list )
  gene_chr_strand = {}
  for ln in open(gtf):
    if ln.startswith('#'): continue
    if 'gene_biotype "protein_coding"' not in ln: continue
    fields = ln.split('\t')
    if fields[2] != 'exon': continue
    start, end = int( fields[3] ), int( fields[4] )
    try:
      gene_name = fields[-1].split('gene_name "')[1].split('"')[0]
    except:
      continue
    dd[gene_name].append(  range(start, end)  )
    chrom, strand = fields[0], fields[6]
    if not chrom.startswith('chr'): chrom = 'chr'+ chrom
    gene_chr_strand[gene_name] = [chrom, strand]

  gene_exon_range_dic = {}
  for gene, range_list in dd.items():
    l = []
    for r in range_list:
      l.extend( list(r) )
    range_list = list( to_ranges(l) )
    if len( range_list) > 0:
      gene_exon_range_dic[gene] = {'chrom':gene_chr_strand[gene][0], 'strand':gene_chr_strand[gene][1], 'exon':range_list  }
  return gene_exon_range_dic

# d = get_gene_merge_exon_dic(r'../ensembl_Homo_sapiens.GRCh38.93.gtf')
# for k in list( d.keys() )[:5]:
#   print(k, d[k] )

def gene_df2dic(df):
  dic = {'+': {}, '-': {}}
  gene_name_list = list( df['gene_name'] )
  chr_list = list( df['chrom'] )
  start_list = list( df['start'] )
  end_list = list( df['end'] )
  strand_list = list( df['strand'] )
  for i, gene in enumerate(gene_name_list):
    strand = strand_list[i]
    dic[strand][gene] = [chr_list[i], start_list[i], end_list[i]]
  return dic