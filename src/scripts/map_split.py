import os, fire
# 输入 sam, 输出 unique_sam
from plot import mapping_ratio_stack

def parse_tags(s):
  tags = {}
  for tag in s.split('\t'):
    name, type_, val = tag.split(':')
    tags[name] = int(val) if type_ == 'i' else val
  return tags

def draw_stack_plot(steps, stackbar_list, output_root):
  x_un, x_uni, x_mul = stackbar_list
  plt.figure(figsize=(10,7)) #设置画布的尺寸
  plt.bar(steps, x_un, label="unmapped reads",edgecolor = 'black')
  plt.bar(steps, x_uni, label="unique mapped reads",edgecolor = 'black', bottom = x_un)
  plt.bar(steps, x_mul, label="non-unique mapped reads",edgecolor = 'black', bottom = [i+j for i, j in zip(x_un,x_uni)])
  plt.legend( loc=3,fontsize=14)  # 设置图例位置
  plt.ylabel('Mapping reads number',fontsize=16)
  plt.xlabel('Mapping reads split method',fontsize=16)

  for x1, y1, y2, y3 in zip(steps, x_un, x_uni, x_mul):
    p1 = y1/(y1+y2+y3)
    p2 = y2/(y1+y2+y3)
    p3 = y3/(y1+y2+y3)
    if p1 >0.05:
      plt.text(x1, y1 * 0.4, '{:.0%}'.format(p1), ha='center',fontsize = 15)
    if p2>0.05:
      plt.text(x1, y1 + (y2)* 0.4,  '{:.0%}'.format(p2), ha='center',fontsize = 15)
    if p3>0.05:
      plt.text(x1, y1 + y2 + (y3)* 0.4, '{:.0%}'.format(p3), ha='center',fontsize = 15)

  plt.savefig(output_root +'imgs/mapping_split.png')
  plt.savefig(output_root +'imgs/mapping_split.pdf')


def main(sam_dir, sam_file, output_root):
  # unmapped: SAM flag includes 4
  # 简单来说: 区分unique和nonredundant，nonredundant是绝对的单一比对，而unique是可以有多比对，但最高分就一个。
  # AS是比对的最高分，XS是第二高分。
  # uniquemapped: include multihit but uniquely determinable reads by alignment score
  # multimapped: multiple aligned reads with same scores
  # nonredundant: aligned exactly 1 time (bowtie2)
  # redundant: aligned >1 times(bowtie2), include uniquely determinable reads by alignment score
  total,unmapped, uniquemapped, multimapped, nonredundant, redundant  = 0, 0, 0, 0, 0, 0
  unmapped_sam = open(sam_dir + 'unmapped.sam', 'w')
  uniquemapped_sam = open(sam_dir + 'uniquemapped.sam', 'w')
  multimapped_sam = open(sam_dir + 'multimapped.sam', 'w')
  nonredundant_sam = open(sam_dir + 'nonredundant.sam', 'w')
  redundant_sam = open(sam_dir + 'redundant.sam', 'w')
  for ln in open(sam_dir + sam_file):
    if ln.startswith('@'):
      unmapped_sam.write(ln)
      uniquemapped_sam.write(ln)
      multimapped_sam.write(ln)
      nonredundant_sam.write(ln)
      redundant_sam.write(ln)
      continue
    _, flag, _, _, _, _, _, _, _, _, _, tags = ln.split('\t', 11)
    total += 1
    if int(flag) & 4:
      unmapped += 1
      unmapped_sam.write(ln)
    else:
      tags = parse_tags(tags)
      if "AS" not in tags.keys():
        continue
      as_ = tags["AS"]
      xs = tags.get("XS")
      if as_ != xs:
        uniquemapped += 1
        uniquemapped_sam.write(ln)
      else:
        multimapped += 1
        multimapped_sam.write(ln)
      if xs is None or (as_ != 0 and xs == 0):
        nonredundant += 1
        nonredundant_sam.write(ln)
      else:
        redundant += 1
        redundant_sam.write(ln)
  unmapped_sam.close()
  uniquemapped_sam.close()
  multimapped_sam.close()
  nonredundant_sam.close()
  redundant_sam.close()
  print('    unmapped--'+str(unmapped) )
  print('    uniquemapped--'+str(uniquemapped) )
  print('    multimapped--'+str(multimapped) )
  print('    nonredundant--'+str(nonredundant) )
  print('    redundant--'+str(redundant) )
  csv_dir = sam_dir.replace('/sam/', '/csv/mapping_report.csv')
  with open(csv_dir, 'w') as f:
    f.write( 'total_mapped,' + str(total)+'\n' )
    f.write( 'unmapped,'+ str(unmapped)+'\n' )
    f.write( 'unique_mapped,' + str(uniquemapped)+'\n' )
    f.write( 'multiple_mapped,' + str(multimapped)+'\n'  )
    f.write( 'nonredundant_mapped,' + str(nonredundant)+'\n')
    f.write( 'redundant_mapped,' + str(redundant)+'\n')

  stackbar_list = [
    [unmapped, unmapped ],
    [uniquemapped, nonredundant],
    [multimapped, redundant]
  ]
  mapping_ratio_data = {'map_ratio_list': stackbar_list,
  'soft': ['bwa-mem', 'bowtie2'] }
  mapping_ratio_stack(mapping_ratio_data, output_root)

if __name__ == '__main__':
  fire.Fire( main )
