import matplotlib.pyplot as plt
import matplotlib.colors as colors
color_list = list(colors._colors_full_map.values())

def boxplot( boxplot_data, title, x_label, y_label, xtick_labels, output_file ):
  # matplotlib比较old， 输入是 二维list
  fig = plt.figure(1, figsize=(12, 7))
  ax = fig.add_subplot(111)
  bp = ax.boxplot(boxplot_data, patch_artist=True)

  for box in bp['boxes']:
    box.set(facecolor='#087E8B', alpha=0.6, linewidth=2)

  for whisker in bp['whiskers']:
    whisker.set(linewidth=2)

  for median in bp['medians']:
    median.set(color='black', linewidth=3)

  ax.set_title( title )
  ax.set_xlabel( x_label)
  ax.set_ylabel( y_label )
  ax.set_xticklabels( xtick_labels )
  fig.savefig(output_file)
  plt.close(fig)

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

def lollipopplot(chrs, lollipop_dic, output_file):
  fig, axes = plt.subplots(len(chrs)*2, 1, sharex=True, sharey=False, figsize=(30, 2*len(chrs) ))
  plt.xticks([])
  plt.subplots_adjust(wspace =0, hspace =0)
  i = 0
  for chr in chrs:
    color = color_list[i]
    try:
      axes[i].stem( lollipop_dic[chr]['forward']['site'], lollipop_dic[chr]['forward']['count'], linefmt=color, markerfmt='.', basefmt=color  )
    except:
      pass
    axes[i].set_ylabel( chr+':+', color=color )
    axes[i].yaxis.tick_right()
    axes[i].yaxis.get_major_ticks()[1].set_visible(False)
    axes[i].yaxis.get_major_ticks()[-2].set_visible(False)
    axes[i].spines['top'].set_visible(False)
    axes[i].spines['left'].set_visible(False)
    axes[i].spines['right'].set_color( color )
    i+=1
    try:
      axes[i].stem( lollipop_dic[chr]['reverse']['site'], [-1*x for x in lollipop_dic[chr]['reverse']['count']], linefmt=color, markerfmt='.', basefmt=color )
    except:
      pass
    axes[i].set_ylabel( chr+':-', color=color)
    axes[i].yaxis.tick_right()
    axes[i].yaxis.get_major_ticks()[1].set_visible(False)
    axes[i].yaxis.get_major_ticks()[-2].set_visible(False)
    axes[i].spines['bottom'].set_visible(False)
    axes[i].spines['left'].set_visible(False)
    axes[i].spines['right'].set_color(color)
    i+=1


  plt.title('Genome wide pausing sites distribution', y=1.1, fontsize=18)
  fig.savefig(output_file)
  plt.close(fig)