import argparse
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

nColors = 5

parser = argparse.ArgumentParser()
parser.add_argument('input',metavar='input_dir',
                    help='Directory and basename of the report files')
parser.add_argument('output',metavar='output_dir',
                    help='Directory to save the heatmaps')
parser.add_argument('nColors',metavar='numberOfColors',
                    help='Number of color regions in heatmap',
                    nargs='?')
args = parser.parse_args()

if args.nColors:
    nColors = int(args.nColors)

print('qLD Heatmap creator running with arguments:')
print('\tinput:\t'+args.input)
print('\toutput:\t'+args.output)
print('\tcolors:\t'+str(nColors)+'\n')

ifile = os.path.splitext(args.input)[0]

for report in glob.glob(args.input+'*'):
    if report.startswith('.'):
        continue
    print('Found '+report)
    df = pd.read_csv(report,sep='\t',usecols=['POS2','POS1','LD'])

    x = df['POS2']
    y = df['POS1']
    z = df['LD']

    df = df.set_index(['POS2','POS1'])
    df = df.unstack()
    df.columns = df.columns.get_level_values(1)

    color = plt.cm.get_cmap("tab20b_r", nColors)

    fig, ax = plt.subplots(figsize=(4, 3), dpi= 300)
    plt.tight_layout()
    ax.autoscale()
    ax.tick_params(axis='both', which='major', labelsize=4)
    sns.heatmap(df, cmap=color, vmin=0, vmax=1, square=True, cbar_kws={"shrink": .8})
    plt.xlabel('Position 1')
    plt.ylabel('Position 2')
    ax.set_title("LD scores\nin "+report, fontsize=11)
    plt.draw()
    plt.savefig(report.replace('.txt','.png'), bbox_inches='tight')
print('Files saved in: '+os.path.dirname(args.input))
