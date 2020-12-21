import argparse
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('input',metavar='input_dir',
                    help='Directory and basename of the report files')
parser.add_argument('output',metavar='output_dir',
                    help='Directory to save the heatmaps')
args = parser.parse_args()

print('Running with args: ')
print('\tinput:\t'+args.input)
print('\toutput:\t'+args.output+'\n')

for report in glob.glob(args.input+'*.txt'):
    if report.startswith('.'):
        continue
    print('Found '+report)
    df = pd.read_csv(report,sep='\t',usecols=['POS2','POS1','LD'])

    x = df['POS2']
    y = df['POS1']
    z = df['LD']

    df = df.set_index(['POS2','POS1'])
    #print(df.index)
    df = df.unstack()
    #df = df.fillna(0)
    df.columns = df.columns.get_level_values(1)
    print(df)

    mask = np.triu(np.ones_like(df, dtype=bool))
    n = 5
    color = plt.cm.get_cmap("tab20b_r", n)

    fig, ax = plt.subplots(figsize=(4, 3), dpi= 300)
    fig.tight_layout()
    ax.autoscale()
    ax.tick_params(axis='both', which='major', labelsize=4)
    sns.heatmap(df, cmap=color, vmin=0, vmax=1, square=True, cbar_kws={"shrink": .8})
    #sns.heatmap(df, mask=mask, cmap=color, vmin=0, vmax=1, square=True, cbar_kws={"shrink": .8})
    plt.xlabel('Position 1')
    plt.ylabel('Position 2')
    ax.set_title("LD scores\nin "+report, fontsize=11)
    plt.draw()
    plt.savefig(report.replace('.txt','.png'), bbox_inches='tight')
print('Ended')
