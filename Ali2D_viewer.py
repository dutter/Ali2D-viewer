#!/usr/bin/env python3

import argparse
import re
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import pandas as pd
import warnings

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",
                    help="Ali2D text output [underscores in names will be turned into spaces]")
parser.add_argument("-o","--output",
                    help="Where to save the .pdf image [default: INPUT.pdf from INPUT.txt]", default="")
parser.add_argument("-w","--wrap",
                    help="How many positions to print per line before wrapping [default: 180]", default=180, type=int)
parser.add_argument("-W", "--width", help="Width of figure to save [default: 12]", default=12, type=float)
parser.add_argument("-H", "--height", help="Height of figure to save [default: 8]", default=8, type=float)
parser.add_argument("-s","--size_aa",
                    help="Font size in pt for position in alignment [default: 5]", default=5)
parser.add_argument("-S","--size_name",
                    help="Font size in pt for names of sequences [default: 7]", default=7)
parser.add_argument("--color_helix", default='0.9,0.4,0.4',
                    help="Comma-separated RGB color for max alpha helix confidence, default is '0.9,0.4,0.4'")
parser.add_argument("--color_sheet", default='0.4,0.4,0.9',
                    help="Comma-separated RGB color for max beta sheet confidence, default is '0.3,0.3,0.9'")

args=parser.parse_args()

if args.output is "" and re.search(r'\....$', args.input):
    args.output = re.sub(r'\....$', '.pdf', args.input)



with open(args.input) as f:
    d= dict(line.strip().split(':')[0:2] for line in f)

aln = pd.DataFrame.from_dict(d, orient = 'index', columns=None)

# Figure out key sequence dimensions
wrap = args.wrap
nResidues = len(aln.loc['seq0',0])
nWraps = int(np.ceil(nResidues/wrap))


#  GET THE BITS TO PLOT, FULL LENGTH
names = aln.loc[[x for x in list(aln.index.values) if re.search('name',x)],0]
names = [re.sub('_',' ', x) for x in names]


conf = aln.loc[[x for x in list(aln.index.values) if re.search(r"^confi", x )],0]
conf = [list(re.sub(' ','', re.sub('-','0', x))) for x in conf]
conf = np.array(conf, dtype='float32').__truediv__(9)  # rescale 0-9:0-1 for alpha

residues = aln.loc[[x for x in list(aln.index.values) if re.search(r"^seq[a-z0-9]*$",x)],0]
residues = np.array([list(re.sub(" ", "", x)) for x in residues])

structure = aln.loc[[x for x in list(aln.index.values) if re.search(r"^seq_struct", x )],0]

# make structure dict to convert the structure type code to the corresponding color value, then apply to each character
# 3rd dimension of length 3 is now RGB value for each residue
structDict={'C':[1,1,1],
            'H':[float(x) for x in args.color_helix.split(',')], 'E':[float(x) for x in args.color_sheet.split(',')]}
structure = [list(map(structDict.get, list(re.sub('-','C', re.sub(' ','', x))))) for x in structure]

# figure out how long many positions are in the alignment
nSeq = len(conf[:,0])

# reshape so is position x sequence x RGB_of_2ndary_type
structure = np.array(structure, dtype='float32').reshape((nSeq, nResidues - 1, 3))

fig, ax = plt.subplots(nWraps, figsize=(args.width,args.height))

for chunk in range(nWraps):

    # subset to be the window of each wrap
    subStruct = structure[:,0:(wrap*chunk + wrap),:]
    subConf = np.expand_dims(conf[:,0:(wrap*chunk + wrap)], axis=2)

    subConf = np.concatenate((subStruct, subConf), axis=2)

    # plot the color squares
    ax[chunk].set_xlim(left=-0.5 + wrap*chunk, right=wrap + wrap*chunk)
    im = ax[chunk].imshow(subConf, aspect='auto')

    # change y axis labels from indices to corresponding names
    ax[chunk].set_yticks(ticks=np.linspace(0,len(names),len(names)+1), minor=False)
    ax[chunk].set_yticklabels(names, fontsize=args.size_name)
    ax[chunk].tick_params(axis='y', length=0)

    ax[chunk].spines['top'].set_visible(False)
    ax[chunk].spines['right'].set_visible(False)
    ax[chunk].spines['bottom'].set_visible(True) # keep the x axis line
    ax[chunk].spines['left'].set_visible(False)

    # overlay the alignment text
    for y in range(nSeq):
        for x in range(chunk*wrap, chunk*wrap + wrap):
            try:
                ax[chunk].text(x, y, residues[y, min(x, nResidues)],
                               ha='center', va='center', color='k', fontsize = args.size_aa)
            except IndexError:
                continue

acb = fig.colorbar(
    plt.cm.ScalarMappable(norm=colors.Normalize(vmin=0,vmax=9),
                          cmap=colors.ListedColormap(
                              np.linspace((1,1,1), tuple(float(x) for x in args.color_helix.split(',')),10))),
             fig.add_axes([0.91,0.15,0.01,0.7]))
bcb = fig.colorbar(
    plt.cm.ScalarMappable(norm=colors.Normalize(vmin=0,vmax=9),
                          cmap=colors.ListedColormap(
                              np.linspace((1,1,1),tuple(float(x) for x in args.color_sheet.split(',')),10))),
             fig.add_axes([0.96,0.15,0.01,0.7]))

# set ticks/labels to the middle of their bar, not bottom:top
tick_locs = (np.arange(10) + 0.5)*(10-1)/10
acb.set_ticks(tick_locs)
bcb.set_ticks(tick_locs)
acb.set_ticklabels(np.arange(10))
bcb.set_ticklabels(np.arange(10))

#titles
acb.ax.set_title(r'$\alpha$')
bcb.ax.set_title(r'$\beta$')

with warnings.catch_warnings():
    warnings.simplefilter("ignore", UserWarning)
    plt.tight_layout(pad=1, h_pad=0.2, rect=(0,0,0.91,1))

plt.savefig(fname=args.output)





