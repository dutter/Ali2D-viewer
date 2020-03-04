# Ali2D-viewer
A simple python script for pretty plotting of Ali2D text output (alignment + predicted 2-dimensional structure)

[Ali2D](https://toolkit.tuebingen.mpg.de/tools/ali2d) is an excellent webtool for predicting secondary structure (alpha helices, beta sheets) from an amino acid alignment. However, saving and reporting Ali2D results in an aesthetic way can be a bit challenging. This simple plotter aims to help speed up the generation of aesthetically-acceptable alignments from Ali2D by parsing and replotting the text output of an Ali2D run.

## Installation:
* Download Ali2D_viewer.py
* Put it somewhere in your $PATH variable; for instance, `/usr/local/bin/` can be an OK place on a Mac
* Or put it somewhere else and directly call it, like `~/Downloads/Ali2D_viewer.py --input ali2d_results.out ...`

## Basic usage:
#### Step 1: Download your Ali2D results
When the run is completed, the rightmost tab "Text Output" is the information you want. Click 'Download' in the pane below, or copy/paste the text output into a new plain text file.

#### Step 2: Run Ali2D_viewer.py with basic options
The basic Ali2D_viewer.py functionality is invoked with `Ali2D_viewer.py -i ali2d_results.out`. This will generate a plot, with default number of positions per line, and save it in a pdf file at `ali2d_results.pdf` with default dimensions.

Some useful flags for controlling the output:
* `--wrap N` sets N as the number of positions per line before wrapping to another chunk. Default is 180 residues/line
* `--width N` and `--height M` specify the width and height of the plot
* `--output OUTFILE` can specify the name of the saved pdf
* Other options for text size and so forth can be found in the help, `Ali2D_viewer.py --help`

#### Some notes for posterity
Underscores in sequence name are translated to spaces

As Ali2D uses the `:` character to separate the key from the object, colons should be avoided in sequence names

The red/blue color scheme was chosen to be visually similar to the Ali2D HTML viewer but can be changed by passing different colors as 0-1 RGB triplets as 'r,g,b' to `--color_helix` and `--color_sheet`, like `--color_helix '0,1,0` would make the alpha helix confidences scale from white to green.
