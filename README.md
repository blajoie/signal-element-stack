<img height=40 src='http://my5C.umassmed.edu/images/3DG.png' title='3D-Genome' />
&nbsp;&nbsp;
<img height=30 src='http://my5C.umassmed.edu/images/dekkerlabbioinformatics.gif' />
&nbsp;&nbsp;
<img height=30 src='http://my5C.umassmed.edu/images/umasslogo.gif' />

# single-element-stack

stack/aggregrate signal around a set of elements/features

```
single-element-stack/
	single-element-stack.py - stack/aggregrate signal around a set of elements/features
```

## Installation

The main requirements are numpy

You can install the dependencies with

for req in $(cat requirements.txt); do pip install $req; done

## Full Documentation

## Communication

- [Bryan Lajoie](https://github.com/blajoie)
- Twitter: [@my5C](https://twitter.com/my5C)

## What does it do?

single-element-stack can stack/aggregrate signal around a set of elements/features

## Usage

```

$ python scripts/signal-element-stack.py --help
usage: signal-element-stack.py [-h] -i SIGNALTRACK -e ELEMENTTRACK -g GENOME
                               [-w PILEUPWIDTH] [-b NUM_BINS]
                               [-ya YAXISRANGE [YAXISRANGE ...]]
                               [-sc SIGNALCOLUMN] [-sm] [-mp] [-v]

signal pile up plot for bedGraph/wiggle tracks

optional arguments:
  -h, --help            show this help message and exit
  -i SIGNALTRACK, --input SIGNALTRACK
                        input bed/wiggle track to be used for the pile up
                        (default: None)
  -e ELEMENTTRACK, --element ELEMENTTRACK
                        element bed track for pile up (default: None)
  -g GENOME, --genome GENOME
                        genome chromosome length file (default: None)
  -w PILEUPWIDTH, --width PILEUPWIDTH
                        width of pile up plot (in BP) (default: 500)
  -b NUM_BINS, --num_bins NUM_BINS
                        number of num_bins for pile up plot width (default:
                        500)
  -ya YAXISRANGE [YAXISRANGE ...], --yaxisrange YAXISRANGE [YAXISRANGE ...]
                        y-axis for aggregrate plot (default: [0, 0])
  -sc SIGNALCOLUMN, --signalcolumn SIGNALCOLUMN
                        signal column number in bed/bedgraph/wig/tsv file
                        (default: 3)
  -sm, --sortMatrix     sort pile up matrix by total row signal, descending
                        (default: False)
  -mp, --midpointMode   use midpoint of element, otherwise will assume TSS
                        (default: False)
  -v, --verbose         Increase verbosity (specify multiple times for more)
                        (default: None)
                        
## Usage Examples

```

## Change Log

## Bugs and Feedback

For bugs, questions and discussions please use the [Github Issues](https://github.com/blajoie/single-element-stack/issues).

## LICENSE

Licensed under the Apache License, Version 2.0 (the 'License');
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

<http://www.apache.org/licenses/LICENSE-2.0>

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an 'AS IS' BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
