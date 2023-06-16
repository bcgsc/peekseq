![Logo](https://github.com/bcgsc/cpkseq/blob/main/cpkseq-logo.png)

# cpkseq
## De novo [protein] coding potential calculator using a k-mer approach
### 2023


## Contents

1. [Description](#description)
2. [Implementation and requirements](#implementation)
3. [Install](#install)
4. [Dependencies](#dep)
5. [Documentation](#docs)
6. [Citing cpkseq](#cite)
7. [Credits](#credits)
8. [Running cpkseq](#run)
9. [Test data](#data)
10. [Output](#output)
11. [Algorithm](#algorithm)
12. [Quick reference](#quickref)
13. [Generating plots](#bplot)
14. [License](#license)


## Description <a name=description></a>

cpkseq systematically processes the k-mers of a supplied DNA sequence and computes regions with coding potential, using target codon translation tables

cpkseq is used to quickly survey regions with [protein] coding potential on all 6 frames


## Implementation and requirements <a name=implementation></a>

cpkseq is developed in PERL and runs on any system where PERL is installed.


## Install <a name=install></a>

Clone and enter the unikseq directory.
<pre>
git clone https://github.com/bcgsc/cpkseq
cd cpkseq
</pre>



## Dependencies <a name=dep></a>

If PERL is installed on your system, you're good to go (no additional libraries needed, nor dependencies).


## Documentation <a name=docs></a>

Refer to the README.md file on how to install and run cpkseq. Read the 


## Citing cpkseq <a name=cite></a>

Thank you for your [![Stars](https://img.shields.io/github/stars/bcgsc/cpkseq.svg)](https://github.com/bcgsc/cpkseq/stargazers) and for using, developing and promoting this free software!

If you use unikseq in your research, please cite: 

TBD

## Credits <a name=credits></a>
<pre>
cpkseq concept, algorithm design, implementation: Rene L Warren
</pre>


## Running cpkseq <a name=run></a>

<pre>
Usage: ./cpkseq.pl VERSION
 -f FASTA (DNA/RNA) file (required)
 -k length (option, default: -k 90)
 -c genetic code translation table id (option, default: -c 1 [standard])
 -s min. reference FASTA region [size] (bp) to output (option, default: -s 270 bp)
 -v output tsv file (option, -v 1==yes -v 0==no [default])
</pre>

Notes:
<pre>


 -f FASTA (DNA/RNA) file (required)

 -k length (option, default: -k 90)

 -c genetic code translation table id (option, default: -c 1 [standard])

 -s min. reference FASTA region [size] (bp) to output (option, default: -s 270 bp)

 -v output tsv file (option, -v 1==yes -v 0==no [default])


</pre>

### Test data <a name=data></a>
---------

TBD

## Output  <a name=output></a>

1) TSV file 


2) FASTA file (*DNA.fa)

   A multi-FASTA 

3) FASTA file (*PROTEIN.fa)


4) LOG file (.log)

   Captures the verbose STDOUT in a log file, showing the progress of the program.    


## Algorithm design and implementation <a name=algorithm></a>

TBD
 
## Quick reference <a name=quickref></a>

TBD

## Generating plots <a name=bplot></a>

TBD

## License <a name=license></a>

Unikseq Copyright (c) 2020-2023 British Columbia Cancer Agency Branch.  All rights reserved.

Unikseq is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

For commercial licensing options, please contact
Patrick Rebstein <prebstein@bccancer.bc.ca>
