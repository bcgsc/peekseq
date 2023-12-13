[![Release](https://img.shields.io/github/release/bcgsc/peekseq.svg)](https://github.com/bcgsc/peekseq/releases)
[![Conda](https://img.shields.io/conda/dn/bioconda/peekseq?label=Conda)](https://anaconda.org/bioconda/peekseq)
[![Issues](https://img.shields.io/github/issues/bcgsc/peekseq.svg)](https://github.com/bcgsc/peekseq/issues)

![Logo](https://github.com/bcgsc/peekseq/blob/main/peekseq-logo.png)

# peekseq
## Protein Estimation systEm in [DNA/RNA] SEQuences, using K-mers
### De novo protein-coding potential predictor


## Contents

1. [Description](#description)
2. [Implementation and requirements](#implementation)
3. [Install](#install)
4. [Dependencies](#dep)
5. [Documentation](#docs)
6. [Citing peekseq](#cite)
7. [Credits](#credits)
8. [Running peekseq](#run)
9. [Test data](#data)
10. [Output](#output)
11. [Algorithm](#algorithm)
12. [Quick reference](#quickref)
13. [Generating plots](#bplot)
14. [License](#license)


## Description <a name=description></a>

peekseq systematically processes the k-mers of a supplied DNA/RNA sequence and computes regions with coding potential, using target codon translation tables

peekseq is used to quickly survey regions with [protein] coding potential on all 6 frames, simultaneously

peekseq is best suited for evaluating the protein-coding potential of simple DNA sequences, including bacterial gene structures (i.e., it will not predict intron/exon boundaries, etc.)


## Implementation and requirements <a name=implementation></a>

peekseq is developed in PERL and runs on any system where PERL is installed.


## Install <a name=install></a>

Clone and enter the peekseq directory.
<pre>
git clone https://github.com/bcgsc/peekseq
cd peekseq
</pre>

or using conda:
<pre>
conda install -c bioconda peekseq
</pre>


## Dependencies <a name=dep></a>

If PERL is installed on your system, you're good to go (no additional libraries needed, nor dependencies).


## Documentation <a name=docs></a>

Refer to the README.md file on how to install and run peekseq. 


## Citing peekseq <a name=cite></a>

Thank you for your [![Stars](https://img.shields.io/github/stars/bcgsc/peekseq.svg)](https://github.com/bcgsc/peekseq/stargazers) and for using, developing and promoting this free software!

If you use peekseq in your research, please cite: 

TBD


## Credits <a name=credits></a>
<pre>
peekseq concept, algorithm design, implementation: Rene L Warren
</pre>


## Running peekseq <a name=run></a>

<pre>
Usage: ./peekseq.pl
 -f FASTA (DNA/RNA) file (required)
 -k length (option, default: -k 90)
 -c genetic code translation table id (option, default: -c 1 [standard])
 -s min. reference FASTA region [size] (bp) to output (option, default: -s 270 bp)
 -v output tsv file (option, -v 1==yes -v 0==no [default])
</pre>

Notes:
<pre>

 -f FASTA (DNA/RNA) file (required)
 Single or Multi-FASTA file to predict open-reading frames from 

 -k length (option, default: -k 90)
 The (nucleotide) k-mer length (nt)

 -c genetic code translation table id (option, default: -c 1 [standard])
 The table to use for translation. Please refer to:
 https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c#top
 <pre> 
    1. The Standard Code
    2. The Vertebrate Mitochondrial Code
    3. The Yeast Mitochondrial Code
    4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
    5. The Invertebrate Mitochondrial Code
    6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
    9. The Echinoderm and Flatworm Mitochondrial Code
    10. The Euplotid Nuclear Code
    11. The Bacterial, Archaeal and Plant Plastid Code
    12. The Alternative Yeast Nuclear Code
    13. The Ascidian Mitochondrial Code
    14. The Alternative Flatworm Mitochondrial Code
    16. Chlorophycean Mitochondrial Code
    21. Trematode Mitochondrial Code
    22. Scenedesmus obliquus Mitochondrial Code
    23. Thraustochytrium Mitochondrial Code
    24. Rhabdopleuridae Mitochondrial Code
    25. Candidate Division SR1 and Gracilibacteria Code
    26. Pachysolen tannophilus Nuclear Code
    27. Karyorelict Nuclear Code
    28. Condylostoma Nuclear Code
    29. Mesodinium Nuclear Code
    30. Peritrich Nuclear Code
    31. Blastocrithidia Nuclear Code
    33. Cephalodiscidae Mitochondrial UAA-Tyr Code
  </pre>

 -s min. reference FASTA region [size] (bp) to output (option, default: -s 270 bp)
  minimum sequence length to output  

 -v output tsv file (option, -v 1==yes -v 0==no [default])
  boolean, 1/0=yes/no, outputs a tsv file with information for plotting the 6-frame codin potential (see below section Output (1) TSV file) 

</pre>

### Test data <a name=data></a>
---------

<pre>
1. Go to ./testdata
(cd testdata)

2. Run peekseq on the provided test data

Example commands for C. maximus mitogenome and SARS-CoV-2 genome:
/usr/bin/time ../peekseq.pl -f CEMA.fa.gz -k 150 -s 200 -c 2 -v 1
/usr/bin/time ../peekseq.pl -f SARS.fa.gz -k 150 -s 270 -c 11 -v 1&

These commands will generate four output files for each genome:
e.g.
peekseq_v0.0.1-f_CEMA.fa.gz-k150-frameKmers.tsv
peekseq_v0.0.1-f_CEMA.fa.gz-k150-s200-c2-codingDNA.fa
peekseq_v0.0.1-f_CEMA.fa.gz-k150-s200-c2-codingPROTEIN.fa
peekseq_v0.0.1-f_CEMA.fa.gz-k150-s200-c2.log

They should take no more than 1 and 2 seconds, respectively, to run on a MacBook Pro (Catalina10.15.7, 2.6 GHz 6-Core Intel Core i7).

If the run is successful, the peekseq_v0.0.1-f_CEMA.fa.gz-k150-s200-c2-codingDNA.fa output should contain 15 sequences and the peekseq_v0.0.1-f_SARS.fa.gz-k150-s270-c2-codingDNA.fa, 30 sequences.


</pre>


## Output  <a name=output></a>

1) TSV file 

   A tab-separated file containing the FASTA header, position, frame and whether the corresponding k-mer encodes a peptide stretch (1 or -1 for plus/minus strands, respectively).
   <pre>
   header	position	frame	coding
   ...
   NZ_CP028101.1   3030    1       0
   NZ_CP028101.1   3030    -1      -1
   NZ_CP028101.1   3031    2       0
   NZ_CP028101.1   3031    -2      0
   NZ_CP028101.1   3032    3       1
   NZ_CP028101.1   3032    -3      0
   NZ_CP028101.1   3033    1       0
   NZ_CP028101.1   3033    -1      -1
   NZ_CP028101.1   3034    2       0
   NZ_CP028101.1   3034    -2      0
   NZ_CP028101.1   3035    3       1
   ...
   </pre>

2) FASTA file (*-codingDNA.fa)

   A multi-FASTA with sequence region for any given frames with a corresponding translation in -codingPROTEIN.fa, excised from the original input FASTA. DNA/RNA outside of the initiation codon are soft-masked (i.e., lower case bases).


3) FASTA file (*-codingPROTEIN.fa)

   A multi-FASTA with protein translation corresponding the regions in -codingDNA.fa. Residues outside of the initiation amino acid are soft-masked (i.e., lower case residues).


4) LOG file (.log)

   Captures the verbose STDOUT in a log file, showing the progress of the program.    


## Algorithm design and implementation <a name=algorithm></a>

Peekseq uses words of length k (k-mers) to quickly assess the coding potential simultaneously over the 6 frames, by sliding over each k frame of the input FASTA.
 

## Quick reference <a name=quickref></a>

TBD

## Generating plots <a name=bplot></a>

### Basking shark mitogenome
![peekseqPlot](https://github.com/bcgsc/peekseq/blob/main/CEMAcpPlot.png)
This example predicted coding potential on C. maximus (basking shark) mitogenome with the following command:
<pre>
/usr/bin/time ./peekseq.pl -f CEMA.fa.gz -k 150 -s 200 -c 2 -v 1
</pre>

Plotted using the coding potential tsv file (-v 1 : -frameKmers.tsv) generated by peekseq.
R code to generate similar plots
<pre>
#coding potential plots

library(ggplot2)
library(ggdark)

df <- read.table('peekseq_v0.0.1-f_CEMA.fa.gz-k150-frameKmers.tsv', sep="\t", header = TRUE)
my_x_title <- expression(paste(italic("Cetorhinus maximus"), " [basking shark] mitogenome (Genbank accession KF597303.1)"))

p<-ggplot(df, aes(y=coding, x=position, fill = factor(frame))) +
  geom_col() + 
  scale_fill_manual(values = c("#b35806","#f1a340","#fee0b6","#d8daeb","#998ec3","#542788")) +
  ylim(-1,1) + labs(x=my_x_title) + ylab("Coding potential") + dark_theme_gray(base_size = 14)  # Default
p + theme(legend.position = "bottom") + guides(fill = guide_legend(title = "Frame", nrow = 1))
</pre>

### SARS-CoV-2 genome
![peekseqPlot](https://github.com/bcgsc/peekseq/blob/main/SARS2cpPlot.png)
This example predicted coding potential regions within the SARS-CoV-2 genome (k150) at the pandemic onset (Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1) with the following command:
<pre>
/usr/bin/time ./peekseq.pl -f SARS.fa.gz -k 150 -s 270 -c 11 -v 1&
</pre>

## License <a name=license></a>

peekseq Copyright (c) 2023-present British Columbia Cancer Agency Branch.  All rights reserved.

peekseq is released under the GNU General Public License v3

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
