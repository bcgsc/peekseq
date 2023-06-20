# peekseq
## Protein Estimation systEm in [DNA/RNA] SEQuences, using K-mers
### De novo protein-coding potential predictor

### Test data
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

Example commands for C. maximus mitogenome and SARS-CoV-2:
/usr/bin/time ../peekseq.pl -f CEMA.fa.gz -k 150 -s 200 -c 2 -v 1&
/usr/bin/time ../peekseq.pl -f SARS.fa.gz -k 150 -s 270 -c 11 -v 1&

</pre>
