#!/usr/bin/env perl

#AUTHOR
#   Rene Warren
#   rwarren at bcgsc.ca

#NAME
# peskseq : Protein evaluation System with a K-mer approach


#SYNOPSIS
# Identifies DNA regions with coding potential de novo, using a k-mer approach

#DOCUMENTATION
#   README.md distributed with this software
#   We hope this code is useful to you -- Please send comments & suggestions to rwarren * bcgsc.ca

#LICENSE
#peskseq Copyright (c) 2023 British Columbia Cancer Agency Branch.  All rights reserved.
#peskseq is released under the GNU General Public License v3

use strict;
use Getopt::Std;

use vars qw($opt_k $opt_f $opt_s $opt_v $opt_c);
getopts('k:f:s:v:c:');

my $version = "v0.0.1";
my ($k,$tsvflag,$code) = (90,0,1);
my $regsz = 3 * $k;
system("clear");

if(! $opt_f){
   print "Usage: $0 $version\n";
   print " -f FASTA (DNA/RNA) file (required)\n";
   print " -k length (option, default: -k $k)\n";
   print " -c genetic code translation table id (option, default: -c $code [standard])\n";
   print " -s min. reference FASTA region [size] (bp) to output (option, default: -s $regsz bp)\n";
   die " -v output tsv file (option, -v 1==yes -v 0==no [default])\n\n";
}

### Fetch options

my $f = $opt_f; #reference
$k = $opt_k if($opt_k);
$regsz = $opt_s if($opt_s);
$tsvflag = $opt_v if($opt_v);
$code = $opt_c if($opt_c);


###Prepare output
#-----
my $fn = "peskseq_" . $version . "-f_" . $f . "-k" . $k;
my $tsv= $fn . "-frameKmers.tsv";

$fn .= "-s" . $regsz . "-c" . $code;

my $dna=$fn . "-codingDNA.fa";
my $pep=$fn . "-codingPROTEIN.fa";

my $log=$fn . ".log";

open(LOG,">$log") || die "Can't write to $log -- fatal.\n";

my $message = "\nRunning: $0 $version\n\t-k $k\n\t-f $f\n\t-s $regsz\n\t-c $code\n\t-v $tsvflag\n";

print $message;
print LOG $message;


###Checking files and options
#-----
if(! -e $f){
   die "Invalid file: $f -- fatal\n";
}

if($k % 3 != 0){
   die "-k $k must be a multiple of 3 -- fatal.\n";
}

&evalCP($f,$k,$dna,$pep,$regsz,$tsv,$code,$tsvflag);

$message = "done.\n";
$message .= "-" x 30, "\n";
$message .= "\nOutput in-frame regions >= $regsz bp in:\n$dna\n$pep\n";
$message .= "\nOutput $k-mers frames:\n$tsv\ndone.\n" if($tsvflag);
print $message;
print LOG $message;

close LOG;

exit;

#--------------------------------
sub evalCP{
  
   my ($f,$k,$dna,$pep,$regsz,$tsv,$code,$tsvflag) = @_;
   
   my ($head,$prevhead,$seq,$preventry) = ("","","","");

   open(DNA,">$dna") || die "Can't write $dna -- fatal.\n";
   open(PEP,">$pep") || die "Can't write $pep -- fatal.\n";

   if($tsvflag){
      open(TSV,">$tsv") || die "Can't write $tsv -- fatal.\n";
      print TSV "header\tposition\tframe\tcoding\n";
   }

   my ($codon2aa,$startcodon) = &initializeAA($code);

   open(IN,$f) || die "Can't read $f -- fatal.\n";
   while(<IN>){
      s/\r\n/\n/g;### DOS to UNIX
      chomp;

      if(/^\>(\S+)/){
         my $entry = $1;
         $head = $entry;#$_;
         if($prevhead ne $head && $prevhead ne "" && $seq ne ""){
            &printOutput($codon2aa,$startcodon,$prevhead,$seq,$preventry,$k,$regsz,$tsvflag);
         }
         $seq = "";
         $prevhead = $head;
         $preventry = $entry;
      }else{
         my $seqstretch = $1 if(/^(\S+)/); ###this prevents DOS new lines from messing up the TSV output
         $seq .= uc($seqstretch);
      }
   }
  
   &printOutput($codon2aa,$startcodon,$prevhead,$seq,$preventry,$k,$regsz,$tsvflag);

   close DNA;
   close PEP;
   close TSV if($tsvflag);   
}

#--------------------------------
sub printOutput{

   my ($codon2aa,$startcodon,$head,$seq,$preventry,$k,$regsz,$tsvflag) = @_;

   my ($initial,$unique,$notunique,$sum,$sumout,$buffer) = (-1,0,0,0,0,0);
   my $s;
   my $lastframe = 0;
   my $lastframerev = 0;
   my $lastcodon = "";
   for(my $pos=0;$pos<=(length($seq)-$k);$pos++){
      my $kmer = substr($seq,$pos,$k);
      my $codon = substr($kmer,($k-3),3);
      $lastcodon = $codon;
      $kmer =~ tr/U/T/; ### handles RNA U>>>T
      my $rckmer = &reverseComplement($kmer);

      my $frame = ($pos % 3) + 1;
      my $framerev = -1 * $frame;
      $lastframe = $frame;
      $lastframerev = $framerev;

      my $aamer = &dna2protein($kmer,$codon2aa);
      my $rcaamer = &dna2protein($rckmer,$codon2aa);

      my $bit=1;
      if($aamer=~/\*/ || $aamer=~/X/){# stop or end
         $bit=0;
         if(length($s->{$frame}{'seq'}) >= $regsz){
            $s->{$frame}{'seq'} .= $codon;
            $s->{$frame}{'end'} = $pos + $k;
            ###
            my $lcdna = lc($s->{$frame}{'seq'});
            my @ntarr = ( $lcdna =~ m/.../g );
            my $flag=0;
            my $newnt="";
            my $newaa="";
            foreach my $cd (@ntarr){
               my $aa = lc(&dna2protein($cd,$codon2aa));
               if(defined $startcodon->{$cd}){$flag=1;}
               if($flag){
                  $newnt .= uc($cd);
                  $newaa .= uc($aa);
               }else{
                  $newnt .= $cd;
                  $newaa .= $aa;
               }
            }
            ###
            my $ntlen = length($newnt);
            print DNA ">$head $s->{$frame}{'start'}-$s->{$frame}{'end'} length:$ntlen frame:$frame\n$newnt\n";
            my $aalen = length($newaa);
            print PEP ">$head $s->{$frame}{'start'}-$s->{$frame}{'end'} length:$aalen frame:$frame\n$newaa\n";
         }
         $s->{$frame}{'start'} = 0;
         $s->{$frame}{'seq'} = ""; 
      }else{
         if($s->{$frame}{'seq'} eq ""){
            $s->{$frame}{'seq'} = $kmer;
            $s->{$frame}{'start'} = $pos;
         }else{
            $s->{$frame}{'seq'} .= $codon;
            $s->{$frame}{'end'} = $pos + $k;
         }##
      }
      my $rcbit=-1;
      if($rcaamer=~/\*/ || $rcaamer=~/X/){###stop or end
         $rcbit=0;
         if(length($s->{$framerev}{'seq'}) >= $regsz){
            my $rc = reverseComplement($s->{$framerev}{'seq'});
            my $lcdna = lc($rc);
            my @ntarr = ( $lcdna =~ m/.../g );
            my $flag=0;
            my $newnt="";
            my $newaa="";
            foreach my $cd (@ntarr){
               my $aa = lc(&dna2protein($cd,$codon2aa));
               if(defined $startcodon->{$cd}){$flag=1;}
               if($flag){
                  $newnt .= uc($cd);
                  $newaa .= uc($aa);
               }else{
                  $newnt .= $cd;
                  $newaa .= $aa;
               }
            }
            my $ntlen = length($newnt);
            print DNA ">$head $s->{$framerev}{'end'}-$s->{$framerev}{'start'} length:$ntlen frame:$framerev\n$newnt\n";
            my $aalen = length($newaa);
            print PEP ">$head $s->{$framerev}{'end'}-$s->{$framerev}{'start'} length:$aalen frame:$framerev\n$newaa\n";
         }
         $s->{$framerev}{'start'} = $pos;
         $s->{$framerev}{'seq'} = $kmer;
      }else{
         if($s->{$framerev}{'seq'} eq ""){
            $s->{$framerev}{'seq'} = $kmer;
            $s->{$framerev}{'start'} = $pos;
         }else{
            $s->{$framerev}{'seq'} .= $codon;
            $s->{$framerev}{'end'} = $pos + $k;
         }
      }
      print TSV "$head\t$pos\t$frame\t$bit\n$head\t$pos\t$framerev\t$rcbit\n" if($tsvflag);
   }###end of sequence
   ###print what is in the buffer (like a stop, for all 3 frames)
   my @frames = (1,2,3);
   foreach my $fr(@frames){
      if(length($s->{$fr}{'seq'}) >= $regsz){
         #$s->{$lastframe}{'seq'} .= $lastcodon;
         #$s->{$lastframe}{'end'} = length($seq);
         ###
         my $lcdna = lc($s->{$fr}{'seq'});
         my @ntarr = ( $lcdna =~ m/.../g );
         my $flag=0;
         my $newnt="";
         my $newaa="";
         foreach my $cd (@ntarr){
            my $aa = lc(&dna2protein($cd,$codon2aa));
            if(defined $startcodon->{$cd}){$flag=1;}
            if($flag){
               $newnt .= uc($cd);
               $newaa .= uc($aa);
            }else{
               $newnt .= $cd;
               $newaa .= $aa;
            }
         }
         ###
         my $ntlen = length($newnt);
         print DNA ">$head $s->{$fr}{'start'}-$s->{$fr}{'end'} length:$ntlen frame:$fr\n$newnt\n";
         my $aalen = length($newaa);
         print PEP ">$head $s->{$fr}{'start'}-$s->{$fr}{'end'} length:$aalen frame:$fr\n$newaa\n";
      }#length >= user-defined threshold
   }#3 frames
   ###handle reverse
   my @framesrec=(-1,-2,-3);
   foreach my $fr(@framesrec){
      if(length($s->{$fr}{'seq'}) >= $regsz){
         my $rc = reverseComplement($s->{$fr}{'seq'});
         my $lcdna = lc($rc);
         my @ntarr = ( $lcdna =~ m/.../g );
         my $flag=0;
         my $newnt="";
         my $newaa="";
         foreach my $cd (@ntarr){
            my $aa = lc(&dna2protein($cd,$codon2aa));
            if(defined $startcodon->{$cd}){$flag=1;}
            if($flag){
               $newnt .= uc($cd);
               $newaa .= uc($aa);
            }else{
               $newnt .= $cd;
               $newaa .= $aa;
            }
         }
         my $ntlen = length($newnt);
         print DNA ">$head $s->{$fr}{'end'}-$s->{$fr}{'start'} length:$ntlen frame:$fr\n$newnt\n";
         my $aalen = length($newaa);
         print PEP ">$head $s->{$fr}{'end'}-$s->{$fr}{'start'} length:$aalen frame:$fr\n$newaa\n";
      }### length>= user-defined threshold
   }### all 3 frames rev
}

#--------------------------------
sub initializeAA{

    my $code=shift;

    my $codon2aa;
    my $startcodon;

    $codon2aa->{"ttt"}= "F"; $codon2aa->{"ttc"}= "F"; $codon2aa->{"tta"}="L";
    $codon2aa->{"ttg"}= "L"; $codon2aa->{"tct"}= "S"; $codon2aa->{"tcc"}= "S";
    $codon2aa->{"tca"}= "S"; $codon2aa->{"tcg"}= "S"; $codon2aa->{"tat"}= "Y";
    $codon2aa->{"tac"}= "Y"; $codon2aa->{"tgt"}= "C"; $codon2aa->{"tgc"}= "C";
    $codon2aa->{"ctt"}= "L"; $codon2aa->{"ctc"}= "L"; $codon2aa->{"cta"}= "L";
    $codon2aa->{"ctg"}= "L"; $codon2aa->{"cct"}= "P"; $codon2aa->{"ccc"}= "P";
    $codon2aa->{"cca"}= "P"; $codon2aa->{"ccg"}= "P"; $codon2aa->{"cat"}= "H";
    $codon2aa->{"cac"}= "H"; $codon2aa->{"caa"}= "Q"; $codon2aa->{"cag"}= "Q";
    $codon2aa->{"cgt"}= "R"; $codon2aa->{"cgc"}= "R"; $codon2aa->{"cga"}= "R";
    $codon2aa->{"cgg"}= "R"; $codon2aa->{"att"}= "I"; $codon2aa->{"atc"}= "I";
    $codon2aa->{"ata"}= "I"; $codon2aa->{"atg"}= "M"; $codon2aa->{"act"}= "T";
    $codon2aa->{"acc"}= "T"; $codon2aa->{"aca"}= "T"; $codon2aa->{"acg"}= "T";
    $codon2aa->{"aat"}= "N"; $codon2aa->{"aac"}= "N"; $codon2aa->{"aaa"}= "K";
    $codon2aa->{"aag"}= "K"; $codon2aa->{"agt"}= "S"; $codon2aa->{"agc"}= "S";
    $codon2aa->{"aga"}= "R"; $codon2aa->{"agg"}= "R"; $codon2aa->{"gtt"}= "V";
    $codon2aa->{"gtc"}= "V"; $codon2aa->{"gta"}= "V"; $codon2aa->{"gtg"}= "V";
    $codon2aa->{"gct"}= "A"; $codon2aa->{"gcc"}= "A"; $codon2aa->{"gca"}= "A";
    $codon2aa->{"gcg"}= "A"; $codon2aa->{"gat"}= "D"; $codon2aa->{"gac"}= "D";
    $codon2aa->{"gaa"}= "E"; $codon2aa->{"gag"}= "E"; $codon2aa->{"ggt"}= "G";
    $codon2aa->{"ggc"}= "G"; $codon2aa->{"gga"}= "G"; $codon2aa->{"ggg"}= "G";
    $codon2aa->{"tgg"}= "W";
    $codon2aa->{"tag"}= "*";
    $codon2aa->{"tga"}= "*";
    $codon2aa->{"taa"}= "*";

    $startcodon->{"atg"} = "M";
    #The standard code currently allows initiation from UUG and CUG in addition to AUG.
    $startcodon->{"ttg"} = "L";
    $startcodon->{"ctg"} = "L";
    $startcodon->{"gtg"} = "V";

    #There are Differences from the Standard Code
    #https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c#top
    #1. The Standard Code
    #2. The Vertebrate Mitochondrial Code
    #3. The Yeast Mitochondrial Code
    #4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
    #5. The Invertebrate Mitochondrial Code
    #6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
    #9. The Echinoderm and Flatworm Mitochondrial Code
    #10. The Euplotid Nuclear Code
    #11. The Bacterial, Archaeal and Plant Plastid Code
    #12. The Alternative Yeast Nuclear Code
    #13. The Ascidian Mitochondrial Code
    #14. The Alternative Flatworm Mitochondrial Code
    #16. Chlorophycean Mitochondrial Code
    #21. Trematode Mitochondrial Code
    #22. Scenedesmus obliquus Mitochondrial Code
    #23. Thraustochytrium Mitochondrial Code
    #24. Rhabdopleuridae Mitochondrial Code
    #25. Candidate Division SR1 and Gracilibacteria Code
    #26. Pachysolen tannophilus Nuclear Code
    #27. Karyorelict Nuclear Code
    #28. Condylostoma Nuclear Code
    #29. Mesodinium Nuclear Code
    #30. Peritrich Nuclear Code
    #31. Blastocrithidia Nuclear Code
    #33. Cephalodiscidae Mitochondrial UAA-Tyr Code

    if($code==2){
       $codon2aa->{"aga"}="*";
       $codon2aa->{"agg"}="*";
       $codon2aa->{"ata"}="M";
       $codon2aa->{"tga"}="W";
       $startcodon->{"ata"}="M";
       $startcodon->{"att"}="I";
       $startcodon->{"atc"}="I";
    }elsif($code==3){
       $codon2aa->{"ata"}="M"; 
       $codon2aa->{"ctt"}="T";
       $codon2aa->{"ctc"}="T";
       $codon2aa->{"cta"}="T";
       $codon2aa->{"ctg"}="T";
       $codon2aa->{"tga"}="W";
       $startcodon->{"gtg"}="V";
    }elsif($code==4){
       $codon2aa->{"tga"}="W";
    }elsif($code==5){
       $codon2aa->{"tga"}="W";
       $codon2aa->{"ata"}="M";
       $codon2aa->{"agg"}="S";
       $codon2aa->{"aga"}="S";
       $startcodon->{"ata"}="M";
       $startcodon->{"att"}="I";
    }elsif($code==6){
       $codon2aa->{"taa"}="Q";
       $codon2aa->{"tag"}="Q";
    }elsif($code==9){
       $codon2aa->{"aaa"}="N";
       $codon2aa->{"aga"}="S";
       $codon2aa->{"agg"}="S";
       $codon2aa->{"tga"}="W";
       $startcodon->{"gtg"}="V";
       $startcodon->{"ttg"}="L";
    }elsif($code==10){
       $codon2aa->{"tga"}="C";
    }elsif($code==12){
       $codon2aa->{"ctg"}="S";
    }elsif($code==13){
       $codon2aa->{"aga"}="G";
       $codon2aa->{"agg"}="G";
       $codon2aa->{"ata"}="M";
       $codon2aa->{"tga"}="W";
       $startcodon->{"gtg"}="V";
       $startcodon->{"ttg"}="L";
       $startcodon->{"ata"}="M";
    }elsif($code==14){
       $codon2aa->{"aaa"}="N";
       $codon2aa->{"aga"}="S";
       $codon2aa->{"agg"}="S";
       $codon2aa->{"taa"}="Y";
       $codon2aa->{"tga"}="W";
    }elsif($code==16){
       $codon2aa->{"tag"}="L";
    }elsif($code==21){
       $codon2aa->{"tga"}="W";
       $codon2aa->{"ata"}="M";
       $codon2aa->{"aga"}="S";
       $codon2aa->{"agg"}="S";
       $codon2aa->{"aaa"}="N";
    }elsif($code==22){
       $codon2aa->{"tca"}="*";
       $codon2aa->{"tag"}="L";
    }elsif($code==24){
       $codon2aa->{"aga"}="S";
       $codon2aa->{"agg"}="K";
       $codon2aa->{"tga"}="W";
    }elsif($code==25){
       $codon2aa->{"tga"}="G";
       $startcodon->{"gtg"}="V";
       $startcodon->{"ttg"}="L";
    }elsif($code==26){
       $codon2aa->{"ctg"}="A";
       $startcodon->{"gtg"}="V";
       $startcodon->{"ttg"}="L";
    }elsif($code==27){
       $codon2aa->{"tag"}="Q";
       $codon2aa->{"taa"}="Q";
       #or $codon2aa->{"tga"}="W";
    }elsif($code==28){
       #or$codon2aa->{"tag"}="Q";
       #or $codon2aa->{"taa"}="Q";
       #or $codon2aa->{"tga"}="W";
    }elsif($code==29){
       $codon2aa->{"tag"}="Y";
       $codon2aa->{"taa"}="Y";
    }elsif($code==30){
       $codon2aa->{"tag"}="E";
       $codon2aa->{"taa"}="E";
    }elsif($code==31){
       $codon2aa->{"tga"}="W";
       #or $codon2aa->{"tag"}="E";
       #or $codon2aa->{"taa"}="E";
    }elsif($code==33){
       $codon2aa->{"aga"}="S";
       $codon2aa->{"taa"}="Y";
       $codon2aa->{"tga"}="W";
       $codon2aa->{"agg"}="K";
    }

    return $codon2aa,$startcodon;
}

#--------------------------------
sub dna2protein{

    my ($seq,$codon2aa) = @_;

    my $seq = lc($seq);

    my @codonlist = ( $seq =~ m/.../g );

    my $aa_seq = "";

    foreach my $codon (@codonlist){
       #print ">>>$codon<<<\n";
       if (defined $codon2aa->{$codon}){
          $aa_seq.=$codon2aa->{$codon};
       }else{
          $aa_seq.="X";
       }
    }
    return $aa_seq;
}

#--------------------------------
sub reverseComplement{
   $_ = shift;
   $_ = uc();
   tr/ATGCYRKMBDHV/TACGRYMKVHDB/;
   return (reverse());
}
