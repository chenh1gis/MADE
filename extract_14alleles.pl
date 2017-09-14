#!/usr/bin/perl

##########################################################################
## Author: Hui Chen
## Created Time: 2017-8-17 09:40:44
## File Name: extract_allele.pl
##########################################################################
#use strict;
#use warnings;

$FILE=$ARGV[0];
$OUTFILE=$ARGV[1];

@CODON=(138,145,156,158,159,160,183,186,190,193,194,219,226,246);

%Genetic_code=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'-','TAG'=>'-','TGC'=>'C','TGT'=>'C','TGA'=>'-','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');

open (FILE,"$FILE");
$seq=<FILE>;
chomp $seq;
if ($seq=~/\>/)
{
	$seq=<FILE>;
	chomp $seq;	
}

$base=substr($seq,0,3);
if ($Genetic_code{$base} != "Q")
{
	print "The sequence of H3 HA1 is not started with \"QKLP\", please double check the sequence!\n";
	exit(0);	
}
$len=length($seq);
print "Length of sequence is $len\n";
if ($len<3*$CODON[$#CODON])
{
        print "The sequence of H3 HA1 is not complete and information of some codons are missing!\n";
        exit(0);
}

else
{
	open (DEMO,">$OUTFILE");
	print DEMO "CODON\tAminoAcid\n";
	foreach $codon(@CODON)	
	{
		$base=substr($seq,3*($codon-1),3);
		print DEMO "$codon\t$Genetic_code{$base}\n";
	}
	close DEMO;
}
