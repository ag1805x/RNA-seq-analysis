#!/usr/bin/perl
use Cwd;
#Perl script to covert multiple SAM files to BAM files using SAMTOOLS

#Array holding all the SRA ids
@srr=(SRR1234567, SRR1234568, SRR1234569);


$cwd = getcwd;
print "\n\n Current directory is $cwd\n\n";


$samtools = 'samtools';


for ($i=0; $i< @srr; $i++)
	{
	 print "\nRUNNING SAM to BAM CONVERSION FOR $srr[$i]\n";
	 
	 system ("$samtools sort -o $srr[$i].bam $srr[$i].sam "); 
   
	};

print "\n\nEND\n\n";
