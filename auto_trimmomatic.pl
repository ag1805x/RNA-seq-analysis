#################################################
# Perl script to trim multiple RNAseq raw reads #
# using trimmomatic                             #
#         -Arindam Ghosh (15 May 2018)          #
#################################################   


#!/usr/bin/perl
use Cwd;


@srr=(SRR1234567, SRR1234568, SRR1234569);


$cwd = getcwd;
print "\n\n Current directory is $cwd\n\n";


$java = 'java';


for ($i=0; $i< @srr; $i++)
	{
	 print "\nRUNNING RAW READ TRIMMING WITH TRIMMOMATIC FOR $srr[$i]\n";
	 mkdir $srr[$i];
         chdir("$cwd/$srr[$i]");


# Command for PAIRED END READS
system ("$java -jar /path/to/trimmomatic-0.36.jar PE -trimlog $srr[$i]_trim_log /path/to/read1/$srr[$i]/$srr[$i]_1.fastq.gz /path/to/read2/$srr[$i]/$srr[$i]_2.fastq.gz $srr[$i]_1P_clean.fq.gz $srr[$i]_1U_clean.fq.gz $srr[$i]_2P_clean.fq.gz $srr[$i]_2U_clean.fq.gz ILLUMINACLIP:/path/to/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:10:30 LEADING:28 TRAILING:28 MINLEN:75");


#Print the executed command to file
open($fh, '>', 'cmd.txt');
print $fh "\n\n$java -jar /path/to/trimmomatic-0.36.jar PE -trimlog $srr[$i]_trim_log /path/to/read1/$srr[$i]/$srr[$i]_1.fastq.gz /path/to/read2/$srr[$i]/$srr[$i]_2.fastq.gz $srr[$i]_1P_clean.fq.gz $srr[$i]_1U_clean.fq.gz $srr[$i]_2P_clean.fq.gz $srr[$i]_2U_clean.fq.gz ILLUMINACLIP:/path/to/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:10:30 LEADING:28 TRAILING:28 MINLEN:75\n\n";
close $fh;

	 $end = date;
	 system("$end > end.txt");
         chdir("$cwd");
	};

print "\n\nEND\n\n";
