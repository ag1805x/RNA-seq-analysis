#!/usr/bin/perl
use Cwd;

#Command to align multiple RNAseq Single read data to the Reference Genome using HISAT2


#This code assumes that the current woking directory contains sub-directories serial numbered 1-16
#Each sub-directory conatins data downloaded from European Nucleotide Archive(ENA)
#The aligned reads will be send to another directory ALIGN and named as per the number/name 
#used in the sub-directory containing the original data


#This code is exactly what I had used for my system and may need editing PATH as per requiement



#Get current working directory and print it
$cwd = getcwd; 
print "\n\nCWD is $cwd\n\n";


#Description of the work
print "\n\nRUNNING READ ALIGNMENT\n\n";


#Initialize variable with the hisat2 command to be used as bash
$hisat='hisat2';


#for loop to change folder
for ($i=5;$i<=16;$i++)
{
	print "\n\nRUNNING READ ALIGNMENT for $i\n\n";	
	#chdir PERL command to change directory	
	chdir("$cwd/$i");	
	#system() is used to execute certain commands from PERL script to the TERMINAL
	#Since hisat2 command cannot be directly executed as system(), a variable was used
	#The necessary arguments for the command are fed in as usual
	#The alignment summary is saved as align_summary.txt file
	system("$hisat -p 12 --dta -x /home/arindam/Work/NGS/ref_gen/Human_84/hisat2_index/grch38_tran/genome_tran -U /home/arindam/Work/NGS/Data/$i/*.fastq.gz -S /home/arindam/Work/NGS/Data/align/$i.sam 2> $cwd/$i/align_summary.txt");
	#report end time of alignment for each data set into a file end.txt	
	$end = date;
	system("$end > $cwd/$i/end.txt");
	chdir("$cwd");
};
#report final end time
$time = localtime();
print "\n\n\nEND TIME  IS $time\n\n\n";
