#!/usr/bin/perl
use Cwd;


#Command to align multiple RNAseq Single read data to the Reference Genome using HISAT2
#created on 6 February 2018


#This code assumes that the current woking directory contains sub-directories
#Each sub-directory conatins data downloaded from European Nucleotide Archive(ENA)
#The aligned reads will be send to another directory ALIGN and named as per the SRR ids


#Create an array containing all the SRR ids
@srr=(SRR1234567, SRR1234568, SRR1234569, SRR1234560);


#Get current working directory and print it
$cwd = getcwd;
print "\n\n Current directory is $cwd\n\n";


#Initialize variable with the hisat2 command to be used as bash
$hisat = 'hisat2';


#for loop to change folder
for ($i=0; $i< @srr; $i++)
	{
   #Description of the work
	 print "\nRUNNING ALIGNMENT FOR $srr[$i]\n";
   
   #chdir PERL command to change directory	
	 chdir("$cwd/$srr[$i]");
   
   #system() is used to execute certain commands from PERL script to the TERMINAL
	 #Since hisat2 command cannot be directly executed as system(), a variable was used
	 #The necessary arguments for the command are fed in as usual
	 system ("$hisat -p 12 --dta -x /home/arindam/Work/NGS/ref_gen/Human_84/hisat2_index/grch38_tran/genome_tran -1 $srr[$i]_1.fastq.gz -2 $srr[$i]_2.fastq.gz  -S $cwd/align/$srr[$i].sam --un-gz $srr[$i]_unaln --summary-file $srr[$i]_sum --met-file $srr[$i]_met");

	 #report end time of alignment for each data set into a file end.txt	
	 $end = date;
	 system("$end > end.txt");
	 chdir("$cwd");	 
	};

print "\n\nEND\n\n";
