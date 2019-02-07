###############################################################
#  Perl script to extract sample source and other             #
#    information from NCBI BioSample database                 #
#    for a given list of sample ids.                          #
#                                                             #
#                     - Arindam Ghosh (7 February 2019)       #
###############################################################

#!/usr/bin/perl
use Cwd;

$cwd = getcwd;

#SampleIds.txt contains a list of sample ids (eg: SAMN04383981)
open (FH1, "SampleIds.txt");
@sampleids = <FH1>;
close FH1;


open (FH2, '>', "output.txt");
print FH2 "Sl.no.\tSample_ID\tSample_Title\tOrganismName\tSample_Source\n";


for($i=1; $i<@sampleids; $i++)
{

chop($sampleids[$i]);

$title_cmd = "efetch -db biosample -id $sampleids[$i] -format xml | xtract -pattern BioSampleSet -division BioSample -block Description -element Title Organism/OrganismName";
$title = qx/$title_cmd/;
chop($title);

$source_name_cmd = "efetch -db biosample -id $sampleids[$i] -format xml | xmllint --xpath \"/BioSampleSet/BioSample/Attributes/Attribute[\@attribute_name='source_name']/text()\" -";
$source_name = qx/$source_name_cmd/;

print FH2 "$i\t$sampleids[$i]\t$title\t$source_name\n";

};

close FH2;




