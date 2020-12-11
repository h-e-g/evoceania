#!/bin/perl
#use warnings;

#Script to convert PLINK ped/bim files to the input file format of INPUT_to_Fst_AMOVA_indLevel.pl by Guillaume Laval (glaval@pasteur.fr)
#The script takes as input PLINK ped/bim files; Ped2Fst.pl takes one argument: the prefix of PLINK files
#The script is most efficient if parallelized. It typically takes as input a genomic window of 100Kb, using e.g. "plink --chr --from-bp $start --to-bp $end"


print "\t Ped2Fst.pl is running\n";

#Ped2Fst.pl expects one argument: the prefix of PLINK ped/bim files
if ($ARGV[0]) {
	@arg=split("/",$ARGV[0]);
	$pedfile=$arg[$#arg].".ped";
	$bimfile=$arg[$#arg].".bim";
	$outfile=$arg[$#arg]."_FST.inp";
	open (PEDFILE, "+<$pedfile") || die "\t    Error: cannot open $pedfile \n";
	open (BIMFILE, "+<$bimfile") || die "\t    Error: cannot open $bimfile \n";
	open (OUTFILE, ">$outfile") || die "\t    Error: cannot write in $outfile \n";
}
else {
	die "\t    Error: Ped2Fst.pl expects one argument: the prefix of PLINK input ped/bim files \n";
}


#Reads bim file
my $nb_SNPs=0;
my @pos=0;
my @rs=0;
while ($sline = <BIMFILE>){  
    chomp($sline);
    @splitted_sline = split("\t",$sline);
    $pos[$nb_SNPs]=$splitted_sline[3]; #SNP position
    $rs[$nb_SNPs]=$splitted_sline[1]; #rsID
	if (length($splitted_sline[4]) > 1) {
		$m[$nb_SNPs]=substr($splitted_sline[4], -1); # minor allele = last base if insertion
	}
	else {
		$m[$nb_SNPs]=$splitted_sline[4]; # minor allele
	}
	if (length($splitted_sline[5]) > 1) {
		$M[$nb_SNPs]=substr($splitted_sline[5], -1); # major allele = last base if insertion
	}
	else {
		$M[$nb_SNPs]=$splitted_sline[5]; # major allele
	}
    $nb_SNPs++;
}

#Writes positions and rs of SNPs, as rows
print "\t    Writing positions and rs of $nb_SNPs SNPs...\n";
print OUTFILE "\t\t\t";
for ($i=0; $i<$nb_SNPs-1;$i++){
    print OUTFILE $pos[$i]."\t";    
}
print OUTFILE $pos[$nb_SNPs-1]."\n\t\t#rs\t";

for ($i=0; $i<$nb_SNPs-1;$i++){
    print OUTFILE $rs[$i]."\t";
}
print OUTFILE $rs[$nb_SNPs-1]."\n\t\tMAF\t";

for ($i=0; $i<$nb_SNPs-1;$i++){
    print OUTFILE $m[$i]."\t";
}
print OUTFILE $m[$nb_SNPs-1]."\n\t\tMAJ\t";

for ($i=0; $i<$nb_SNPs-1;$i++){
    print OUTFILE $M[$i]."\t";
}
print OUTFILE $M[$nb_SNPs-1]."\n\t\tMAF\t";

for ($i=0; $i<$nb_SNPs-1;$i++){
    print OUTFILE $m[$i]."\t";
}
print OUTFILE $m[$nb_SNPs-1]."\n";

#Writes genotypes
print "\t    Writing individual genotypes...\n";
$b=0;
while ($bline = <PEDFILE>) {
       
	chomp($bline);
	@splitted_bline = split(" ",$bline);

	print OUTFILE $splitted_bline[0]."\t".$splitted_bline[4]."\t".$splitted_bline[1]."\t";
       
	for ($i=0; $i<$nb_SNPs-1; $i++){
		$s=2*$i+6;
		if (length($splitted_bline[$s]) > 1) {
			print OUTFILE substr($splitted_bline[$s], -1)."\t";
		}
		else {
			print OUTFILE $splitted_bline[$s]."\t";
		}
	}
	$s=2*$nb_SNPs+4;
	if (length($splitted_bline[$s]) > 1) {
		print OUTFILE substr($splitted_bline[$s], -1)."\n";
	}
	else {
		print OUTFILE $splitted_bline[$s]."\n";
	}
	
	print OUTFILE $splitted_bline[0]."\t".$splitted_bline[4]."\t".$splitted_bline[1]."\t";

	for ($i=0; $i<$nb_SNPs-1; $i++){
		$s=2*$i+7;
		if (length($splitted_bline[$s]) > 1) {
			print OUTFILE substr($splitted_bline[$s], -1)."\t";
		}
		else {
			print OUTFILE $splitted_bline[$s]."\t";
		}
	}
	$s=2*$nb_SNPs+5;
	if (length($splitted_bline[$s]) > 1) {
		print OUTFILE substr($splitted_bline[$s], -1)."\n";
	}
	else {
		print OUTFILE $splitted_bline[$s]."\n";
	}	
}
	
close OUTFILE;
close PEDFILE;
close BIMFILE;
close CONTROLFILE;
   
print "\t Ped2Fst.pl completed\n";
