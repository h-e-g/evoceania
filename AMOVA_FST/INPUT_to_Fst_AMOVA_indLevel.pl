#!/bin/perl

#AMOVA-based FST computation by Guillaume Laval (glaval@pasteur.fr)
#The script takes as input a custom file format generated from PLINK files with the Ped2Fst.pl script
#The script is most efficient if parallelized. It typically takes as input a genomic window of 100Kb, using e.g. "plink --chr --from-bp $start --to-bp $end"
#The scripts takes three arguments: name of input file, name of output file, and chromosome


###################
#Function that checks if genotypes only include A/C/G/T or 0 bases

sub check_if_base {
	if( ($_[0] ne "A") && ($_[0] ne "T")  && ($_[0] ne "G")  && ($_[0] ne "C")  && ($_[0] ne "0")   ){
		print "check data:: unknown base (should be A, T, G, C or 0) \n";
		print "corrupted field is '" . $_[0] . "':  abort...";
		exit;
	}
}

###################
#Function that returns the number of alleles of a polymorphic site

sub check_numall {
	$p=$_[0]; $s=$_[1];
	$numall=0;
		if($nA_0{$p}{$s}{"A"}) { $numall++ };
		if($nA_0{$p}{$s}{"T"}) { $numall++ };
		if($nA_0{$p}{$s}{"G"}) { $numall++ };
		if($nA_0{$p}{$s}{"C"}) { $numall++ };	
	return $numall;
}

###################
#Function that computes sample size

sub compute_n {
	$p=$_[0];$s=$_[1];
	$res=$nA_0{$p}{$s}{"A"}+$nA_0{$p}{$s}{"T"}+$nA_0{$p}{$s}{"G"}+$nA_0{$p}{$s}{"C"};
	if(!$res){
		print "ERROR: empty SNP, pop=$p snp=$s\n";
	}
	return $res;
}

###################
#Function that computes minor allele frequency

sub find_MAF {
	$p=$_[0];$s=$_[1];	
	if( &check_numall($p,$s) == 1) {
		print "snp=$s monomorphic for the population $p\n";
		return (0,"?");
	}
	$maf=1000000;
	for($ii=0;$ii<4;$ii++){
		$curr=$nA_0{$p}{$s}{$base_def[$ii]};
		if( $curr ){
			if( $curr < $maf ) {
				$maf=$curr; $maf_all=$base_def[$ii];	
			}
		}
	}
	return ($maf,$maf_all);
}

###################
#Function that computes global maf (in merged populations)

sub find_freq_of_GLOBAL_MAF {
	$p=$_[0];$s=$_[1];	
	$prev=1000000;
	for($ii=0;$ii<4;$ii++){
		$curr=$nA_0{"ALL"}{$s}{$base_def[$ii]};
		if( $curr ){
			if( $curr < $prev ) {
				$global_maf_all=$base_def[$ii];
				$prev=$curr;		
			}
		}
	}
	return $nA_0{$p}{$s}{$global_maf_all};
}

###################
#Function that computes heterozygosity

sub compute_Het {
	$p=$_[0];$s=$_[1];$homo=0;
	if( &compute_n($p,$s) ){
		$corr= &compute_n($p,$s) / ( &compute_n($p,$s) - 1 );
		for($ii=0;$ii<4;$ii++){
			$homo+=( $nA_0{$p}{$s}{$base_def[$ii]} / &compute_n($p,$s) )**2;
		}
		$het=$corr*(1-$homo);
	return $het;
	}
	else{
		return "NC";
	}
}

###################

sub check_statu_SNP {
	if( ($_[0] ne "SN") && ($_[0] ne "IN")  && ($_[0] ne "0")  && ($_[0] ne "A")  && ($_[0] ne "T")   && ($_[0] ne "G")  && ($_[0] ne "C")  ){
		print "check data:: unknown SNP status (should be SN, IN or 0) \n";
		print "corrupted field is '" . $_[0] . "': abort...";
		exit;
	}
}

###################
#Function that translates numbers into bases

sub base_translator {
	if( $_[0] eq 0){
		return "0";
	}
	elsif( $_[0] eq 1){
		return "A";
	}
	elsif( $_[0] eq 2){
		return "T";
	}
	elsif( $_[0] eq 3){
		return "G";
	}
	elsif( $_[0] eq 4){
		return "C";
	}
}

###################
#Function that writes the header of the output file

sub Header_statfile {
	print STATFILE  "CHR\trs\tSNP_pos\t";
	
	#foreach my $pop (keys(%hashPOP)) { #total number of chromosomes
	#	print STATFILE  "$pop"."_n\t";
	#}
	#foreach my $pop (keys(%hashPOP)) { #population minor allele count
	#	print STATFILE  "$pop"."_nA\t";
	#}
	
	$nbpop=0; @pop_order=0;
	foreach my $pop (keys(%hashPOP)) { #population minor allele frequency
		print STATFILE  "$pop"."_fA\t";
		$pop_order[$nbpop]=$pop;
		$nbpop++;
	}
	
	#foreach my $pop (keys(%hashPOP)) { #heterozygosity MAF
	#	print STATFILE  "$pop"."_HetA\t";
	#}
	#print STATFILE  "GLOB_n\t"; #total number of chromosomes
	#print STATFILE  "GLOB_MAF\t"; #global MAF
	
	print STATFILE  "GLOB_minor_allele\t"; #global minor allele	
	print STATFILE  "GLOB_MAF\t"; #global MAF
	
	#print STATFILE  "GLOB_HetA\t"; #global MAF	
	#foreach my $pop (keys(%hashPOP)) { #  Hst == Population Fst based on Laurent Excoffier's method 
	#	print STATFILE  "$pop"."_Hst\t"; #NOT IMPLEMENTED
	#}
	print STATFILE  "GLOB_Fst\t"; #GLOBAL Fst
	
	for ($i=0;$i<$nbpop;$i++) {
		for ($j=$i+1;$j<$nbpop;$j++) {
			print STATFILE  $pop_order[$i]."_".$pop_order[$j]."_Fst\t";
		}
	}
	
	#print STATFILE  "Country_Fst\t"; #GLOBAL FST computed using country definition
	#foreach my $pop (keys(%hashPOP)) { # population Fst (continent here) where subdivion is by country
	#	print STATFILE  "$pop"."_Fst\t";
	#}

	print STATFILE  "\n";
}

###################
# Function that outputs the row of computed statistics

sub result_line {

	$iSNP=$_[0]; $chr=$_[1];
	$snp=$rsSNP[$iSNP]; $location=$pos[$iSNP];
	
	$format_double="%8.6f";
	print STATFILE  "$chr\t$snp\t$location\t";
	
	#foreach my $pop (keys(%hashPOP)) { #total number of chromosomes
	#	print STATFILE  &compute_n($pop,$snp) . "\t";
		#print STATFILE  $n{$pop}{$snp} . "\t"; 
	#}
	#foreach my $pop (keys(%hashPOP)) { #Frequency in population of the Global MAF
	#	print STATFILE  &find_freq_of_GLOBAL_MAF($pop,$snp) . "\t";
		#($a,$b)=&find_MAF($pop,$snp);
		#print STATFILE  "$b\t";
	#}
	
	$nbpop=0; @pop_order=0;
	foreach my $pop (keys(%hashPOP)) { #relative frequency in population of the Global MAF
		$pop_order[$nbpop]=$pop;
		$nbpop++;	
		if( !&compute_n($pop,$snp) ){ #empty continent(remove from calculation)
			print STATFILE "NC\t";
		}
		else{
			printf STATFILE "$format_double\t",  &find_freq_of_GLOBAL_MAF($pop,$snp) / &compute_n($pop,$snp);
		}
	}
	
	#foreach my $pop (keys(%hashPOP)) { #heterozgosity MAF
	#	printf STATFILE "$format_double\t",  &compute_Het($pop,$snp) ;
	#}
	#print STATFILE  $n{"ALL"}{$snp} . "\t"; #total number of chromosomes
	#print STATFILE  &compute_n("ALL",$snp) . "\t"; #total number of chromosomes
	
	($a,$b)=&find_MAF("ALL",$snp);	
	#print STATFILE  "$a\t"; #global MAF
	#printf STATFILE "$format_double\t",  $a/$n{"ALL"}{$snp}; #global MAF
	print STATFILE  "$b\t"; #global MAF
	printf STATFILE "$format_double\t",  $a/&compute_n("ALL",$snp); #global MAF
	#printf STATFILE "$format_double\t",  &compute_Het("ALL",$snp) ; #global Het	
	#foreach my $pop (keys(%hashPOP)) { #Population Fst
	#	print STATFILE  " \t";
	#}
	
	printf STATFILE "$format_double\t",  &multiPOP_Fst($snp,\%hashPOP) ; #GLOBAL Fst

	for ($i=0;$i<$nbpop;$i++) {
		for ($j=$i+1;$j<$nbpop;$j++) {
			$temphashPOP{$pop_order[$i]}=$pop_order[$i];
			$temphashPOP{$pop_order[$j]}=$pop_order[$j];
			printf STATFILE "$format_double\t",  &multiPOP_Fst($snp,\%temphashPOP);
			%temphashPOP=();
		}
	}

	#printf STATFILE "$format_double\t",  &multiPOP_Fst($snp,\%hashCOUNTRY) ; #GLOBAL FST computed using country def	
	#foreach my $pop (keys(%hashPOP)) { # population Fst (continent here) where subdivion is by country
	#	printf STATFILE "$format_double\t",  &multiPOP_Fst($snp,\%{$hashPOP_COUNTRY{$pop}}) ;
	#}
	
	print STATFILE  "\n";
}

###################

sub combi_freq {
	$p_snp=$_[0]; $ref_hashPOP=$_[1]; ($p_maf,$p_maf_all)=&find_MAF("ALL",$snp);

	$n=0; $n_all=0;$phet=0;	
	foreach my $pop (keys(%$ref_hashPOP)) { 
		$n+=&compute_n($pop,$snp); $n_all+=$nA_0{$pop}{$snp}{$p_maf_all};
	}	
	$pp=$n_all/$n;
	$corr= $n / ( $n - 1 );$phet=2*$pp*(1-$pp);
	$phet=$corr*($phet);
	return ($pp,$phet);
}

###################
#Function that computes AMOVA FST (multiple populations)

sub multiPOP_Fst {         #AMOVA  #Checked vs Arlequin
	#Passing arguments table or hash table by reference \%hashPOP \@listPOP
	$p_snp=$_[0]; $ref_hashPOP=$_[1]; ($p_maf,$p_maf_all)=&find_MAF("ALL",$snp);
	
	#SST all populations 
	$SST=0; $npop=0; $n=0; $n_1=0; $n_square=0; $n_A1=0;
	foreach my $pop (keys(%$ref_hashPOP)) { 
		$n_i=&compute_n($pop,$snp); $n_A1+=$nA_0{$pop}{$snp}{$p_maf_all};
		$n+= $n_i; $n_1+= ( $n_i - 1 ); $n_square+= $n_i*$n_i;		
	}
	if( !$n ){ #in case all countries from the same continent are empty
		return "NC";
	}
	
	$SST= ( $n_A1 * ( $n - $n_A1 ) ) / $n;
	
	$SSWP=0; $n_A1=0;
	foreach my $pop (keys(%$ref_hashPOP)) { 
		if( &compute_n($pop,$snp) ){ #empty continent or country (remove from calculation)
			$n_i=&compute_n($pop,$snp); $n_A1=$nA_0{$pop}{$snp}{$p_maf_all};		
			$tempSSWP=( $n_A1 * ( $n_i - $n_A1 ) ) / $n_i;		
			$SSWP+=$tempSSWP;		
			
			$npop+=1; 
		}
	}	
		
	if($npop <= 1){
		return "NC";
	}
	
	$n_c=( $n  -  ( $n_square / $n ) ) / ( $npop - 1 ) ;
	$b=$n_c*($npop-1);	
	
	$N_hetero=0;
	foreach my $pop (keys(%$ref_hashPOP)) { 
		$N_hetero+=$P_hetero{$pop}{$snp};
	}
	
	
	$SSAP=$SST - $SSWP; $SSWI=$N_hetero / 2; $SSAIWP=$SSWP - $SSWI;	
	
	$MSA = $SSAP / ( $npop - 1  );

	$MSB = $SSAIWP / ( ( $n / 2 ) - $npop );
	
	$MSC = $SSWI / ( $n / 2);

	$sigmaC = $MSC;
	$sigmaB = ( $MSB - $MSC ) * 0.5;
	$sigmaA = ( $MSA - $MSB ) / $n_c ; 
	$sigmaT = $sigmaA + $sigmaB + $sigmaC;
	
	
	$numerator=( $sigmaA ); $denominator=( $sigmaT  );
	if( $denominator != 0 ){
		$fst= $numerator / $denominator ;
	}
	else{
		$fst="NC";
	}
	
	return $fst;
}

###################
#Function that computes AMOVA FST / not used

sub pop_Fst {	
	$first_allele=&get_A($_[0],$iSNP); #param ipop and iSNP
	
	#SSWP all populations
	$SST=0;         $npop=4;     $n=0; $n_1=0; $n_square=0; $n_A1=0;
	for($ipop=0 ;$ipop<$npop;$ipop++){
		$n_i= ( $nA[$ipop][$iSNP][0] + $nA[$ipop][$iSNP][1] );
		$n_A1+=&get_n_A($first_allele,$ipop,$iSNP);
		$n+= $n_i; $n_1+= ( $n_i - 1 ); $n_square+= $n_i*$n_i;
	}	
	$SST+= ( $n_A1 * ( $n - $n_A1 ) ) / $n;
	
	#SSWP
	$SSWP=0; $n_A1=0;
	for($ipop=0 ;$ipop<$npop;$ipop++){
		
		$tempSSWP=0;
		$n_i= ( $nA[$ipop][$iSNP][0] + $nA[$ipop][$iSNP][1] );
		
		$n_A1=&get_n_A($first_allele,$ipop,$iSNP);
		
		$tempSSWP=( $n_A1 * ( $n_i - $n_A1 ) ) / $n_i;
		if($ipop == $_[0] ){
			$SSWP_i=( $n_A1 * ( $n_i - $n_A1 ) ) / $n_i;
		}		
		$SSWP+=$tempSSWP;
	}
	
	$n_c=( $n  -  ( $n_square / $n ) ) / ( $npop - 1 ) ;
	$b=$n_c*($npop-1);
	
	$SSAP=$SST - $SSWP;
	
	$MSA = $SSAP / ( $npop - 1  );
	$MSB = $SSWP / ( $n - $npop );
	
	$sigmaB = $MSB;
	$sigmaA = ( $MSA - $MSB ) / $n_c ; 
	$sigmaT = $sigmaA + $sigmaB ;
	
	if( $rsSNP[$iSNP] eq $ckeckedSNP ) {
		print "\nSST='$SST', SSAP='$SSAP', SSWP='$SSWP'";
	}
	
	$n_i= ( $nA[$_[0]][$iSNP][0] + $nA[$_[0]][$iSNP][1] );
	$numerator=( $SSAP / ( $npop - 1 )   - ( 1 / $n_i )*( $n / ( $n - $npop ) )*$SSWP_i  ) / ( $n_c );
	$denominator=( $sigmaT );
	if( $denominator != 0 ){
		$fst= $numerator / $denominator ;
	}
	else{
		$fst="NC";
	}
	

}

###################
#Reads input files

open (FILE, "+<$ARGV[0]") || die "cannot open $ARGV[0] ";
	
	$base_def[0]="A";$base_def[1]="T";$base_def[2]="G";$base_def[3]="C";

	$line 	= <FILE>;	chomp($line);  @cell =split(/\t/, $line); #name of marker
	$NumSNP	= (@cell-3); #number of markers
	$line1 	= <FILE>;	chomp($line1); @cell1=split(/\t/,$line1); #status of marker SN or IN/DEL
	$line2 	= <FILE>;	chomp($line2); @cell2=split(/\t/,$line2); #code 0 
	$line3 	= <FILE>;	chomp($line3); @cell3=split(/\t/,$line3); #code 1
	$line4 	= <FILE>;	chomp($line4); @cell4=split(/\t/,$line4); #ancestral polymorphism	
	for($i=3;$i<($#cell+1);$i++){	
			$pos[$i]	=  $cell[$i] ;    chomp($pos[$i]);
			$rsSNP[$i]	= $cell1[$i]; chomp($rsSNP[$i]);
			$glob_MAF[$i]	= uc($cell2[$i]); chomp($glob_MAF[$i]); &check_if_base($glob_MAF[$i]);
			$glob_MAJ[$i]	= uc($cell3[$i]); chomp($glob_MAJ[$i]); &check_if_base($glob_MAJ[$i]);
			$ancestor[$i]	= uc($cell4[$i]); chomp($ancestor[$i]); &check_if_base($ancestor[$i]);	
	}
	
	for($i=0;$i<($#pos);$i++){
		if( $pos[$i] > $pos[$i+1]  ){
			print "check data:: SNPs are not in the right order \n"; print "abort..."; exit;
		}		
	}	
	print "\t\tExtracting data ... \n";
	$prev="start"; $diploid=1; 
	while($genot_line = <FILE>){
		chomp($genot_line); @cell=split(/\t/, $genot_line);
		$pop=$cell[0]; chomp($pop);
		$country="ctry_" . $cell[1]; chomp($country);
		$ind=$cell[2]; chomp($ind);
		
		$hashPOP{$pop}=$pop;
		$nchr{$pop}++; #ideal number of chromosomes
		$hashIND{$pop}{$ind}=$ind;
		
		
		$hashCOUNTRY{$country}=$country;
		$hashPOP_COUNTRY{$pop}{$country}=$country;
		$hashIND_country{$country}{$ind}=$ind;		
		
		if($ind eq $prev){
			$all="q"; $diploid=1; 
		}
		else{
			if( !$diploid ){
				print "\t WARNING :: individual $prev is not diploid, abort ....\n";exit;
			}
			$all="p"; $diploid=0; 
		}
		
		
		for($i=3;$i<($#cell+1);$i++){
				$tempall=$cell[$i] ; chomp($tempall); &check_if_base($tempall);				
				$snp=$rsSNP[$i];				
				$allele{$pop}{$ind}{$snp}{$all}=$tempall;
				#print "'$pop'\t'$snp'\t'$ind'\t$all\t'$tempall'\t";
				#print "'$allele{$pop}{$ind}{$snp}{$all}'\n";
				
				$allele_country{$country}{$ind}{$snp}{$all}=$tempall;
				
							
		}
		$prev=$ind;
	}
	#print "\t\tdone\n";
	#print "NumSNP=$NumSNP\n";
close FILE;
	
	for($isnp=3 ;$isnp<$NumSNP+3;$isnp++){
		$snp=$rsSNP[$isnp];
		foreach my $pop (keys(%hashPOP)) {
			$P_hetero{$pop}{$snp}=0;
			$nA_0{$pop}{$snp}{"A"}=0;$nA_0{$pop}{$snp}{"T"}=0;$nA_0{$pop}{$snp}{"G"}=0;
			$nA_0{$pop}{$snp}{"C"}=0;$nA_0{$pop}{$snp}{"0"}=0;			
		}
		$P_hetero{"ALL"}{$snp}=0;
		$nA_0{"ALL"}{$snp}{"A"}=0;$nA_0{"ALL"}{$snp}{"T"}=0;$nA_0{"ALL"}{$snp}{"G"}=0;
		$nA_0{"ALL"}{$snp}{"C"}=0;$nA_0{$pop}{$snp}{"0"}=0;
	}
	
	for($isnp=3 ;$isnp<$NumSNP+3;$isnp++){
		$snp=$rsSNP[$isnp];
		foreach my $pop (keys(%hashPOP)) {
			foreach my $ind (keys(  %{$hashIND{$pop}}  ) ) {	
				$base1=$allele{$pop}{$ind}{$snp}{"q"};			
				$base2=$allele{$pop}{$ind}{$snp}{"p"};
				$nA_0{$pop}{$snp}{$base1}++;	
				$nA_0{$pop}{$snp}{$base2}++;	
				$nA_0{"ALL"}{$snp}{$base1}++;	
				$nA_0{"ALL"}{$snp}{$base2}++;	
				if($base1 ne $base2){
					$P_hetero{$pop}{$snp}++;
					$P_hetero{"ALL"}{$snp}++;
				}
			}
		}
		
		#TO COMPUTE FST BETWEEN COUNTRIES (INDEPENDENT COMPUTATION)
		foreach my $country (keys(%hashCOUNTRY)) {
			foreach my $ind (keys(  %{$hashIND_country{$country}}  ) ) {	
				$base1=$allele_country{$country}{$ind}{$snp}{"q"};			
				$base2=$allele_country{$country}{$ind}{$snp}{"p"};
				$nA_0{$country}{$snp}{$base1}++;	
				$nA_0{$country}{$snp}{$base2}++;	
				if($base1 ne $base2){
					$P_hetero{$country}{$snp}++;
				}
			}
		}
		
		
	}
	
	print "\t\tComputing summary statistics for $NumSNP SNPs .... \n";
	
#Writes output file
open (STATFILE, "+>$ARGV[1]") || die "cannot open $ARGV[1] ";
	$chr=$ARGV[2];
	&Header_statfile();
	$numMonomorphicSNP=0;$numTrimorphicSNP=0;
	for($isnp=3 ;$isnp<$NumSNP+3;$isnp++){
		$snp=$rsSNP[$isnp];
		$n=&check_numall("ALL",$snp);
		if($n == 0){
			$fulllEMPTY++,
		}
		elsif($n == 1){
			$numMonomorphicSNP++,
		}
		elsif($n == 3){
			$numTrimorphicSNP++,
		}				
		elsif( $n == 2){
			&result_line($isnp,$chr);
		}
	}
	print "\t\tAnalysis of " . ($NumSNP-($numMonomorphicSNP + $numTrimorphicSNP))  . " SNPs: ";
	print "$fulllEMPTY fully empty, $numMonomorphicSNP monomorphic and $numTrimorphicSNP trimorphic positions were excluded \n";
close STATFILE;	
