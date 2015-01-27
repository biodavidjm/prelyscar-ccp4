#!/usr/bin/perl -w

use POSIX;
use strict;
use warnings;


# David Jimenez-Morales
# @ davidjm scripts
# Release date: 3/1/2014

# ------Description ---------------------------------------#
# prelyscarATccp4.pl
#
# Predictor of Lysine Carboxylation
#
# Adapted for CCP4
# 3 inputs required: PDB + CHAIN + PRIOR
#
# ---------------------------------------------------------#

my $command = "prelyscarATccp4.pl";
if (@ARGV != 3 )
{
	print "ERROR !!!!!!!!!!!!!!!\n";
	print $command." requires 3 arguments (PDB file + PROTEIN_CHAIN + PRIOR )\n";
	exit;
}

#INPUT 1, the PDB
my $filename = $ARGV[0];
chomp($filename);

my $PDB = '';
if ( ($filename =~ /^(\w{4}).pdb$/) )
{
	$PDB = $1;
}
else
{
	print "ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	print "The first argument must be a PDB.pdb\n";
	print $command." requires 3 arguments (PDB file + PROTEIN_CHAIN + PRIOR )\n";
	exit;
}

# INPUT 2, the protein chain
my $chain = $ARGV[1];
chomp($chain);
unless ($chain =~ /^\w$/)
{
	print "ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	print "The protein chain is missed\n";
	print $command." requires 3 arguments (PDB file + PROTEIN_CHAIN + PRIOR )\n";
	exit;
}

my $prior = $ARGV[2];
chomp($prior);
unless ($prior =~ /^0\.?\d+$/)
{
	print "ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	print "The prior probability must be a number between 0 and 0.99999\n";
	print $command." requires 3 arguments (PDB file + PROTEIN_CHAIN + PRIOR )\n";
	exit;
}

#print "Arguments: ".$PDB." ".$chain." ".$prior."\n";

print "\n************************************\n";
print "P r e L y s C a r \n";
print "____________________________________\n";
print "Predictor of Lysine Carboxylation\n";
print "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-\n";

my $pdb = $PDB;
my $aatofind = "LYS";
my $distance = 5;

# Do we need to INCREASE DISTANCE for testing purpuses?
# This is the place...
my $extradistance = $distance;

# Loading PDB file
my @opdb = get_file_data($filename);

# Checking number of residues in the protein
my %pdbnumres = ();

# Control: the input is only one PDB chain
my $initialchain = ''; 
my $c = 0;

foreach my $pdbline (@opdb)
{
	chomp($pdbline);
	# GETTING THE PROTEIN SEQUENCE
	if ($pdbline =~ /^ATOM/)
	{
		my $ithaschain = substr($pdbline,21,1);
		
		if($chain eq $ithaschain)
		{
			$c++;
			my $resnum = substr ($pdbline,22,4);
			$resnum =~ s/\s//g;
			
			my $aa3 = substr($pdbline,17,3);
			
			if(!$pdbnumres{$resnum})
			{
				$pdbnumres{$resnum} = $aa3
			}
		}
	} #if it's an ATOM
	
}#foreach

# Control
unless ($c != 0)
{
	print "\nThis chain does not contain amino acids\n";
	print "Try another one\n\n";
	exit;
}

# PROTEIN LENGTH: All the proteins known
my $protlen = keys %pdbnumres;
print "\nPDB ID: ".$pdb.", chain:$chain length: ".$protlen." amino acids\n";

my $protlencode = 0;

if ($protlen >= 200)
{
	$protlencode = 2;
#	print "Protein Length: ".$protlen." (more than 200)\n\n";
}
elsif ($protlen < 199)
{
	$protlencode = 1;
#	print "Protein Length: ".$protlen." (less than 200)\n\n";
}
else
{
	# Control: if there is no length, there is no protein
	print "\nFATAL ERROR !!!!!!!!!!!!!!!\n";
	print "Error:$command cannot determine the length of the protein in this chain\n\n";
	exit;
}


# LOADPDB COORDINATES FOR THE SELECTED AA + ALL THE OTHER
my ($aacoor, $allothercoor) = get_pdbsidechain_coordinateskcxcase($filename,$aatofind,$chain);
my %aacoordinates = %$aacoor;
my %alltheothercoordinates = %$allothercoor;


my $naa = keys %aacoordinates;
my $not = keys %alltheothercoordinates;

if ( ($naa == 0) | ($not == 0) )
{
	print "\nThere is no $aatofind in this subunit ($pdb)\n";
	print "For any question/suggestion/comment contact us at prelyscar\@gmail.com\n\n";
	exit;
}

# GETTING THE ATOMIC COORDINATES OF IONS AND WATER
my ($ionhetatmcoor, $watercoor) = Get_PDBHETATM_Coordinates_Ion_Water($filename,$chain);
my %ionhetatmcoor = %$ionhetatmcoor;
my %watercoordinates = %$watercoor;

my %listofaasandnames = ();
# $listofaasandnames{AANUMBER} = AANAME;
my %listofnearaminoacidstothechoseone = ();
# $listofnearaminoacidstothechoseone{CHOSEN AA NUMBER}{UNIQUE AA CLOSE BY} = BFACTOR;
my %averageeuclideans = ();
# $averageeuclideans{AANUMBER} = AVEEUCLIDEAN

for my $thechosenaa (sort {$a<=>$b} keys %aacoordinates)
{
	my $euclideanaverage = 0;
	my $eutotal = 0;
	for my $atomoftheone (keys %{$aacoordinates{$thechosenaa}})
	{
		#@tmp = ($resnum,$aa3,$x,$y,$z,$bfactor,$atomname);
		my $aasN = $aacoordinates{$thechosenaa}{$atomoftheone}[0];
		my $aa3 = $aacoordinates{$thechosenaa}{$atomoftheone}[1];
		my $kx = $aacoordinates{$thechosenaa}{$atomoftheone}[2];
		my $ky = $aacoordinates{$thechosenaa}{$atomoftheone}[3];
		my $kz = $aacoordinates{$thechosenaa}{$atomoftheone}[4];
		my $atomname = $aacoordinates{$thechosenaa}{$atomoftheone}[6];
		
		#FIRST GO AROUND ALL THE OTHER RESIDUES IN THE PROTEIN
		for my $theotheraanumber (sort {$a<=>$b} keys %alltheothercoordinates)
		{
			for my $theotheratom (keys %{$alltheothercoordinates{$theotheraanumber}})
			{
				my $oaan = $alltheothercoordinates{$theotheraanumber}{$theotheratom}[0];
				my $oaa3 = $alltheothercoordinates{$theotheraanumber}{$theotheratom}[1];
				my $ox = $alltheothercoordinates{$theotheraanumber}{$theotheratom}[2];
				my $oy = $alltheothercoordinates{$theotheraanumber}{$theotheratom}[3];
				my $oz = $alltheothercoordinates{$theotheraanumber}{$theotheratom}[4];
				my $bfactor = $alltheothercoordinates{$theotheraanumber}{$theotheratom}[5];
				my $otheratom = $alltheothercoordinates{$theotheraanumber}{$theotheratom}[6];
				my $eudistance = EuclideanDistance($kx,$ky,$kz,$ox,$oy,$oz);
				if ($atomname eq "PCX")
				{
					if ($eudistance <= $extradistance)
					{
						if (!$listofnearaminoacidstothechoseone{$thechosenaa}{$theotheraanumber})
						{			
							$listofnearaminoacidstothechoseone{$thechosenaa}{$theotheraanumber} = $bfactor;
							if (!$listofaasandnames{$oaan})
							{
								$listofaasandnames{$oaan} = $oaa3;
							}
							$euclideanaverage += $eudistance;
							$eutotal++;
							last;	
						}
					}
				}
				else
				{
					if ($eudistance <= $distance)
					{
						if (!$listofnearaminoacidstothechoseone{$thechosenaa}{$theotheraanumber})
						{
							$listofnearaminoacidstothechoseone{$thechosenaa}{$theotheraanumber} = $bfactor;
							if (!$listofaasandnames{$oaan})
							{
								$listofaasandnames{$oaan} = $oaa3;
							}
							$euclideanaverage += $eudistance;
							$eutotal++;
							last;	
						}
					}
				}

			}
		}#END OF CHECKING AROUND "ALL" THE OTHER RESIDUES
		
		#CHECK AROUND YOUR OWN RESIDUES
		for my $thechosenaa2 (sort {$a<=>$b} keys %aacoordinates)
		{
			if ($thechosenaa2 == $thechosenaa)
			{
				# We don't want to compare the same Chosen Amino Acid
				next;
			}
			else
			{
				for my $atomoftheone2 (keys %{$aacoordinates{$thechosenaa2}})
				{
					my $aasN2 = $aacoordinates{$thechosenaa2}{$atomoftheone2}[0];
					my $aa32 = $aacoordinates{$thechosenaa2}{$atomoftheone2}[1];
					my $kx2 = $aacoordinates{$thechosenaa2}{$atomoftheone2}[2];
					my $ky2 = $aacoordinates{$thechosenaa2}{$atomoftheone2}[3];
					my $kz2 = $aacoordinates{$thechosenaa2}{$atomoftheone2}[4];
					my $bfactor = $aacoordinates{$thechosenaa2}{$atomoftheone2}[5];
					my $otheratom = $aacoordinates{$thechosenaa2}{$atomoftheone2}[6];
					
					my $eudistance = EuclideanDistance($kx,$ky,$kz,$kx2,$ky2,$kz2);

						if ($atomname eq "NZ")
						{
							if ($eudistance <= $extradistance)
							{
								if (!$listofnearaminoacidstothechoseone{$thechosenaa}{$thechosenaa2})
								{			
									$listofnearaminoacidstothechoseone{$thechosenaa}{$thechosenaa2} = $bfactor;
									if (!$listofaasandnames{$aasN2})
									{
										$listofaasandnames{$aasN2} = $aa32;
									}
									$euclideanaverage += $eudistance;
									$eutotal++;
									last;	
								}
							}
						}
						else
						{
							if ($eudistance <= $distance)
							{
								if (!$listofnearaminoacidstothechoseone{$thechosenaa}{$thechosenaa2})
								{
									$listofnearaminoacidstothechoseone{$thechosenaa}{$thechosenaa2} = $bfactor;
									if (!$listofaasandnames{$aasN2})
									{
										$listofaasandnames{$aasN2} = $aa32;
									}
									$euclideanaverage += $eudistance;
									$eutotal++;
									last;	
								}
							}					
						}
				
				
				} # for my $atomoftheone2 (keys %{$aacoordinates{$thechosenaa2}})
			}
		} # END OF CHECKING AROUND SIMILAR RESIDUES
		
	
		#-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
		# CHECKING WATER MOLECULES AROUND:
		for my $theotheraanumber (sort {$a<=>$b} keys %watercoordinates)
		{
			for my $theotheratom (keys %{$watercoordinates{$theotheraanumber}})
			{
				my $oaan = $watercoordinates{$theotheraanumber}{$theotheratom}[0];
				my $oaa3 = $watercoordinates{$theotheraanumber}{$theotheratom}[1];
				my $ox = $watercoordinates{$theotheraanumber}{$theotheratom}[2];
				my $oy = $watercoordinates{$theotheraanumber}{$theotheratom}[3];
				my $oz = $watercoordinates{$theotheraanumber}{$theotheratom}[4];
				my $bfactor = $watercoordinates{$theotheraanumber}{$theotheratom}[5];
				my $otheratom = $watercoordinates{$theotheraanumber}{$theotheratom}[6];
				my $eudistance = EuclideanDistance($kx,$ky,$kz,$ox,$oy,$oz);
				if ($atomname eq "PCX")
				{
					if ($eudistance <= $extradistance)
					{
						if (!$listofnearaminoacidstothechoseone{$thechosenaa}{$theotheraanumber})
						{			
							$listofnearaminoacidstothechoseone{$thechosenaa}{$theotheraanumber} = $bfactor;
							if (!$listofaasandnames{$oaan})
							{
								$listofaasandnames{$oaan} = $oaa3;
							}
							last;	
						}
					}
				}
				else
				{
					if ($eudistance <= $distance)
					{
						if (!$listofnearaminoacidstothechoseone{$thechosenaa}{$theotheraanumber})
						{
							$listofnearaminoacidstothechoseone{$thechosenaa}{$theotheraanumber} = $bfactor;
							if (!$listofaasandnames{$oaan})
							{
								$listofaasandnames{$oaan} = $oaa3;
							}
							last;	
						}
					}
				}
			}
		}#END OF CHECKING AROUND "ALL" THE WATER MOLECULES!!
		#-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
		
		
		#-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
		# CHECKING AROUND IONS:
		for my $theotheraanumber (sort {$a<=>$b} keys %ionhetatmcoor)
		{
			for my $theotheratom (keys %{$ionhetatmcoor{$theotheraanumber}})
			{
				my $oaan = $ionhetatmcoor{$theotheraanumber}{$theotheratom}[0];
				my $oaa3 = $ionhetatmcoor{$theotheraanumber}{$theotheratom}[1];
				my $ox = $ionhetatmcoor{$theotheraanumber}{$theotheratom}[2];
				my $oy = $ionhetatmcoor{$theotheraanumber}{$theotheratom}[3];
				my $oz = $ionhetatmcoor{$theotheraanumber}{$theotheratom}[4];
				my $bfactor = $ionhetatmcoor{$theotheraanumber}{$theotheratom}[5];
				my $otheratom = $ionhetatmcoor{$theotheraanumber}{$theotheratom}[6];
				my $eudistance = EuclideanDistance($kx,$ky,$kz,$ox,$oy,$oz);
				if ($atomname eq "PCX")
				{
					if ($eudistance <= $extradistance)
					{
						if (!$listofnearaminoacidstothechoseone{$thechosenaa}{$theotheraanumber})
						{			
							$listofnearaminoacidstothechoseone{$thechosenaa}{$theotheraanumber} = $bfactor;
							if (!$listofaasandnames{$oaan})
							{
								$listofaasandnames{$oaan} = $oaa3;
							}
							last;	
						}
					}
				}
				else
				{
					if ($eudistance <= $distance)
					{
						if (!$listofnearaminoacidstothechoseone{$thechosenaa}{$theotheraanumber})
						{
							$listofnearaminoacidstothechoseone{$thechosenaa}{$theotheraanumber} = $bfactor;
							if (!$listofaasandnames{$oaan})
							{
								$listofaasandnames{$oaan} = $oaa3;
							}
							last;	
						}
					}
				}

			}
		}#END OF CHECKING AROUND "ALL" THE IONSSSSS!!
		#-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
	}#for my $atomoftheone (keys %{$aacoordinates{$thechosenaa}})

	if ( ($euclideanaverage>0) & ($eutotal>0) )
	{
		my $euavevalue = ($euclideanaverage/$eutotal);
		$averageeuclideans{$thechosenaa} = $euavevalue;
	}
}

#  p r e l y s c a r oooooooooooooooooooooooooooooooout
# OUTPUT FILE: Pymol script
my $outpymol = $pdb.".chain".$chain.".euclidean".$distance."-prelyscar-".$aatofind.".pml";
open(PYMOL,">".$outpymol);
#
## OUTPUT FILE: AMINO ACID VECTOR
my $outeuclideanoutvector = $pdb.".chain".$chain.".euclidean".$distance."-prelyscar-".$aatofind.".csv";
open (VECTOR,">".$outeuclideanoutvector);
#  p r e l y s c a r oooooooooooooooooooooooooooooooout

# VECTOR LABELS:
# AAS AA# PDB NAA ASP GLU NEG LYS ARG HIS POS PyH NON CHA ST NQ POL OTH ARO HYD k EC1 RES
#print VECTOR $aatofind." AA# PDB NAA ASP GLU NEG LYS ARG HIS POS PyH NON CHA ST NQ POL OTH ARO HYD k EC1 RES BFA\n";

my %pdbtopredict = ();
my @lysres = ();
my $s = 0;

##General stuff
print PYMOL "bg_color white\n";
print PYMOL "hide all\n";
print PYMOL "show cartoon, chain $chain\n";
print PYMOL "set cartoon_transparency, 0.7\n";
print PYMOL "color gray85, chain $chain\n";
print PYMOL "select heteroatoms, (hetatm and not resn HOH) AND chain $chain\n";
print PYMOL "show sticks, heteroatoms\n";
print PYMOL "color magenta, heteroatoms\n";

for my $thechosenaa (sort {$a<=>$b} keys %listofnearaminoacidstothechoseone)
{
#	print "Seq around $aatofind ".$thechosenaa." ";
	##PRINTING OUT RESIDUES AROUND THE MODLYS
	print PYMOL "select $aatofind$thechosenaa, resi $thechosenaa AND chain $chain\n";
	print PYMOL "select around$aatofind$thechosenaa, resi ";
	my $number = keys %{$listofnearaminoacidstothechoseone{$thechosenaa}};
#	print "(we have $number):\n";
	my $sumbfactor = 0;
	my $sequencearound = '';
	for my $theotheraanumber (sort {$a<=>$b} keys %{$listofnearaminoacidstothechoseone{$thechosenaa}})
	{
		print PYMOL $theotheraanumber."+";
		my $aa3 = $listofaasandnames{$theotheraanumber};
		my $aa1 = iubjust3to1($aa3);
		if ($aa1)
		{
			$sequencearound .= $aa1;
		}
#		print $theotheraanumber."($aa3 - $aa1) ".$listofnearaminoacidstothechoseone{$thechosenaa}{$theotheraanumber}." $sequencearound\n";
		$sumbfactor += $listofnearaminoacidstothechoseone{$thechosenaa}{$theotheraanumber};
	}
	
	my($averagebfactor) =  sprintf("%.4f",($sumbfactor/$number));
#	print "average b-factor ".$averagebfactor."\n";
	my ($R,$D,$E,$H,$K,$pNoncharged,$pCharged,$ST,$NQ,$pPolar,$pOthers,$pAromatics,$pHydrophobics, $pk,$i,$o) = AminoAcidFrequenciesBayesian($sequencearound);
	my $positives = $R+$K;
	my $posandhis = $R+$K+$H;
	my $negatives = $D+$E;
	
	my $seqNOwaterandions = $sequencearound;
	$seqNOwaterandions =~ s/o//g;
	$seqNOwaterandions =~ s/i//g;
	my $seqlen = length($sequencearound);
	my $seqlennowater = length($seqNOwaterandions);
	
	my $aveeuclideanprint = sprintf("%.4f",($averageeuclideans{$thechosenaa}));
	my $pdbchain = $PDB.".chain".$chain;
	# LABELS
	# AA NUM PDB LEN D E NEG K R H POS POSHIS NONCH CHA ST NQ POL OTH ARO HYD KCX EC RES BFA AVEEUCLIDEAN ION WAT PROTLENCODE PROTLEN
	print VECTOR $aatofind." $thechosenaa $pdbchain $seqlennowater $D $E $negatives $K $R $H $positives $posandhis $pNoncharged $pCharged $ST $NQ $pPolar $pOthers $pAromatics $pHydrophobics $i $o $protlencode $protlen\n";
	my @tempall = ();
	push (@tempall, ($aatofind,$thechosenaa,$pdbchain,$seqlennowater,$D,$E,$negatives,$K,$R,$H,$positives,$posandhis,$pNoncharged,$pCharged,$ST,$NQ,$pPolar,$pOthers,$pAromatics,$pHydrophobics,$i,$o,$protlencode,$protlen));
	push (@lysres, $thechosenaa);

	# ppppppppppppppppp
#	print "@tempall ";
#	print "\n";
	# ppppppppppppppppp

	$pdbtopredict{$s} = [@tempall];
	$s++;
	print PYMOL " AND chain $chain\n";
#	print $sequencearound." AveEuclidean:".$aveeuclideanprint."\n";
} #for my $thechosenaa
	
print PYMOL "show sticks, $aatofind* AND chain $chain\n";
print PYMOL "color red, $aatofind* AND chain $chain\n";
print PYMOL "color tv_blue, around* AND chain $chain\n";

print PYMOL "select histidines, resn his AND chain $chain\n";
print PYMOL "show sticks, histidines\n";
print PYMOL "show spheres, heteroatoms\n";
print PYMOL "deselect\n";

close PYMOL;
close VECTOR;

print "\nLysine residues found in this protein chain: \n";
my $d = 1;
foreach my $lys (@lysres)
{
	print "LYS".$lys." ";
	if ($d % 7 == 0)
	{
		print "\n";
	}
	$d++;
}

# BAYESIAN PREDICTION METHOD
# Adapted from: bayesianPreLysCar.pl > 
# Calculate Posterior probability of residues being carboxylated and non-carboxylated.
# ---------------------------------------------------------

# MAP of FEATURES 
# 3  4 5  6  7 8 9  10   11     12   13 14 15  16  17  18  19 20  21  22          23                             
#LEN D E NEG K R H POS POSHIS NONCH CHA ST NQ POL OTH ARO HYD ION WAT PROTLENCODE PROTLEN

# Using the Features giving the BEST performance in the leave-one-out cross validation test with the training data set.
my @chosenfeatures = (3,6,11,14,15,18,19,20,21); #BEST performance

# TAKE SELECTED FEATURES, transfer to hash
my %ilikehash = ();
foreach my $ele (@chosenfeatures)
{
	$ilikehash{$ele} = 1;
}

# LOADING TRAINING DATA SETS:
# -----------------------------------------------
# KCX set
my ($kcxhashofarrays_a) = loading_kcx_training();
my %kcxhashofarrays = %$kcxhashofarrays_a;
my $kcxtotalnumber = keys %kcxhashofarrays;

# LYS set
my ($lyshashofarrays_a) = loading_lys_training();
my %lyshashofarrays = %$lyshashofarrays_a;
my $lystotalnumber = keys %lyshashofarrays;


########################################################
#  FREQUENCIES FOR EACH AMINO ACID
########################################################

# THE FREQUENCIES FOR LYS RESIDUES
my %frequencieslys = Calculate_Probability_Event(\%lyshashofarrays,\%ilikehash);
#pppppppppppppppppppppppppppppppppppppppppppppppppppppppppp	
#	for my $bins (sort {$a<=>$b} keys %frequencieslys)
#	{
#		for my $feature (sort keys %{$frequencieslys{$bins}})
#		{
#			my $valor = $frequencieslys{$bins}{$feature};
#			print $bins." ".$feature." ".$valor."\n";
#		}
#		print "\n";
#	}
#pppppppppppppppppppppppppppppppppppppppppppppppppppppppppp	

#THE FREQUENCIES FOR KCX RESIDUES
my %frequencieskcx = Calculate_Probability_Event(\%kcxhashofarrays,\%ilikehash);
#pppppppppppppppppppppppppppppppppppppppppppppppppppppppppp	
#	for my $bins (sort {$a<=>$b} keys %frequencieskcx)
#	{
#		for my $feature (sort keys %{$frequencieskcx{$bins}})
#		{
#			my $valor = $frequencieskcx{$bins}{$feature};
#			print $bins." ".$feature." ".$valor."\n";
#		}
#		print "\n";
#	}
#pppppppppppppppppppppppppppppppppppppppppppppppppppppppppp	

# Control:
my $hits = 0;

# GRAVING THE VECTORS
my @allkeys= keys %pdbtopredict; 

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# LENCONTROL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Control: this number only can be 0, 1 or 2
my $sizecheck = '';
$sizecheck = $pdbtopredict{0}[22];
if ($sizecheck !~ /[0,1,2]/)
{
	print "FATAL ERROR PROBLEM!!!!!!!!!\n";
	print $command." error 139\n";
	print "Please, contact us at prelyscar\@gmail.com\n";

	exit;
}
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# LENCONTROL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

foreach my $therandomkey (@allkeys)
{
	# EXTRACTING ONE AT RANDOM
	my %onevectorrandom = ();
	my %alltheothervectors = ();
	

	my @yep = @{$pdbtopredict{$therandomkey} };
	$onevectorrandom{$therandomkey} = [@yep];
	
	my $aan = $onevectorrandom{$therandomkey}[1];
	my $pdb = $onevectorrandom{$therandomkey}[2];
	
	# Get the chain ID:
	my $chainid = '';
	
	if ($pdb =~ /.*\.(chain\w)/)
	{
		$chainid = $1;	
	}
	
	# I am doing this just because of the original lysloocv version
	my %kckminusone = %pdbtopredict;

	
	# STRUCTURE OF THE REDUCED RANDOM VECTOR:
	# $reducedrandomvector{FEATURE} = NUMBEROFTIMES (also call, the BINS of the prob matrices);
	my %reducedrandomvector = ();
	for my $features (sort {$a<=>$b} keys (%ilikehash))
	{
		$reducedrandomvector{$features} = $onevectorrandom{$therandomkey}[$features];
#		print $features." ".$onevectorrandom{$therandomkey}[$features]." ".$reducedrandomvector{$features}."\n";
	}
	
	my ($posteriorC1,$posteriorC2) = Calculate_Posterior_Probability($prior,\%reducedrandomvector,\%frequencieskcx,\%frequencieslys);

	if ($posteriorC1>$posteriorC2)
	{
		$hits++;

		# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		# LENCONTROL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if ($sizecheck == 1)
		{
			print "\n____________________________________________________________________\n";
			print "********************************************************************\n";
			print "Lysine $aan in $chainid predicted as carboxylated (probability=";
			printf ("%.5f ",$posteriorC1);
			print ")\n";
			print "... although this is a SMALL protein chain (less than 200 residues!)\n";
			print "therefore there is a higher probability of being a false positive\n";
			print "********************************************************************\n\n";
		}
		else
		{
			print "\n____________________________________________________________________\n";
			print "********************************************************************\n";
			print "Lysine $aan in $chainid predicted as carboxylated (probability=";
			printf ("%.5f ",$posteriorC1);
			print ")\n";
			print "********************************************************************\n\n";
		}
	}
}

if ($hits == 0)
{
	print "\n____________________________________________________________________\n";
	print "********************************************************************\n";
	print "\nNO carboxylated lysine residues predicted in this protein chain\n\n";
	print "____________________________________________________________________\n";
	print "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n\n";
	print "Please, cite: Jimenez-Morales D, et.al. (2014). Acta Cryst. D 70, 48-57.\n";
	print "For any question/suggestion/comment contact us at prelyscar\@gmail.com\n";
	print "\nThanks for using PreLysCar\n\n";
}
else
{
	print "Please, cite: Jimenez-Morales D, et.al. (2014). Acta Cryst. D 70, 48-57.\n";
	print "For any question/suggestion/comment contact us at prelyscar\@gmail.com\n";
	print "\nThanks for using PreLysCar\n\n";
}


exit;

###############################################################
# # ~ ~ s u b r o u t i n e s ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # #
###############################################################

sub get_file_data 
{
    my($filename) = @_;

    # Initialize variables
    my @filedata = ();
    
    unless( open(GET_FILE_DATA, $filename) ) 
    {
        print "Cannot open file \"$filename\"";
        exit;
    }
	@filedata = <GET_FILE_DATA>;
	close GET_FILE_DATA;
    return @filedata;
}

sub EuclideanDistance
{
	my ($x1,$y1,$z1,$x2,$y2,$z2) = @_;
	
	my $x = ($x1-$x2)**2;
	my $y = ($y1-$y2)**2;
	my $z = ($z1-$z2)**2;
	my $sum = $x+$y+$z;
	my $sqrtsum = sqrt($sum);
	
	return $sqrtsum;
	# my $euclidean = EuclideanDistance($x,$y,$z,$x2,$y2,$z2);
}

sub get_pdbsidechain_coordinateskcxcase
{
	my ($pdbfile,$aatofind,$subunit) = @_;
	
	my @openfile = get_file_data($pdbfile);
	
	my %aacoordinates = ();
	my %alltheother = ();
	# $COORDINTATES{AANUMBER}{ATOM} = [AA3][X][Y][X][BFACTOR]

	# Control
	my $initialchain = ''; #Just to make sure we only have one chain;

	foreach my $line (@openfile)
	{
		chomp($line);
		
		if ($line =~ /^ATOM/)
		{
			my $chainhere = substr($line,21,1);
			if ($chainhere ne $subunit)
			{
				next;
			}
			
			my $resnum = substr ($line,22,4);
			$resnum =~ s/\s//g;
			
			my $aa3 = substr($line,17,3);
			
			my $atomnumber = substr($line,6,5);
			$atomnumber =~ s/\s//g;
			
			my $atomname = substr($line,12,4); 
			$atomname =~ s/\s//g;
			
			my $bfactor = substr($line,60,5);
			$bfactor =~ s/\s//g;
			
			my $x = substr($line,30,8);
			$x =~ s/\s//g;
			my $y = substr($line,38,8);
			$y  =~ s/\s//g;
			my $z = substr($line,46,8);
			$z  =~ s/\s//g;
		
			my @tmp = ();	
			@tmp = ($resnum,$aa3,$x,$y,$z,$bfactor,$atomname);
			
			if ($aa3 eq $aatofind)
			{		
				# NOT TAKING atoms from the BACKBONE AND CARBOXYL GROUP FROM KCX
				unless ( ($atomname =~ /^N$/) | ($atomname =~ /^CA$/) | ($atomname =~ /^C$/) | ($atomname =~ /^O$/) | ($atomname =~ /^CX$/) | ($atomname =~ /^OQ1$/) | ($atomname =~ /^OQ2$/) ) 
				{
					$aacoordinates{$resnum}{$atomnumber} = [@tmp];
					# print "\tATOM:".$atomnumber." AA#:".$resnum." AA3:".$aa3." X:$x Y:$y Z:$z BFACTOR:$bfactor\n";
				}
				if ($atomname =~ /^NZ$/)
				{
					my ($xcxx,$ycxx,$zcxx) = GiveMeTheLysTip($pdbfile,$resnum);
					my $bfactortip = 0;
					my $atomcxname = "PCX"; #PROYECTED PCX
					my $atomnumber = 9999;
					my @temp = ($resnum,$aa3,$xcxx,$ycxx,$zcxx,$bfactortip,$atomcxname);
					$aacoordinates{$resnum}{$atomnumber} = [@tmp];
					#print "\t$resnum,$aa3,$xcxx,$ycxx,$zcxx... in!\n"
				}
			}
			# ALL THE OTHER AMINO ACIDS, TAKE THE WHOLE THING
			else
			{
				$alltheother{$resnum}{$atomnumber} = [@tmp];
			} #else
		} #if it's an ATOM
	}
	
	return (\%aacoordinates,\%alltheother);
	#	It is sent as reference, therefore write in the main script:
	#	my ($aacoor, $allothercoor) = get_pdbsidechain_coordinates($filename,$theaa);
	#	my %aacoordinates = %$aacoor;
	#	my %alltheothercoordinates = %$allothercoor;

} # End of sub get_pdbsidechain_coordinateskcxcase

# SAME AS get_pdbsidechain_coordinateskcxcase, BUT JUST FOR WATER AND IONS (HETATOM)
sub Get_PDBHETATM_Coordinates_Ion_Water
{
	my ($pdbfile,$subunit) = @_;
	
	my @openfile = get_file_data($pdbfile);
	
	# Loading the list of accepted ions

	my @oionfile = zmetalioncenters();
	my %listofions = ();
	foreach my $line (@oionfile)
	{
		chomp($line);
		if (!$listofions{$line})
		{
			$listofions{$line} = 1;
		}
	}
	
	my %ionhetatmcoordinates = ();
	my %watercoordinates = ();
	# $COORDINTATES{AANUMBER}{ATOM} = [AA3][X][Y][X][BFACTOR]

	# Control
	my $initialchain = ''; #Just to make sure we only have one chain;
	my $c = 0;

	foreach my $line (@openfile)
	{
		chomp($line);
		
		if ($line =~ /^HETATM/)
		{	

			my $chainhere = substr($line,21,1);
			if ($chainhere ne $subunit)
			{
				next;
			}

			my $resnum = substr ($line,22,4);
			$resnum =~ s/\s//g;
			
			my $aa3 = substr($line,17,3);
			$aa3 =~ s/\s//g;
	
			my $atomnumber = substr($line,6,5);
			$atomnumber =~ s/\s//g;
			
			my $atomname = substr($line,12,4); 
			$atomname =~ s/\s//g;
			
			my $bfactor = substr($line,60,5);
			$bfactor =~ s/\s//g;
			
			my $x = substr($line,30,8);
			$x =~ s/\s//g;
			my $y = substr($line,38,8);
			$y  =~ s/\s//g;
			my $z = substr($line,46,8);
			$z  =~ s/\s//g;
		
			my @tmp = ();	
			@tmp = ($resnum,$aa3,$x,$y,$z,$bfactor,$atomname);

			if ($listofions{$aa3})
			{		
				$ionhetatmcoordinates{$resnum}{$atomnumber} = [@tmp];
			#	print "\tHETATM:".$atomnumber." AA#:".$resnum." AA3:".$aa3." X:$x Y:$y Z:$z BFACTOR:$bfactor\n";
			}
			
			if ($aa3 =~/HOH/)
			{
				$watercoordinates{$resnum}{$atomnumber} = [@tmp];
			#	print "\tHETATM:".$atomnumber." AA#:".$resnum." AA3:".$aa3." X:$x Y:$y Z:$z BFACTOR:$bfactor\n";
			}
		} #if it's an ATOM
	}
	
	return (\%ionhetatmcoordinates,\%watercoordinates);
#	It has to be sent as reference: write in the main script:
#	my ($ionhetatmcoor, $watercoor) = Get_PDBHETATM_Coordinates_Ion_Water($pdbfile);
#	my %ionhetatmcoor = %$ionhetatmcoor;
#	my %watercoordinates = %$watercoor;

} #END OF sub Get_PDBHETATM_Coordinates_Ion_Water

# PROYECTING A "CARBOXYL GROUP" ON TOP OF EVERY LYSINE RESIDUE
# STEP 1
sub GiveMeTheLysTip
{
	my ($pdbfile,$aanumber) = @_;

	my @oinput = get_file_data($pdbfile);
	
	my ($xca,$xcb,$xcg,$xcd,$xce,$xnz,$xcx) = '';
	my ($yca,$ycb,$ycg,$ycd,$yce,$ynz,$ycx) = '';
	my ($zca,$zcb,$zcg,$zcd,$zce,$znz,$zcx) = '';
	
	foreach my $linea (@oinput)
	{
		chomp($linea);
		
		if ($linea =~ /^ATOM/)
		{
			my $aa3 = substr($linea,17,3);
			my $aan = substr($linea,22,4);
			$aan =~ s/\s//g;
			
			my $atomserialnumber = substr($linea,6,5);
			$atomserialnumber =~ s/\s//g;
			
			my $atomname = substr($linea,12,4); 
			$atomname =~ s/\s//g;
		
			my $x = substr($linea,30,8);
			$x =~ s/\s//g;
			my $y = substr($linea,38,8);
			$y  =~ s/\s//g;
			my $z = substr($linea,46,8);
			$z  =~ s/\s//g;
			
			if ( ($linea =~ /^ATOM/) & ($aan eq $aanumber) & ($atomname eq "CA"))
			{
		#		print $x." ".$y." ".$z."\n";
				$xca = $x;
				$yca = $y;
				$zca = $z;
			}
			elsif ( ($linea =~ /^ATOM/) & ($aan eq $aanumber) & ($atomname eq "CB"))
			{
		#		print $x." ".$y." ".$z."\n";
				$xcb = $x;
				$ycb = $y;
				$zcb = $z;
			}
			elsif ( ($linea =~ /^ATOM/) & ($aan eq $aanumber) & ($atomname eq "CG"))
			{
		#		print $x." ".$y." ".$z."\n";
				$xcg = $x;
				$ycg = $y;
				$zcg = $z;
			}
			elsif ( ($linea =~ /^ATOM/) & ($aan eq $aanumber) & ($atomname eq "CD"))
			{
		#		print $x." ".$y." ".$z."\n";
				$xcd = $x;
				$ycd = $y;
				$zcd = $z;
			}
			elsif ( ($linea =~ /^ATOM/) & ($aan eq $aanumber) & ($atomname eq "CE"))
			{
		#		print $x." ".$y." ".$z."\n";
				$xce = $x;
				$yce = $y;
				$zce = $z;
			}
			elsif ( ($linea =~ /^ATOM/) & ($aan eq $aanumber) & ($atomname eq "NZ"))
			{
		#		print $x." ".$y." ".$z."\n";
				$xnz = $x;
				$ynz = $y;
				$znz = $z;
			}
		}# if ($line)
	} #foreach my $linea (@oinput)
	
	#Trying to get here an average direction vector
	my ($dx1,$dy1,$dz1) = VectorDirection($xca,$yca,$zca,$xcb,$ycb,$zcb);
	my ($dx2,$dy2,$dz2) = VectorDirection($xcb,$ycb,$zcb,$xcg,$ycg,$zcg);
	my ($dx3,$dy3,$dz3) = VectorDirection($xcg,$ycg,$zcg,$xcd,$ycd,$zcd);
	my ($dx4,$dy4,$dz4) = VectorDirection($xcd,$ycd,$zcd,$xce,$yce,$zce);
	my ($dx5,$dy5,$dz5) = VectorDirection($xce,$yce,$zce,$xnz,$ynz,$znz);
	
	my $avedx = ($dx1+$dx2+$dx3+$dx4+$dx5)/5;
	my $avedy = ($dy1+$dy2+$dy3+$dy4+$dy5)/5;
	my $avedz = ($dz1+$dz2+$dz3+$dz4+$dz5)/5;
	
	# GET DIRECTION OF A HYPOTHETICAL TIP ON NZ
	my $xcxx = $xnz+($avedx);
	my $ycxx = $ynz+($avedy);
	my $zcxx = $znz+($avedz);

	return ($xcxx,$ycxx,$zcxx);
#	my ($xcxx,$ycxx,$zcxx) = GiveMeTheLysTip($pdbfile,$aanumber);
} # End of GiveMeTheLysTip

# PROYECTING A "CARBOXYL GROUP" ON TOP OF EVERY LYSINE RESIDUE
# STEP 2
sub VectorDirection 
{
	my($x1,$y1,$z1,$x2,$y2,$z2) = @_;
	
	my $dx = $x2-$x1;
	my $dy = $y2-$y1;
	my $dz = $z2-$z1;
	
	return ($dx,$dy,$dz);
	#Main body:	my ($dx1,$dy1,$dz1) = VectorDirection($x1,$y1,$z1,$x2,$y2,$z2);
}

sub iubjust3to1 {

    my($input) = @_;
    
    my %three2one = (
		'ALA' => 'A',
		'VAL' => 'V',
		'LEU' => 'L',
		'ILE' => 'I',
		'PRO' => 'P',
		'TRP' => 'W',
		'PHE' => 'F',
		'MET' => 'M',
		'GLY' => 'G',
		'SER' => 'S',
		'THR' => 'T',
		'TYR' => 'Y',
		'CYS' => 'C',
		'ASN' => 'N',
		'GLN' => 'Q',
		'LYS' => 'K',
		'ARG' => 'R',
		'HIS' => 'H',
		'ASP' => 'D',
		'GLU' => 'E',
		'ASX' => 'B',
		'GLX' => 'Z',
		'XLE' => 'J',
		'XAA' => 'X',
		'UNK' => 'X',
		'MSE' => 'm',
		# MODIFIED RESIDUES
		'HYP' => 'p',	# PRO  4-HYDROXYPROLINE                                                       
		'KCX' => 'k',	# LYS  LYSINE NZ-CARBOXYLIC ACID                                            
		'SMC' => 'c',	# CYS  S-METHYLCYSTEINE                                   
		'MME' => 'm',	# MET  N-METHYL METHIONINE            
		'PTR' => 'y',	# TYR  O-PHOSPHOTYROSINE       
#		'LLP' => 'k',	# 2-LYSINE(3-HYDROXY-2-METHYL-5-PHOSPHONOOXYMETHYL
		'OCS' => 'c',	# CYS  CYSTEINE SULFONIC ACID      
#		'ALY' => 'k',	# LYS  MODIFIED LYSINE
		'LLP' => 'x',   # LLP 2-LYSINE(3-HYDROXY-2-METHYL-5-PHOSPHONOOXYMETHYL-PYRIDIN-4-YLMETHANE) N'-PYRIDOXYL-LYSINE-5'-MONOPHOSPHATE LLP    4(C14 H24 N3 O7 P)
		'CME' => 'c',	# S,S-(2-HYDROXYETHYL)THIOCYSTEINE
      

		# WATER MOLECULES
		'HOH' => 'o',	# WATER MOLECULE
		# IONS
		'MG'  => 'i', 	# MAGNESIUM ION
		'SO4' => 'i',	# SULFATE ION
		'CL'  => 'i',	# CHLORIDE ION
		'NI'  => 'i',	# NICKEL(II) ION
		'MN'  => 'i',	# MANGANESE(II) ION
		'ZN'  => 'i',	# ZINC ION
		'CO3' => 'i',	# CARBONATE ION
		'PO4' => 'i',	# PHOSPHATE ION
		'CD'  => 'i',	# CADMIUM ION
		'NA'  => 'i',	# SODIUM ION
		'CA'  => 'i',	# CALCIUM ION
		'CO'  => 'i',	# COBALT(II) ION
		'MLI' => 'i',	# MALONATE ION
		'CAC' => 'i',	# CACODYLATE ION
		'ACT' => 'i',	# ACETATE ION
		'FE2' => 'i',	# FE (II) ION
		'FE'  => 'i',	# FE(III) ION
		'AZI' => 'i',	# AZIDE ION
		'SO3' => 'i',	# SULFITE ION
		'BR'  => 'i',	# BROMIDE ION 		
    ); #END OF THE ARRAY
    
    my $seq = '';
	if(not defined $three2one{$input}) 
	{
		print "\t[SUB iubjust3to1] Code $input not defined\n";
	}
	else
	{
		$seq = $three2one{$input};
	}
    return $seq;
}

# FREQUENCY OF EACH RESIDUE (INCLUDING DOUBLES, UNKNOWN, AND MODIFIED RESIDUES)
sub AminoAcidFrequenciesBayesian
{
	my($protein) = @_;
	
	#---------------------------------------------------------------------------
	my($A) = ($protein =~ tr/A//);	# Ala  	  Alanine
	my($R) = ($protein =~ tr/R//);	# Arg 	  Arginine
	my($N) = ($protein =~ tr/N//);	# Asn 	  Asparagine
	my($D) = ($protein =~ tr/D//);	# Asp 	  Aspartic acid
	my($C) = ($protein =~ tr/C//);	# Cys 	  Cysteine
	my($Q) = ($protein =~ tr/Q//);	# Gln 	  Glutamine
	my($E) = ($protein =~ tr/E//);	# Glu 	  Glutamic acid
	my($G) = ($protein =~ tr/G//);	# Gly 	  Glycine
	my($H) = ($protein =~ tr/H//);	# His 	  Histidine
	my($I) = ($protein =~ tr/I//);	# Ile 	  Isoleucine
	my($L) = ($protein =~ tr/L//);	# Leu 	  Leucine
	my($K) = ($protein =~ tr/K//);	# Lys 	  Lysine
	my($M) = ($protein =~ tr/M//);	# Met 	  Methionine
	my($F) = ($protein =~ tr/F//);	# Phe 	  Phenylalanine
	my($P) = ($protein =~ tr/P//);	# Pro 	  Proline
	my($S) = ($protein =~ tr/S//);	# Ser 	  Serine
	my($T) = ($protein =~ tr/T//);	# Thr 	  Threonine
	my($W) = ($protein =~ tr/W//);	# Trp 	  Tryptophan
	my($Y) = ($protein =~ tr/Y//);	# Tyr 	  Tyrosine
	my($V) = ($protein =~ tr/V//);	# Val 	  Valine
	my($B) = ($protein =~ tr/B//);	# Asx 	  Aspartic acid or Asparagine
	my($Z) = ($protein =~ tr/Z//);	# Glx 	  Glutamic acid or Glutamine
	my($X) = ($protein =~ tr/X//);	# Xaa 	  Any amino acid
	my($d) = ($protein =~ tr/-//);  # Dash  		NUMBER OF DASH (-)
	my($m) = ($protein =~ tr/m//);	# 'MSE' => 'm' || 'MME' => 'm',	# MET  N-METHYL METHIONINE     
	my($p) = ($protein =~ tr/p//);	# 'HYP' => 'p',	# PRO  4-HYDROXYPROLINE		
	my($k) = ($protein =~ tr/k//);	# 'KCX' => 'k',	# LYS  LYSINE NZ-CARBOXYLIC ACID | 'LLP' => 'k',	# 2-LYSINE(3-HYDROXY-2-METHYL-5-PHOSPHONOOXYMETHYL |'ALY' => 'k',	# LYS  MODIFIED LYSINE         
	my($c) = ($protein =~ tr/c//);	# 'SMC' => 'c',	# CYS  S-METHYLCYSTEINE | 		'OCS' => 'c',	# CYS  CYSTEINE SULFONIC ACID  
	my($y) = ($protein =~ tr/y//);	# 'PTR' => 'y',	# TYR  O-PHOSPHOTYROSINE
	my($i) = ($protein =~ tr/i//);  # 'ION' => 'i',	# ION. IT CAN BE ANYTHING. CHECK iubjust3to1
	my($o) = ($protein =~ tr/o//);  # 'HOH' => 'w',	# Wather Molecule
	#---------------------------------------------------------------------------
	
	my ($hydrophobics) = $M+$I+$L+$V;
	my ($others) = $C+$G+$P+$A;
	my ($aromatics) = $F+$W+$Y;
	my ($polar) = $S+$T+$Q+$N;
	my ($ST) = $S+$T;
	my ($QN) = $Q+$N;

	my ($charged) = $R+$D+$E+$H+$K+$B+$Z;
	my ($noncharged) = $A+$I+$L+$V+$M+$C+$G+$P+$F+$W+$Y+$S+$T+$Q+$N+$k;

	my ($total) = $A+$I+$L+$V+$M+$C+$G+$P+$F+$W+$Y+$S+$T+$Q+$N+$R+$D+$E+$H+$K+$B+$Z+$k;
	
	return ($R,$D,$E,$H,$K,$noncharged,$charged,$ST,$QN,$polar,$others,$aromatics,$hydrophobics,$k,$i,$o);
}


#################################################
# CALCULTATE PROBABILITY OF AN EVENT
sub Calculate_Probability_Event
{
	my ($allaaafeatures_h,$hashoffeatures_h) = @_;
	
	my %allaaafeatures = %$allaaafeatures_h;
	my %hashoffeatures = %$hashoffeatures_h;

	##pppppppppppppppppppppppppppppppppppppppppppppppppppppppppp	
#	while ( my($joder,$suputa) = each (%hashoffeatures) )
#	{
#		print $joder." and ".$suputa."\n";
#	}
#	
#	for $family (sort {$a<=>$b} keys %allaaafeatures ) 
#	{
#		print "$family: @{ $allaaafeatures{$family} }\n";
#	}
	##pppppppppppppppppppppppppppppppppppppppppppppppppppppppppp

	########################################################
	#FREQUENCIES FOR EACH AMINO ACID
	# $aaafeaturesfrequencies{PREINDEX FROM 0 to 15?}{feature1} = number of observations for this value for this feature
	my %aaafeaturesfrequencies = ();
	
	# INITIALIZE THE HASH %aaafeaturesfrequencies
	# I need to create the frequencies indexes (I choose 15 by now):
	# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	for (my $i = 0;$i<=30;$i++)
	{
		for my $v (sort {$a<=>$b} keys %hashoffeatures) 
		{
#			print $i." ".$v." = 0\n";
			$aaafeaturesfrequencies{$i}{$v} = 0;
		}
	}
	
	##pppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
#	for my $joder (sort {$a<=>$b} keys %aaafeaturesfrequencies)
#	{
#		print $joder."\n";
#		for my $laputa (keys %{$aaafeaturesfrequencies{$joder}})
#		{
#			print "\t".$laputa." ".$aaafeaturesfrequencies{$joder}{$laputa}."\n";
#		}
#	}
	##pppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
		
	#FILLING UP THE HASH FOR FREQUENCIES:
	for my $k (sort {$a<=>$b} keys %allaaafeatures ) 
	{
		for my $v (sort {$a<=>$b} keys %hashoffeatures) 
		{
#			print "Vector:".$k." ".$v.": has the value ";
			my $featurevalue = $allaaafeatures{$k}[$v];
#			print $featurevalue."\n";
			# SUMMING UP THIS ELEMENT FOR THIS FEATURE
			if (!$aaafeaturesfrequencies{$featurevalue}{$v})
			{
				$aaafeaturesfrequencies{$featurevalue}{$v} = 1;
#				print "\t".$featurevalue." ".$v." ".$aaafeaturesfrequencies{$featurevalue}{$v}."\n";
			}
			else
			{
				$aaafeaturesfrequencies{$featurevalue}{$v}++;
#				print "\t".$featurevalue." ".$v." ".$aaafeaturesfrequencies{$featurevalue}{$v}."\n";
			}
		}
#	    print "\n";
	}

	##ppppppppppppppppppepppppppppppppppppppppppppppppppppppppppp
	# PRINTING OUT FREQUENCIES
#	for my $joder (sort {$a<=>$b} keys %aaafeaturesfrequencies)
#	{
#		print $joder."\n";
#		for my $laputa (keys %{$aaafeaturesfrequencies{$joder}})
#		{
#			print "\t".$laputa." ".$aaafeaturesfrequencies{$joder}{$laputa}."\n";
#		}
#	}
	##pppppppppppppppppppppppppppppppppppppppppppppppppppppppppp

	# ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
	# TAKE THE SUM PER COLUM, and at the same time, CHECK POINT
	# It's a check point because I make sure that all the columns 
	# has exactly the same value, all of them, because it does not make
	# sense otherwise. This is because in my training data I don't miss
	# even one single value, therefore all the feature must sum up to 
	# the exact same value for all of them
	my $totalcheckpoint = 0;
	my $ctcp = 1;
	for my $venga (keys %hashoffeatures)
	{
		my $totalhere = 0;
		for my $vamos (sort {$a<=>$b} keys %aaafeaturesfrequencies)
		{
#			print $venga." ".$vamos." ".$aaafeaturesfrequencies{$vamos}{$venga}."\n";
			$totalhere += $aaafeaturesfrequencies{$vamos}{$venga};
		}
			if ($totalhere == 0)
			{
				print "Oh no!\n";
				exit;
			}
#		print "--------\n";
#		print "total = $totalhere\n";
		if ($ctcp == 1)
		{
			$totalcheckpoint = $totalhere;
		}
		else
		{
			if ($totalcheckpoint != $totalhere)
			{ 
				print "<br/>FATAL ERROR PROBLEM<br/>\n";
				print $command." error 385<br/>\n";
				print "<small>Please, contact us at <a href='mailto:prelyscar\@gmail.com?Subject=Fatal error $command 385' target='_top'>prelyscar\@gmail.com</a></small> <br/>\n";
				exit;
			}
		}
		$ctcp++;
	}
	
	# FREQUENCY MATRIX (Fx|C)
	my %frequenciesaaa = ();
	
	for my $bins (sort {$a<=>$b} keys %aaafeaturesfrequencies)
	{
	#	print $oh." = ";
		for my $feature (sort keys %{$aaafeaturesfrequencies{$bins}})
		{
			my $valor = $aaafeaturesfrequencies{$bins}{$feature};
#			print $bins." ".$feature." ".$valor." ";
#			print "\nTOTALCHECKPOINT:".$totalcheckpoint." ";
			my ($probability) =  sprintf("%.4f",(($valor/$totalcheckpoint)));
		
			if ($probability == 0)
			{
				$frequenciesaaa{$bins}{$feature} = 0.0001;
			}
			else
			{
				$frequenciesaaa{$bins}{$feature} = $probability;
#				print $probability."\n";
			}
		}
#		print "\n";
	}
	
	##pppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
	# PRINTING OUT PROBABILITIES
#	for my $bins (sort {$a<=>$b} keys %frequenciesaaa)
#	{
#		for my $feature (sort keys %{$frequenciesaaa{$bins}})
#		{
#			my $valor = $frequenciesaaa{$bins}{$feature};
#			print $bins." ".$feature." ".$valor."\n";
#		}
#		print "\n";
#	}
	##pppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
#	exit;
	return(%frequenciesaaa);
#	my %frequenciesaaa = Calculate_Probability_Event(%aaafeaturesfrequencies);
} #Calculate_Probability_Event



#################################################
# CALCULTATE PROBABILITY OF AN EVENT
sub Calculate_Posterior_Probability
{
	my ($prior,$vector_h,$prob_matrix_a_h,$prob_matrix_b_h) = @_;

	my $otherprior = (1-$prior);
	my %vector = %$vector_h;
	my %prob_matrix_a = %$prob_matrix_a_h;
	my %prob_matrix_b = %$prob_matrix_b_h;

#pppppppppppppppppppppppppppppppppppppppppppppppppppppppppp	
#	for my $bins (sort {$a<=>$b} keys %prob_matrix_a)
#	{
#		for my $feature (sort keys %{$prob_matrix_a{$bins}})
#		{
#			my $valor = $prob_matrix_a{$bins}{$feature};
#			print $bins." ".$feature." ".$valor."\n";
#		}
#		print "\n";
#	}
#pppppppppppppppppppppppppppppppppppppppppppppppppppppppppp	
		
	my $probabilityFC1 = 1;
	my $probabilityFC2 = 1;
	for my $f (sort {$a<=>$b} keys (%vector) )
	{
		my $bins = $vector{$f};
#		print $bins." ".$f." ";
		
		#Now KCX
		my $P_f1_if_C1 = $prob_matrix_a{$bins}{$f};
#		print " = ".$P_f1_if_C1." and ";
		$probabilityFC1 *= $P_f1_if_C1;

		#Now LYS
		my $P_f1_if_C2 = $prob_matrix_b{$bins}{$f};
#		print $P_f1_if_C2."\n";
		$probabilityFC2 *= $P_f1_if_C2;
	}
#	print "Likelihood KCX ".$probabilityFC1."\n";
#	print "Likelihood LYS ".$probabilityFC2."\n";
	
	my $priorMlikelihoodC1 = $prior*$probabilityFC1;
	my $priorMlikelihoodC2 = $otherprior*$probabilityFC2;
	my $evidence = $priorMlikelihoodC1+$priorMlikelihoodC2;
	my $posteriorC1 = ($priorMlikelihoodC1/$evidence);
	my $posteriorC2 = ($priorMlikelihoodC2/$evidence);
	
#	print "Now posteriors: C1: $posteriorC1 and $posteriorC2\n";
	return ($posteriorC1,$posteriorC2);
}


sub open_file_in_hash_of_arrays
{
	my($filename) = @_;
	my @ofilename = get_file_data($filename);
	
	my $sizecontrol = 0;
	my $c = 0;
	
	my %hashofarrays = ();
	
	foreach my $line (@ofilename)
	{
		chomp($line);
#		print $line."\n";
		my @temp = split(" ",$line);
		my $size = @temp;
		#Size control for the first one
		if ($c == 0)
		{
			$sizecontrol = $size;
		}
		
		#All the others - -- - - - -- ---
		if($size != $sizecontrol)
		{
			print "\n\nFATAL ERROR PROBLEM\n";
			print $command." error 518\n";
			print "Please, contact us at prelyscar\@gmail.com\n\n";
			exit;
		}
		#- - - - - - -- - -end of control
		
		$hashofarrays{$c} = [@temp];
		
		$c++;
	}
	#ppppppppppppppppppppppppppppppppppppppppppppppppp
#	for my $k (sort {$a<=>$b} keys %hashofarrays ) 
#	{
#		print "$k=>";
#		for my $v ( 0 .. $#{ $hashofarrays{$k} } ) 
#		{
#			print " $v:$hashofarrays{$k}[$v] ";
#		}
#	    print "\n";
#	}
	#ppppppppppppppppppppppppppppppppppppppppppppppppp
	return (\%hashofarrays);
} #open_file_in_hash_of_arrays


# LOADING DATA----------------------------------------------
# List of metal ions known involved in lysine carboxylation
sub zmetalioncenters 
{
	my @listmetal = ("ZN","MG","CO","FE2","FE","NI","MN","CD");
	return(@listmetal)
}

# Loading DATA
# KCX training data
sub loading_kcx_training
{
	my $kcxtraining = "KCX 201 1AA1.subB 12 2 1 3 0 0 3 0 3 6 6 2 1 3 1 1 1 1 2 2 438
KCX 201 1BWV.subA 12 2 1 3 0 0 3 0 3 6 6 1 1 2 1 1 2 1 2 2 471
KCX 70 1E4D.subA 12 0 0 0 1 0 0 1 1 11 1 3 0 3 2 3 3 0 2 2 243
KCX 224 1E8C.subA 16 0 1 1 1 0 0 1 1 14 2 4 1 5 2 4 3 0 1 2 496
KCX 219 1E9Y.subB 14 0 0 0 0 0 4 0 4 10 4 2 0 2 2 2 4 2 0 2 569
KCX 219 1E9Z.subB 13 0 0 0 0 0 3 0 3 10 3 2 0 2 2 2 4 2 0 2 569
KCX 1217 1EF2.subA 13 0 0 0 0 0 3 0 3 10 3 2 0 2 3 1 4 2 0 2 566
KCX 1217 1EJR.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 3 1 4 2 1 2 553
KCX 1217 1EJS.subC 13 0 0 0 0 0 2 0 2 11 2 2 1 3 3 1 4 2 1 2 566
KCX 1217 1EJT.subC 13 0 0 0 0 0 2 0 2 11 2 2 1 3 3 1 4 2 1 2 566
KCX 1217 1EJU.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 3 1 4 2 1 2 553
KCX 1217 1EJV.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 3 1 4 2 2 2 553
KCX 1217 1EJW.subC 14 0 0 0 0 0 4 0 4 10 4 2 0 2 3 1 4 2 2 2 566
KCX 1217 1EJX.subC 14 0 0 0 0 0 4 0 4 10 4 2 0 2 3 1 4 2 1 2 556
KCX 129 1EPV.subA 11 0 0 0 0 0 2 0 2 9 2 2 0 2 1 2 4 0 3 2 382
KCX 169 1EYW.subA 13 0 0 0 0 0 3 0 3 10 3 4 0 4 1 1 4 2 1 2 330
KCX 169 1EZ2.subA 15 0 0 0 0 0 3 0 3 12 3 4 0 4 2 1 5 2 0 2 329
KCX 129 1FTX.subA 11 0 0 0 0 0 2 0 2 9 2 2 0 2 1 2 4 0 3 2 380
KCX 217 1FWA.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 3 1 4 2 2 2 566
KCX 217 1FWB.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 3 1 4 2 2 2 566
KCX 217 1FWC.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 3 1 4 2 2 2 566
KCX 217 1FWD.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 3 1 4 2 2 2 566
KCX 217 1FWE.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 3 1 4 2 0 2 552
KCX 217 1FWF.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 3 1 4 2 1 2 551
KCX 217 1FWG.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 3 1 4 2 2 2 566
KCX 217 1FWH.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 3 1 4 2 1 2 566
KCX 217 1FWI.subC 12 0 0 0 0 0 2 0 2 10 2 2 0 2 3 1 4 1 2 2 551
KCX 217 1FWJ.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 3 1 4 2 2 2 566
KCX 201 1GK8.subA 13 2 1 3 0 0 3 0 3 7 6 2 1 3 1 1 2 1 1 2 469
KCX 150 1GKP.subA 12 0 0 0 0 0 3 0 3 9 3 1 0 1 2 2 4 2 1 2 458
KCX 150 1GKQ.subA 12 0 0 0 0 0 3 0 3 9 3 1 0 1 2 2 4 2 0 2 458
KCX 147 1GKR.subA 15 0 0 0 0 0 3 0 3 12 3 1 0 1 2 3 6 2 0 2 451
KCX 33 1H01.subA 13 1 0 1 1 0 0 1 1 11 2 1 0 1 2 2 6 0 1 2 287
KCX 338 1HL8.subA 7 1 2 3 0 0 0 0 0 4 3 0 1 1 1 0 2 0 0 2 426
KCX 338 1HL9.subA 8 1 2 3 0 0 0 0 0 5 3 0 1 1 1 0 3 0 0 2 421
KCX 220 1IE7.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 2 1 5 2 0 2 570
KCX 201 1IR1.subA 13 2 1 3 0 0 3 0 3 7 6 2 1 3 1 1 2 1 1 2 464
KCX 201 1IR2.subA 13 2 1 3 0 0 4 0 4 6 7 2 1 3 1 1 1 1 1 2 468
KCX 102 1J79.subA 14 0 1 1 0 0 3 0 3 10 4 1 0 1 3 2 4 2 0 2 343
KCX 185 1JBV.subA 13 0 1 1 0 0 1 0 1 11 2 2 0 2 4 0 5 0 3 2 404
KCX 185 1JBW.subA 13 0 1 1 0 0 1 0 1 11 2 2 0 2 5 0 4 0 2 2 410
KCX 169 1JGM.subA 13 0 0 0 0 0 3 0 3 10 3 4 0 4 1 1 4 2 2 2 333
KCX 150 1K1D.subA 13 0 0 0 0 0 4 0 4 9 4 0 0 0 0 2 7 2 0 2 460
KCX 70 1K38.subA 13 0 0 0 1 0 0 1 1 12 1 5 0 5 1 3 3 0 1 2 239
KCX 70 1K4E.subA 12 0 0 0 1 0 0 1 1 11 1 3 0 3 2 3 3 0 1 2 245
KCX 70 1K4F.subA 12 0 0 0 1 0 0 1 1 11 1 3 0 3 2 3 3 0 1 2 244
KCX 70 1K54.subA 12 0 0 0 1 0 0 1 1 11 1 3 0 3 2 3 3 0 1 2 244
KCX 70 1K55.subA 12 0 0 0 1 0 0 1 1 11 1 3 0 3 2 3 3 0 1 2 245
KCX 70 1K56.subA 12 0 0 0 1 0 0 1 1 11 1 3 0 3 2 3 3 0 1 2 243
KCX 70 1K57.subA 11 0 0 0 1 0 0 1 1 10 1 3 0 3 2 2 3 0 2 2 243
KCX 70 1K6S.subA 12 0 0 0 1 0 0 1 1 11 1 3 0 3 2 3 3 0 1 2 244
KCX 217 1KRB.subC 13 0 0 0 0 0 2 0 2 11 2 2 0 2 4 1 4 2 0 2 566
KCX 129 1L6F.subA 12 0 0 0 0 1 2 1 3 9 3 2 0 2 1 2 4 0 3 2 382
KCX 129 1L6G.subA 12 0 0 0 0 1 2 1 3 9 3 2 0 2 1 2 4 0 3 2 382
KCX 70 1M6K.subA 13 0 0 0 1 0 0 1 1 12 1 4 0 4 0 3 5 0 0 2 250
KCX 148 1NFG.subA 12 0 0 0 0 0 3 0 3 9 3 0 0 0 0 3 6 2 0 2 457
KCX 129 1NIU.subA 12 0 0 0 0 1 2 1 3 9 3 2 0 2 1 2 4 0 2 2 382
KCX 33 1OIR.subA 13 1 0 1 1 0 0 1 1 11 2 1 0 1 2 2 6 0 1 2 287
KCX 162 1ONW.subA 14 0 0 0 0 0 3 0 3 11 3 1 0 1 4 1 5 2 2 2 376
KCX 162 1ONX.subA 14 0 0 0 0 0 3 0 3 11 3 1 0 1 4 1 5 2 1 2 389
KCX 169 1P6B.subA 13 0 0 0 0 0 3 0 3 10 3 4 0 4 1 1 4 2 2 2 330
KCX 169 1P6C.subA 15 0 0 0 0 0 3 0 3 12 3 4 0 4 2 1 5 2 2 2 331
KCX 162 1PO9.subA 15 0 0 0 0 0 3 0 3 12 3 1 0 1 5 1 5 2 2 2 373
KCX 162 1POJ.subA 14 0 0 0 0 0 3 0 3 11 3 1 0 1 4 1 5 2 0 2 388
KCX 162 1POK.subA 15 0 0 0 0 0 3 0 3 12 3 1 0 1 4 2 5 2 0 2 377
KCX 205 1PU6.subA 10 0 0 0 1 1 0 2 2 8 2 0 0 0 1 2 5 0 3 2 217
KCX 205 1PU7.subA 10 0 0 0 1 1 0 2 2 8 2 0 0 0 1 2 5 0 3 2 215
KCX 205 1PU8.subA 10 0 0 0 1 1 0 2 2 8 2 0 0 0 1 2 5 0 1 2 216
KCX 169 1QW7.subA 14 0 0 0 0 0 3 0 3 11 3 4 0 4 2 1 4 2 2 2 336
KCX 122 1RCQ.subA 13 0 0 0 0 1 2 1 3 10 3 1 0 1 1 1 7 0 2 2 357
KCX 184 1RQB.subA 11 1 0 1 0 1 1 1 2 7 3 1 0 1 2 1 3 1 5 2 472
KCX 184 1RQH.subA 13 1 0 1 1 1 1 2 3 9 4 1 0 1 2 2 4 1 3 2 471
KCX 184 1RR2.subA 13 2 0 2 1 1 1 2 3 8 5 1 0 1 2 1 4 1 3 2 472
KCX 201 1RXO.subB 12 2 1 3 0 0 3 0 3 6 6 2 1 3 1 1 1 0 1 2 455
KCX 184 1S3H.subA 11 1 0 1 0 1 1 1 2 8 3 1 0 1 2 1 4 1 1 2 472
KCX 220 1S3T.subC 14 0 0 0 0 0 4 0 4 10 4 2 0 2 2 1 5 2 1 2 570
KCX 184 1U5J.subA 11 1 0 1 0 1 1 1 2 8 3 1 0 1 2 1 4 1 2 2 471
KCX 198 1UAG.subA 15 0 0 0 1 1 0 2 2 13 2 4 1 5 2 1 5 0 2 2 428
KCX 220 1UBP.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 2 1 5 2 0 2 570
KCX 201 1UPM.subB 12 2 1 3 0 0 3 0 3 6 6 2 1 3 1 1 1 0 1 2 467
KCX 201 1UPP.subA 12 2 1 3 0 0 3 0 3 6 6 2 1 3 1 1 1 0 0 2 467
KCX 201 1UW9.subA 13 2 1 3 0 0 4 0 4 6 7 2 1 3 1 1 1 1 1 2 465
KCX 201 1UWA.subA 12 2 1 3 0 0 3 0 3 6 6 2 1 3 1 1 1 1 1 2 465
KCX 201 1UZD.subA 12 2 1 3 0 0 3 0 3 6 6 2 1 3 1 1 1 1 0 2 465
KCX 201 1UZH.subA 13 2 1 3 0 0 4 0 4 6 7 2 1 3 1 1 1 1 1 2 465
KCX 129 1VFH.subA 12 0 0 0 0 1 1 1 2 10 2 3 1 4 3 1 2 0 4 2 382
KCX 129 1VFS.subA 12 0 0 0 0 1 1 1 2 10 2 3 2 5 2 1 2 0 3 2 383
KCX 129 1VFT.subA 10 0 0 0 0 0 1 0 1 9 1 3 1 4 2 1 2 0 2 2 382
KCX 188 1W78.subA 12 0 1 1 0 0 0 0 0 11 1 1 0 1 5 0 5 0 4 2 414
KCX 188 1W7K.subA 12 0 1 1 0 0 0 0 0 11 1 1 0 1 5 0 5 0 4 2 411
KCX 201 1WDD.subA 13 2 1 3 0 0 3 0 3 7 6 2 1 3 1 1 2 1 1 2 465
KCX 102 1XGE.subA 13 0 1 1 0 0 2 0 2 10 3 1 0 1 3 2 4 2 1 2 343
KCX 129 1XQK.subA 13 1 0 1 0 1 2 1 3 9 4 2 0 2 1 2 4 0 3 2 382
KCX 129 1XQL.subA 13 1 0 1 0 1 2 1 3 9 4 2 0 2 1 2 4 0 2 2 382
KCX 162 1YBQ.subB 14 0 0 0 0 0 3 0 3 11 3 1 0 1 4 1 5 2 1 2 389
KCX 162 2AQO.subA 14 0 0 0 0 0 3 0 3 11 3 1 0 1 4 1 5 2 1 2 378
KCX 162 2AQV.subA 14 0 0 0 0 0 3 0 3 11 3 1 0 1 4 1 5 2 1 2 375
KCX 102 2E25.subA 12 0 1 1 0 0 2 0 2 9 3 1 0 1 3 1 4 2 0 2 343
KCX 102 2EG6.subA 14 0 1 1 0 0 3 0 3 10 4 1 0 1 3 2 4 2 1 2 343
KCX 102 2EG7.subA 14 0 1 1 0 0 3 0 3 10 4 1 0 1 3 2 4 2 0 2 343
KCX 102 2EG8.subA 13 0 1 1 0 0 3 0 3 9 4 1 0 1 3 1 4 2 1 2 343
KCX 158 2FTW.subA 13 0 0 0 0 0 3 0 3 10 3 0 1 1 2 3 4 2 0 2 484
KCX 167 2FTY.subA 14 0 0 0 0 0 5 0 5 9 5 0 0 0 0 2 7 2 1 2 532
KCX 167 2FVK.subA 14 0 0 0 0 0 5 0 5 9 5 0 0 0 0 2 7 2 0 2 532
KCX 167 2FVM.subA 12 0 0 0 0 0 3 0 3 9 3 0 0 0 0 2 7 2 0 2 532
KCX 185 2GC5.subA 14 0 1 1 0 0 2 0 2 11 3 2 0 2 5 0 4 0 1 2 408
KCX 185 2GC6.subA 14 0 1 1 0 0 2 0 2 11 3 2 0 2 5 0 4 0 2 2 407
KCX 149 2GWN.subA 13 0 0 0 0 0 2 0 2 10 2 1 1 2 2 3 3 2 2 2 451
KCX 154 2ICS.subA 13 0 0 0 0 1 3 1 4 8 4 1 1 2 3 0 3 2 1 2 368
KCX 229 2J6V.subA 14 1 0 1 0 3 1 3 4 8 5 1 0 1 1 1 5 0 1 2 281
KCX 198 2JFF.subA 14 0 0 0 0 1 0 1 1 13 1 4 1 5 2 1 5 0 2 2 433
KCX 198 2JFG.subA 15 0 0 0 1 1 0 2 2 13 2 4 1 5 2 1 5 0 2 2 440
KCX 198 2JFH.subA 13 0 0 0 0 1 0 1 1 12 1 4 1 5 2 1 4 0 2 2 431
KCX 217 2KAU.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 3 1 4 2 0 2 566
KCX 169 2O4M.subA 12 0 0 0 0 0 3 0 3 9 3 4 0 4 1 1 3 2 2 2 330
KCX 169 2O4Q.subA 13 0 0 0 0 0 3 0 3 10 3 4 0 4 1 1 4 2 1 2 331
KCX 169 2OB3.subA 13 0 0 0 0 0 3 0 3 10 3 4 0 4 1 1 4 2 4 2 330
KCX 122 2ODO.subA 14 0 0 0 0 1 2 1 3 10 3 1 0 1 1 1 7 0 0 2 356
KCX 173 2OEK.subA 13 2 1 3 0 0 1 0 1 9 4 1 0 1 1 0 7 1 2 2 412
KCX 173 2OEL.subA 13 2 1 3 0 0 1 0 1 9 4 1 0 1 1 0 7 1 2 2 412
KCX 173 2OEM.subA 13 2 1 3 0 0 1 0 1 9 4 1 1 2 1 0 6 1 1 2 408
KCX 175 2OGJ.subA 13 0 0 0 0 1 3 1 4 8 4 1 1 2 2 1 3 2 0 2 379
KCX 169 2OQL.subA 13 0 0 0 0 0 3 0 3 10 3 4 0 4 1 1 4 2 3 2 330
KCX 315 2P9V.subA 12 0 1 1 1 0 1 1 2 9 3 3 0 3 1 3 2 0 4 2 355
KCX 718 2QF7.subA 12 1 0 1 0 2 1 2 3 8 4 0 0 0 2 1 5 1 3 2 1076
KCX 166 2QPX.subA 13 1 0 1 0 1 3 1 4 8 5 1 1 2 0 3 3 2 2 2 376
KCX 169 2R1K.subA 14 0 0 0 0 0 3 0 3 11 3 4 0 4 1 1 5 2 1 2 329
KCX 169 2R1L.subA 13 0 0 0 0 0 3 0 3 10 3 4 0 4 1 1 4 2 1 2 329
KCX 169 2R1M.subA 14 0 0 0 0 0 3 0 3 11 3 4 0 4 2 1 4 2 1 2 328
KCX 169 2R1N.subA 14 0 0 0 0 0 3 0 3 11 3 4 0 4 2 1 4 2 2 2 328
KCX 169 2R1P.subA 14 0 0 0 0 0 3 0 3 11 3 4 0 4 2 1 4 2 1 2 328
KCX 122 2RJG.subA 13 0 0 0 0 1 2 1 3 10 3 1 0 1 2 1 6 0 3 2 359
KCX 122 2RJH.subA 13 0 0 0 0 1 2 1 3 10 3 1 0 1 2 1 6 0 2 2 359
KCX 70 2RL3.subA 12 0 0 0 1 0 1 1 2 10 2 3 0 3 2 2 3 0 3 2 246
KCX 129 2SFP.subA 13 1 0 1 0 1 2 1 3 9 4 2 0 2 1 2 4 0 3 2 379
KCX 198 2UAG.subA 15 0 0 0 1 1 0 2 2 13 2 4 1 5 2 1 5 0 2 2 429
KCX 220 2UBP.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 2 1 5 2 0 2 570
KCX 198 2UUO.subA 15 0 0 0 1 1 0 2 2 13 2 4 1 5 2 1 5 0 2 2 430
KCX 198 2UUP.subA 14 0 0 0 0 1 0 1 1 13 1 4 1 5 2 1 5 0 2 2 440
KCX 58 2UYN.subA 9 0 1 1 0 1 1 1 2 6 3 0 0 0 1 1 4 0 2 1 127
KCX 201 2V63.subA 12 2 1 3 0 0 3 0 3 6 6 2 1 3 1 1 1 1 1 2 466
KCX 201 2V67.subA 12 2 1 3 0 0 3 0 3 6 6 2 1 3 1 1 1 1 1 2 466
KCX 201 2V68.subA 12 2 1 3 0 0 3 0 3 6 6 2 1 3 1 1 1 1 1 2 466
KCX 201 2V69.subA 12 2 1 3 0 0 3 0 3 6 6 2 1 3 1 1 1 1 0 2 459
KCX 201 2V6A.subA 12 2 1 3 0 0 3 0 3 6 6 2 1 3 1 1 1 1 1 2 467
KCX 137 2VC5.subA 14 0 0 0 0 0 3 0 3 11 3 2 0 2 2 2 5 2 1 2 314
KCX 137 2VC7.subA 14 0 0 0 0 0 3 0 3 11 3 2 0 2 2 2 5 2 0 2 314
KCX 201 2VDH.subA 12 2 1 3 0 0 3 0 3 6 6 3 1 4 0 1 1 1 1 2 465
KCX 201 2VDI.subA 13 2 1 3 0 0 4 0 4 6 7 2 1 3 1 1 1 1 1 2 465
KCX 198 2VTD.subA 14 0 0 0 0 1 0 1 1 13 1 4 1 5 2 1 5 0 2 2 435
KCX 198 2VTE.subA 14 0 0 0 1 1 0 2 2 12 2 4 1 5 2 1 4 0 1 2 434
KCX 70 2WGW.subA 13 0 0 0 2 0 0 2 2 11 2 4 0 4 2 3 2 0 3 2 246
KCX 198 2WJP.subA 15 0 0 0 1 1 0 2 2 13 2 4 1 5 2 1 5 0 3 2 435
KCX 262 2WTZ.subA 14 0 0 0 1 0 0 1 1 13 1 4 0 4 3 3 3 0 0 2 504
KCX 70 2X01.subA 13 0 0 0 2 0 0 2 2 11 2 2 0 2 3 3 3 0 3 2 244
KCX 70 2X02.subA 12 0 0 0 1 0 0 1 1 11 1 3 0 3 2 3 3 0 1 2 245
KCX 198 2X5O.subA 15 0 0 0 1 1 0 2 2 13 2 4 1 5 2 1 5 0 2 2 439
KCX 262 2XJA.subA 15 0 0 0 1 0 0 1 1 14 1 4 0 4 3 3 4 0 0 2 502
KCX 198 2XPC.subA 15 0 0 0 1 1 0 2 2 13 2 4 1 5 2 1 5 0 2 2 435
KCX 198 2Y66.subA 15 0 0 0 1 1 0 2 2 13 2 4 1 5 2 1 5 0 2 2 434
KCX 198 2Y67.subA 15 0 0 0 1 1 0 2 2 13 2 4 1 5 2 1 5 0 3 2 435
KCX 198 2Y68.subA 14 0 0 0 0 1 0 1 1 13 1 4 1 5 2 1 5 0 2 2 435
KCX 102 2Z24.subA 14 0 1 1 0 0 3 0 3 10 4 1 0 1 3 2 4 2 1 2 343
KCX 102 2Z25.subA 14 0 1 1 0 0 3 0 3 10 4 1 0 1 3 2 4 2 1 2 343
KCX 102 2Z26.subA 13 0 1 1 0 0 2 0 2 10 3 1 0 1 3 2 4 2 1 2 344
KCX 102 2Z27.subA 14 0 1 1 0 0 3 0 3 10 4 1 0 1 3 2 4 2 1 2 343
KCX 102 2Z28.subA 14 0 1 1 0 0 3 0 3 10 4 1 0 1 3 2 4 2 1 2 343
KCX 102 2Z29.subA 14 0 1 1 0 0 3 0 3 10 4 1 0 1 3 2 4 2 1 2 343
KCX 102 2Z2A.subA 14 0 1 1 0 0 3 0 3 10 4 1 0 1 3 2 4 2 1 2 343
KCX 102 2Z2B.subA 13 0 1 1 0 0 2 0 2 10 3 1 0 1 3 2 4 2 1 2 331
KCX 243 2ZC1.subA 14 0 0 0 0 0 3 0 3 11 3 2 0 2 5 1 3 2 1 2 333
KCX 189 3A12.subA 13 2 1 3 0 0 3 0 3 7 6 0 1 1 0 2 4 1 0 2 437
KCX 189 3A13.subA 13 2 1 3 0 0 3 0 3 7 6 0 1 1 0 2 4 1 0 2 436
KCX 169 3A3W.subA 14 0 0 0 0 0 3 0 3 11 3 4 0 4 2 1 4 2 2 2 329
KCX 169 3A3X.subA 14 0 0 0 0 0 3 0 3 11 3 4 0 4 2 1 4 2 2 2 329
KCX 169 3A4J.subA 13 0 0 0 0 0 3 0 3 10 3 4 0 4 1 1 4 2 2 2 328
KCX 122 3B8T.subA 13 0 0 0 0 1 2 1 3 10 3 1 0 1 2 1 6 0 2 2 359
KCX 122 3B8U.subA 13 0 0 0 0 1 2 1 3 10 3 1 0 1 2 1 6 0 2 2 359
KCX 122 3B8V.subA 13 0 0 0 0 1 2 1 3 10 3 1 0 1 2 1 6 0 3 2 359
KCX 122 3B8W.subA 12 0 0 0 0 1 2 1 3 9 3 0 0 0 2 1 6 0 2 2 359
KCX 741 3BG3.subA 13 1 0 1 0 2 1 2 3 9 4 1 0 1 2 1 5 1 0 2 680
KCX 229 3C0Q.subA 13 1 0 1 0 3 1 3 4 8 5 1 0 1 1 1 5 0 1 2 277
KCX 229 3C0S.subA 14 1 0 1 0 3 1 3 4 8 5 1 0 1 1 1 5 0 1 2 281
KCX 169 3C86.subA 14 0 0 0 0 0 3 0 3 11 3 4 0 4 2 1 4 2 2 2 328
KCX 169 3CAK.subA 13 0 0 0 0 0 3 0 3 10 3 4 0 4 1 1 4 2 1 2 330
KCX 169 3CS2.subA 13 0 0 0 0 0 3 0 3 10 3 4 0 4 1 1 4 2 1 2 331
KCX 147 3DC8.subA 13 0 0 0 0 0 5 0 5 8 5 0 0 0 1 3 4 2 1 2 483
KCX 182 3DUG.subA 10 0 1 1 0 0 3 0 3 6 4 1 0 1 2 1 2 2 0 2 394
KCX 169 3E3H.subA 13 0 0 0 0 0 3 0 3 10 3 4 0 4 2 1 3 2 1 2 336
KCX 146 3E74.subA 12 0 0 0 0 0 3 0 3 8 3 0 0 0 4 2 2 2 2 2 429
KCX 243 3F4C.subA 12 0 0 0 0 0 2 0 2 10 2 2 0 2 4 1 3 2 2 2 324
KCX 145 3F4D.subA 13 0 0 0 0 0 3 0 3 10 3 2 0 2 4 1 3 2 0 2 324
KCX 143 3FDK.subA 14 0 0 0 0 0 3 0 3 11 3 2 0 2 5 1 3 2 1 2 322
KCX 84 3FV7.subA 12 0 0 0 1 0 0 1 1 11 1 3 0 3 2 3 3 0 1 2 244
KCX 84 3FYZ.subA 13 0 0 0 1 0 0 1 1 12 1 3 1 4 2 3 3 0 1 2 244
KCX 84 3FZC.subA 12 0 0 0 1 0 0 1 1 11 1 3 0 3 2 3 3 0 1 2 244
KCX 84 3G4P.subA 12 0 0 0 1 0 0 1 1 11 1 3 0 3 2 3 3 0 1 2 244
KCX 277 3GTF.subA 14 0 0 0 0 0 3 0 3 11 3 2 0 2 5 1 3 2 1 2 331
KCX 243 3GTH.subA 14 0 0 0 0 0 3 0 3 11 3 2 0 2 5 1 3 2 1 2 331
KCX 243 3GTI.subA 13 0 0 0 0 0 3 0 3 10 3 2 0 2 4 1 3 2 0 2 331
KCX 243 3GTX.subA 14 0 0 0 0 0 3 0 3 11 3 2 0 2 5 1 3 2 1 2 333
KCX 243 3GU1.subA 15 0 0 0 0 0 3 0 3 12 3 3 0 3 5 1 3 2 1 2 331
KCX 243 3GU2.subA 12 0 0 0 0 0 3 0 3 9 3 2 0 2 4 0 3 2 1 2 330
KCX 243 3GU9.subA 13 0 0 0 0 0 3 0 3 10 3 2 0 2 4 1 3 2 1 2 331
KCX 73 3HBR.subA 13 0 0 0 1 0 0 1 1 12 1 4 0 4 1 3 4 0 2 2 240
KCX 243 3HTW.subA 13 0 0 0 0 0 3 0 3 10 3 2 0 2 5 0 3 2 0 2 322
KCX 294 3ICJ.subA 13 0 1 1 0 0 3 0 3 9 4 1 0 1 1 2 5 2 1 2 468
KCX 248 3IF6.subA 10 2 0 2 0 1 0 1 1 7 3 0 2 2 1 1 3 0 2 2 228
KCX 70 3ISG.subA 13 0 0 0 1 0 0 1 1 12 1 4 0 4 0 3 5 0 0 2 249
KCX 103 3JZE.subA 13 0 1 1 0 0 2 0 2 10 3 1 0 1 3 2 4 2 1 2 344
KCX 189 3KDN.subA 13 2 1 3 0 0 3 0 3 7 6 0 1 1 0 2 4 1 0 2 437
KCX 189 3KDO.subA 13 2 1 3 0 0 3 0 3 7 6 0 1 1 0 2 4 1 0 2 436
KCX 302 3KZC.subA 12 0 0 0 1 1 2 2 4 8 4 2 0 2 2 1 3 0 0 2 328
KCX 302 3KZK.subA 12 0 0 0 1 1 2 2 4 8 4 2 0 2 2 1 3 0 1 2 334
KCX 302 3KZM.subA 12 0 0 0 1 1 2 2 4 8 4 2 0 2 2 1 3 0 1 2 332
KCX 302 3KZN.subA 12 0 0 0 1 1 2 2 4 8 4 2 0 2 2 1 3 0 2 2 332
KCX 302 3KZO.subA 12 0 0 0 1 1 2 2 4 8 4 2 0 2 2 1 3 0 2 2 332
KCX 302 3L02.subA 12 0 0 0 1 1 2 2 4 8 4 2 0 2 2 1 3 0 2 2 332
KCX 302 3L04.subA 12 0 0 0 1 1 2 2 4 8 4 2 0 2 2 1 3 0 2 2 332
KCX 302 3L05.subA 12 0 0 0 1 1 2 2 4 8 4 2 0 2 2 1 3 0 1 2 332
KCX 302 3L06.subA 12 0 0 0 1 1 2 2 4 8 4 2 0 2 2 1 3 0 1 2 332
KCX 490 3LA4.subA 13 0 0 0 0 0 3 0 3 10 3 3 1 4 2 1 3 2 1 2 837
KCX 70 3LCE.subA 12 0 0 0 1 0 0 1 1 11 1 3 0 3 2 3 3 0 1 2 245
KCX 302 3M4J.subA 12 0 0 0 1 1 2 2 4 8 4 2 0 2 2 1 3 0 1 2 332
KCX 84 3MBZ.subA 12 0 0 0 1 0 0 1 1 11 1 3 0 3 2 3 3 0 1 2 244
KCX 102 3MJM.subA 14 0 1 1 0 0 3 0 3 10 4 1 0 1 3 2 4 2 1 2 344
KCX 191 3MKV.subA 13 0 1 1 0 0 3 0 3 9 4 1 0 1 3 0 5 2 1 2 414
KCX 211 3MTW.subA 12 0 1 1 0 0 3 0 3 8 4 1 0 1 4 0 3 2 3 2 403
KCX 188 3N2C.subA 14 0 1 1 1 0 3 1 4 9 5 1 0 1 3 0 5 2 0 2 408
KCX 195 3NWR.subA 13 2 1 3 0 0 3 0 3 7 6 1 0 1 1 1 4 0 3 2 404
KCX 145 3OJG.subA 14 0 0 0 0 0 3 0 3 11 3 2 0 2 4 2 3 2 0 2 323
KCX 169 3OOD.subA 14 0 0 0 0 0 3 0 3 11 3 4 0 4 2 1 4 2 2 2 329
KCX 169 3OQE.subA 13 0 0 0 0 0 3 0 3 10 3 4 0 4 1 1 4 2 2 2 329
KCX 145 3ORW.subA 13 0 0 0 0 0 3 0 3 10 3 2 0 2 4 1 3 2 0 2 324
KCX 153 3OVG.subA 13 1 0 1 0 0 3 0 3 8 4 2 0 2 3 1 2 2 1 2 351
KCX 92 3PNU.subA 12 0 1 1 0 0 2 0 2 9 3 1 0 1 1 1 6 2 0 2 338
KCX 154 3PNZ.subA 13 0 0 0 0 0 4 0 4 9 4 2 1 3 3 2 1 2 1 2 329
KCX 62 3Q7V.subB 12 0 0 0 1 0 0 1 1 11 1 2 1 3 0 3 5 0 2 2 252
KCX 218 3QGA.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 2 2 4 2 1 2 564
KCX 218 3QGK.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 2 2 4 2 0 2 564
KCX 70 3QNB.subA 12 0 0 0 1 0 0 1 1 11 1 3 0 3 2 3 3 0 2 2 242
KCX 70 3QNC.subA 12 0 0 0 1 0 0 1 1 11 1 3 0 3 2 3 3 0 3 2 242
KCX 129 3S46.subA 12 1 0 1 0 1 1 1 2 9 3 2 0 2 2 1 4 0 3 2 367
KCX 169 3SO7.subA 14 0 0 0 0 0 3 0 3 11 3 4 0 4 2 1 4 2 1 2 329
KCX 718 3TW6.subA 12 1 0 1 0 2 1 2 3 8 4 0 0 0 2 1 5 1 3 2 1007
KCX 718 3TW7.subA 12 1 0 1 0 2 1 2 3 8 4 0 0 0 2 1 5 1 0 2 1004
KCX 198 3UAG.subA 15 0 0 0 1 1 0 2 2 13 2 4 1 5 2 1 5 0 2 2 430
KCX 220 3UBP.subC 13 0 0 0 0 0 3 0 3 10 3 2 0 2 2 1 5 2 0 2 570
KCX 198 4UAG.subA 14 0 0 0 0 1 0 1 1 13 1 4 1 5 2 1 5 0 2 2 429
KCX 220 4UBP.subC 14 0 0 0 0 0 4 0 4 10 4 2 0 2 2 1 5 2 0 2 570
KCX 201 8RUC.subA 13 2 1 3 0 0 3 0 3 7 6 2 1 3 1 1 2 1 1 2 467
";
	my $sizecontrol = 0;
	my $c = 0;
	
	my %hashofarrays = ();

	my @data = split('\n',$kcxtraining);

	foreach my $line (@data)
	{
		chomp($line);
		my @temp = split(" ",$line);
		my $size = @temp;
		#Size control for the first one
		if ($c == 0)
		{
			$sizecontrol = $size;
		}
		
		#All the others - -- - - - -- ---
		if($size != $sizecontrol)
		{
			print "\n\nFATAL ERROR PROBLEM\n";
			print $command." error 518\n";
			print "Please, contact us at prelyscar\@gmail.com\n\n";
			exit;
		}
		#- - - - - - -- - -end of control
		
		$hashofarrays{$c} = [@temp];
		
		$c++;
	}
	return (\%hashofarrays);
}

# Loading DATA
# LYS training data
sub loading_lys_training
{
	my $lystraining = "LYS 21 1AA1.subB 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 0 2 438
LYS 81 1AA1.subB 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 0 2 438
LYS 128 1AA1.subB 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 0 2 438
LYS 146 1AA1.subB 5 0 0 0 0 0 0 0 0 5 0 1 0 1 1 0 3 0 3 2 438
LYS 161 1AA1.subB 6 1 1 2 0 0 1 0 1 3 3 0 0 0 1 0 2 0 3 2 438
LYS 164 1AA1.subB 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 2 2 438
LYS 175 1AA1.subB 6 0 0 0 0 0 0 0 0 6 0 1 0 1 3 0 2 0 1 2 438
LYS 177 1AA1.subB 7 1 1 2 0 0 0 0 0 5 2 0 1 1 1 0 3 1 5 2 438
LYS 183 1AA1.subB 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 3 2 438
LYS 227 1AA1.subB 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 2 0 0 3 2 438
LYS 236 1AA1.subB 12 1 1 2 0 1 0 1 1 9 3 1 1 2 3 0 4 0 1 2 438
LYS 252 1AA1.subB 5 1 1 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 4 2 438
LYS 305 1AA1.subB 6 0 1 1 0 2 0 2 2 3 3 0 3 3 0 0 0 0 2 2 438
LYS 316 1AA1.subB 11 1 0 1 0 2 0 2 2 8 3 0 1 1 2 1 4 0 1 2 438
LYS 356 1AA1.subB 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 2 2 438
LYS 450 1AA1.subB 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 0 2 438
LYS 463 1AA1.subB 3 0 0 0 0 0 0 0 0 3 0 0 0 0 1 2 0 0 0 2 438
LYS 9 1BWV.subA 8 0 1 1 0 1 0 1 1 6 2 1 1 2 0 2 2 0 3 2 471
LYS 22 1BWV.subA 4 0 0 0 0 0 0 0 0 4 0 0 0 0 2 0 2 0 0 2 471
LYS 32 1BWV.subA 8 2 1 3 0 0 0 0 0 5 3 1 1 2 0 1 2 0 3 2 471
LYS 83 1BWV.subA 11 1 0 1 0 1 0 1 1 9 2 0 1 1 4 3 1 0 0 2 471
LYS 86 1BWV.subA 5 1 0 1 0 1 0 1 1 3 2 0 0 0 0 2 1 0 2 2 471
LYS 128 1BWV.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 1 0 0 1 2 471
LYS 131 1BWV.subA 7 0 1 1 0 0 0 0 0 6 1 1 1 2 3 0 1 0 0 2 471
LYS 146 1BWV.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 2 0 2 0 2 2 471
LYS 164 1BWV.subA 10 2 1 3 1 2 0 3 3 4 6 0 0 0 0 1 3 0 1 2 471
LYS 175 1BWV.subA 8 1 0 1 1 0 0 1 1 6 2 2 0 2 3 0 1 1 1 2 471
LYS 177 1BWV.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 1 2 471
LYS 183 1BWV.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 1 1 0 0 2 2 471
LYS 194 1BWV.subA 8 0 1 1 0 0 0 0 0 7 1 1 1 2 3 1 1 0 0 2 471
LYS 227 1BWV.subA 5 0 1 1 0 0 0 0 0 4 1 0 1 1 2 1 0 0 0 2 471
LYS 236 1BWV.subA 13 1 1 2 1 1 0 2 2 9 4 2 0 2 3 0 4 0 1 2 471
LYS 258 1BWV.subA 8 0 1 1 0 0 0 0 0 7 1 0 2 2 2 2 1 0 0 2 471
LYS 282 1BWV.subA 6 1 0 1 0 0 0 0 0 5 1 1 2 3 1 1 0 0 3 2 471
LYS 305 1BWV.subA 9 2 0 2 0 2 0 2 2 5 4 3 2 5 0 0 0 0 1 2 471
LYS 316 1BWV.subA 13 1 0 1 0 2 0 2 2 10 3 1 0 1 1 2 6 0 1 2 471
LYS 334 1BWV.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 2 2 471
LYS 347 1BWV.subA 9 0 1 1 0 1 0 1 1 7 2 1 1 2 1 1 3 0 0 2 471
LYS 353 1BWV.subA 6 1 2 3 0 0 0 0 0 3 3 0 0 0 1 0 2 0 1 2 471
LYS 373 1BWV.subA 8 0 0 0 0 2 0 2 2 6 2 0 0 0 1 0 5 0 1 2 471
LYS 450 1BWV.subA 5 0 1 1 0 2 0 2 2 2 3 1 0 1 1 0 0 0 0 2 471
LYS 463 1BWV.subA 5 2 0 2 0 0 0 0 0 3 2 0 0 0 0 1 2 0 0 2 471
LYS 30 1E4D.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 1 0 0 0 2 243
LYS 45 1E4D.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 1 1 1 0 1 2 243
LYS 49 1E4D.subA 3 0 0 0 0 0 0 0 0 3 0 3 0 3 0 0 0 0 3 2 243
LYS 61 1E4D.subA 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 4 2 243
LYS 84 1E4D.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 4 2 243
LYS 91 1E4D.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 5 2 243
LYS 95 1E4D.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 6 2 243
LYS 100 1E4D.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 2 2 243
LYS 134 1E4D.subA 5 0 0 0 1 1 0 2 2 3 2 0 1 1 0 1 1 0 8 2 243
LYS 137 1E4D.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 5 2 243
LYS 138 1E4D.subA 4 0 0 0 1 0 0 1 1 3 1 1 0 1 0 2 0 0 2 2 243
LYS 152 1E4D.subA 7 1 1 2 0 0 0 0 0 5 2 0 1 1 2 1 1 0 2 2 243
LYS 177 1E4D.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 0 0 2 0 9 2 243
LYS 182 1E4D.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 0 1 0 5 2 243
LYS 189 1E4D.subA 10 0 1 1 0 0 1 0 1 8 2 0 2 2 0 2 4 0 1 2 243
LYS 205 1E4D.subA 11 0 0 0 0 0 0 0 0 11 0 5 0 5 1 0 4 0 3 2 243
LYS 228 1E4D.subA 10 0 4 4 0 0 0 0 0 6 4 0 0 0 1 2 3 0 1 2 243
LYS 246 1E4D.subA 8 0 1 1 0 0 0 0 0 7 1 1 1 2 1 0 4 0 9 2 243
LYS 251 1E4D.subA 6 0 0 0 0 1 0 1 1 5 1 3 0 3 1 1 0 0 8 2 243
LYS 256 1E4D.subA 5 0 1 1 0 0 0 0 0 4 1 2 0 2 0 0 2 0 8 2 243
LYS 67 1E8C.subA 4 1 1 2 0 1 0 1 1 1 3 0 0 0 1 0 0 0 6 2 496
LYS 119 1E8C.subA 11 0 1 1 0 0 0 0 0 10 1 5 0 5 2 0 2 0 3 2 496
LYS 150 1E8C.subA 4 0 0 0 0 0 0 0 0 4 0 0 0 0 1 0 3 0 2 2 496
LYS 197 1E8C.subA 6 1 0 1 0 0 1 0 1 4 2 1 0 1 1 1 1 0 2 2 496
LYS 251 1E8C.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 1 2 1 0 3 2 496
LYS 274 1E8C.subA 7 0 0 0 0 0 0 0 0 7 0 3 0 3 1 1 2 0 1 2 496
LYS 330 1E8C.subA 5 1 0 1 0 1 0 1 1 3 2 1 0 1 1 0 1 0 1 2 496
LYS 350 1E8C.subA 7 0 0 0 0 0 0 0 0 6 0 0 0 0 4 0 2 0 0 2 496
LYS 366 1E8C.subA 7 1 1 2 0 0 0 0 0 5 2 0 2 2 2 1 0 0 4 2 496
LYS 378 1E8C.subA 6 2 0 2 0 0 0 0 0 4 2 0 0 0 1 1 2 0 1 2 496
LYS 391 1E8C.subA 6 2 1 3 0 2 0 2 2 1 5 0 0 0 1 0 0 0 0 2 496
LYS 393 1E8C.subA 9 2 0 2 0 2 1 2 3 3 5 1 0 1 2 0 0 0 3 2 496
LYS 438 1E8C.subA 8 0 0 0 0 0 1 0 1 7 1 0 1 1 2 1 3 0 0 2 496
LYS 455 1E8C.subA 5 1 1 2 0 0 0 0 0 3 2 0 2 2 1 0 0 0 3 2 496
LYS 465 1E8C.subA 11 2 1 3 0 1 0 1 1 7 4 0 2 2 2 2 1 0 2 2 496
LYS 2 1E9Y.subB 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 2 569
LYS 3 1E9Y.subB 2 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 2 0 0 2 569
LYS 7 1E9Y.subB 3 0 1 1 0 1 0 1 1 1 2 1 0 1 0 0 0 0 0 2 569
LYS 20 1E9Y.subB 8 2 1 3 0 0 0 0 0 5 3 0 1 1 2 0 2 0 0 2 569
LYS 44 1E9Y.subB 9 0 2 2 1 0 0 1 1 6 3 1 1 2 0 2 2 0 1 2 569
LYS 49 1E9Y.subB 7 0 0 0 1 0 0 1 1 6 1 1 0 1 3 1 1 0 0 2 569
LYS 63 1E9Y.subB 3 0 1 1 0 0 0 0 0 2 1 1 0 1 1 0 0 0 0 2 569
LYS 83 1E9Y.subB 14 1 0 1 0 0 1 0 1 12 2 2 1 3 4 3 2 0 0 2 569
LYS 89 1E9Y.subB 7 2 2 4 0 0 0 0 0 3 4 0 0 0 1 0 2 0 0 2 569
LYS 92 1E9Y.subB 8 1 0 1 0 0 0 0 0 7 1 0 1 1 4 0 2 0 2 2 569
LYS 98 1E9Y.subB 14 2 2 4 0 0 0 0 0 10 4 1 1 2 4 1 3 0 2 2 569
LYS 102 1E9Y.subB 5 1 0 1 0 0 0 0 0 4 1 0 3 3 0 0 1 0 0 2 569
LYS 109 1E9Y.subB 4 0 0 0 0 0 0 0 0 4 0 0 2 2 0 0 2 0 0 2 569
LYS 180 1E9Y.subB 6 0 0 0 0 3 0 3 3 3 3 0 0 0 0 1 2 0 0 2 569
LYS 198 1E9Y.subB 7 0 0 0 0 0 0 0 0 7 0 0 2 2 3 1 1 0 0 2 569
LYS 240 1E9Y.subB 4 2 0 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 0 2 569
LYS 288 1E9Y.subB 8 1 1 2 0 0 0 0 0 6 2 0 0 0 0 0 6 0 1 2 569
LYS 324 1E9Y.subB 9 2 1 3 0 0 2 0 2 4 5 0 0 0 1 1 2 0 0 2 569
LYS 326 1E9Y.subB 2 1 0 1 0 0 0 0 0 1 1 1 0 1 0 0 0 0 0 2 569
LYS 329 1E9Y.subB 5 0 1 1 0 0 0 0 0 4 1 1 1 2 0 0 2 0 0 2 569
LYS 382 1E9Y.subB 8 1 0 1 0 0 0 0 0 7 1 2 2 4 2 0 1 0 1 2 569
LYS 384 1E9Y.subB 11 3 0 3 1 2 0 3 3 5 6 0 1 1 2 0 2 0 0 2 569
LYS 385 1E9Y.subB 5 3 1 4 1 0 0 1 1 0 5 0 0 0 0 0 0 0 0 2 569
LYS 391 1E9Y.subB 2 0 1 1 0 0 0 0 0 1 1 0 0 0 0 0 1 0 1 2 569
LYS 394 1E9Y.subB 3 0 2 2 0 0 0 0 0 1 2 0 0 0 1 0 0 0 1 2 569
LYS 403 1E9Y.subB 12 1 1 2 0 1 0 1 1 9 3 2 1 3 1 3 2 0 1 2 569
LYS 408 1E9Y.subB 14 1 0 1 0 1 1 1 2 11 3 4 2 6 0 3 2 0 0 2 569
LYS 430 1E9Y.subB 12 0 2 2 0 0 0 0 0 10 2 2 0 2 2 1 5 0 1 2 569
LYS 445 1E9Y.subB 7 0 0 0 0 0 0 0 0 7 0 1 0 1 1 3 2 0 1 2 569
LYS 451 1E9Y.subB 12 1 0 1 0 0 0 0 0 11 1 1 2 3 3 2 3 0 0 2 569
LYS 484 1E9Y.subB 4 0 1 1 0 0 0 0 0 3 1 0 0 0 2 1 0 0 0 2 569
LYS 486 1E9Y.subB 11 0 2 2 0 0 2 0 2 7 4 1 1 2 3 2 0 0 0 2 569
LYS 501 1E9Y.subB 7 1 1 2 0 1 0 1 1 4 3 0 0 0 3 0 1 0 0 2 569
LYS 504 1E9Y.subB 8 0 2 2 0 1 0 1 1 5 3 0 1 1 1 0 3 0 0 2 569
LYS 517 1E9Y.subB 8 1 1 2 0 0 0 0 0 6 2 0 2 2 2 1 1 0 1 2 569
LYS 524 1E9Y.subB 16 0 1 1 0 0 1 0 1 14 2 4 2 6 5 1 2 0 0 2 569
LYS 525 1E9Y.subB 3 1 0 1 0 0 0 0 0 2 1 1 1 2 0 0 0 0 2 2 569
LYS 550 1E9Y.subB 4 1 1 2 0 0 0 0 0 2 2 0 0 0 1 0 1 0 0 2 569
LYS 555 1E9Y.subB 3 0 1 1 0 0 0 0 0 2 1 1 0 1 1 0 0 0 0 2 569
LYS 559 1E9Y.subB 2 0 0 0 0 0 0 0 0 2 0 0 1 1 0 0 1 0 0 2 569
LYS 2 1E9Z.subB 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 2 569
LYS 3 1E9Z.subB 2 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 2 0 0 2 569
LYS 7 1E9Z.subB 3 0 1 1 0 1 0 1 1 1 2 1 0 1 0 0 0 0 0 2 569
LYS 20 1E9Z.subB 8 2 1 3 0 0 0 0 0 5 3 0 1 1 2 0 2 0 0 2 569
LYS 44 1E9Z.subB 9 0 2 2 1 0 0 1 1 6 3 1 1 2 0 2 2 0 1 2 569
LYS 49 1E9Z.subB 7 0 0 0 1 0 0 1 1 6 1 1 0 1 3 1 1 0 1 2 569
LYS 63 1E9Z.subB 2 0 1 1 0 0 0 0 0 1 1 1 0 1 0 0 0 0 0 2 569
LYS 83 1E9Z.subB 12 1 0 1 0 0 1 0 1 10 2 1 1 2 4 2 2 0 0 2 569
LYS 89 1E9Z.subB 7 2 2 4 0 0 0 0 0 3 4 0 0 0 1 0 2 0 0 2 569
LYS 92 1E9Z.subB 8 1 0 1 0 0 0 0 0 7 1 0 1 1 4 0 2 0 2 2 569
LYS 98 1E9Z.subB 14 1 2 3 1 0 0 1 1 10 4 1 1 2 4 1 3 0 2 2 569
LYS 102 1E9Z.subB 5 1 0 1 0 0 0 0 0 4 1 0 3 3 0 0 1 0 0 2 569
LYS 109 1E9Z.subB 8 0 0 0 1 0 0 1 1 7 1 0 2 2 3 0 2 0 0 2 569
LYS 180 1E9Z.subB 6 0 0 0 0 3 0 3 3 3 3 0 0 0 0 1 2 0 0 2 569
LYS 198 1E9Z.subB 8 0 0 0 0 0 0 0 0 8 0 0 2 2 4 1 1 0 0 2 569
LYS 240 1E9Z.subB 4 2 0 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 0 2 569
LYS 288 1E9Z.subB 8 1 1 2 0 0 0 0 0 6 2 0 0 0 0 0 6 0 1 2 569
LYS 326 1E9Z.subB 6 2 0 2 0 0 0 0 0 4 2 1 0 1 0 0 3 0 0 2 569
LYS 329 1E9Z.subB 4 0 1 1 0 0 0 0 0 3 1 0 1 1 0 0 2 0 1 2 569
LYS 382 1E9Z.subB 10 1 0 1 1 0 0 1 1 8 2 2 2 4 3 0 1 0 0 2 569
LYS 384 1E9Z.subB 10 3 0 3 0 2 0 2 2 5 5 0 1 1 2 0 2 0 1 2 569
LYS 385 1E9Z.subB 6 2 1 3 1 0 0 1 1 2 4 1 1 2 0 0 0 0 0 2 569
LYS 391 1E9Z.subB 2 0 1 1 0 0 0 0 0 1 1 0 0 0 0 0 1 0 1 2 569
LYS 394 1E9Z.subB 3 0 2 2 0 0 0 0 0 1 2 0 0 0 1 0 0 0 1 2 569
LYS 403 1E9Z.subB 11 1 1 2 0 1 0 1 1 8 3 1 1 2 1 3 2 0 1 2 569
LYS 408 1E9Z.subB 13 0 0 0 0 1 1 1 2 11 2 4 2 6 0 3 2 0 0 2 569
LYS 430 1E9Z.subB 13 0 2 2 0 0 0 0 0 11 2 2 0 2 2 2 5 0 1 2 569
LYS 445 1E9Z.subB 6 0 0 0 0 0 0 0 0 6 0 1 0 1 1 3 1 0 1 2 569
LYS 451 1E9Z.subB 13 1 0 1 0 0 0 0 0 12 1 1 2 3 3 2 4 0 0 2 569
LYS 484 1E9Z.subB 4 0 1 1 0 0 0 0 0 3 1 0 0 0 2 1 0 0 1 2 569
LYS 486 1E9Z.subB 11 0 2 2 0 0 2 0 2 7 4 1 1 2 3 2 0 0 0 2 569
LYS 501 1E9Z.subB 7 1 1 2 0 1 0 1 1 4 3 0 0 0 3 0 1 0 0 2 569
LYS 504 1E9Z.subB 8 0 2 2 0 1 0 1 1 5 3 0 1 1 1 0 3 0 0 2 569
LYS 517 1E9Z.subB 7 1 1 2 0 0 0 0 0 5 2 0 1 1 2 1 1 0 1 2 569
LYS 524 1E9Z.subB 16 0 1 1 0 0 1 0 1 14 2 4 2 6 5 1 2 0 0 2 569
LYS 525 1E9Z.subB 3 1 0 1 0 0 0 0 0 2 1 1 1 2 0 0 0 0 2 2 569
LYS 550 1E9Z.subB 5 0 1 1 0 0 0 0 0 4 1 0 0 0 1 1 2 0 0 2 569
LYS 555 1E9Z.subB 3 0 1 1 0 0 0 0 0 2 1 1 0 1 1 0 0 0 0 2 569
LYS 559 1E9Z.subB 2 0 0 0 0 0 0 0 0 2 0 0 1 1 0 0 1 0 2 2 569
LYS 1020 1EF2.subA 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 1 2 566
LYS 1044 1EF2.subA 9 0 2 2 1 0 0 1 1 6 3 0 1 1 0 2 3 0 1 2 566
LYS 1049 1EF2.subA 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 0 2 566
LYS 1083 1EF2.subA 13 2 0 2 0 0 0 0 0 11 2 1 1 2 4 2 3 0 4 2 566
LYS 1089 1EF2.subA 10 3 0 3 0 1 0 1 1 6 4 0 0 0 3 1 2 0 1 2 566
LYS 1098 1EF2.subA 14 1 2 3 0 0 0 0 0 11 3 3 2 5 3 0 3 0 2 2 566
LYS 1124 1EF2.subA 10 0 1 1 0 0 0 0 0 9 1 1 0 1 5 1 2 0 0 2 566
LYS 1196 1EF2.subA 7 0 0 0 0 0 0 0 0 7 0 0 2 2 4 1 0 0 0 2 566
LYS 1382 1EF2.subA 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 0 2 566
LYS 1401 1EF2.subA 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 1 2 566
LYS 1406 1EF2.subA 14 1 0 1 0 1 1 1 2 11 3 3 1 4 2 3 2 0 0 2 566
LYS 1428 1EF2.subA 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 1 2 566
LYS 1443 1EF2.subA 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 566
LYS 1449 1EF2.subA 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 0 2 566
LYS 1515 1EF2.subA 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 3 2 566
LYS 1522 1EF2.subA 17 0 0 0 0 0 1 0 1 16 1 2 3 5 7 1 3 0 3 2 566
LYS 1020 1EJR.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 1 2 553
LYS 1044 1EJR.subC 8 0 1 1 1 0 0 1 1 6 2 0 1 1 0 2 3 0 4 2 553
LYS 1049 1EJR.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 3 2 553
LYS 1083 1EJR.subC 13 2 0 2 0 0 0 0 0 11 2 1 1 2 4 2 3 0 4 2 553
LYS 1089 1EJR.subC 9 2 0 2 0 1 0 1 1 6 3 0 0 0 3 1 2 0 5 2 553
LYS 1098 1EJR.subC 14 1 2 3 0 0 0 0 0 11 3 3 2 5 3 0 3 0 5 2 553
LYS 1124 1EJR.subC 10 0 1 1 0 0 0 0 0 9 1 1 0 1 5 1 2 0 5 2 553
LYS 1196 1EJR.subC 8 0 0 0 0 0 0 0 0 8 0 0 2 2 4 1 0 0 0 2 553
LYS 1382 1EJR.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 0 2 553
LYS 1401 1EJR.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 4 2 553
LYS 1406 1EJR.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 2 3 2 0 0 2 553
LYS 1428 1EJR.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 4 2 553
LYS 1443 1EJR.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 553
LYS 1449 1EJR.subC 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 1 2 553
LYS 1515 1EJR.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 4 2 553
LYS 1522 1EJR.subC 16 0 0 0 0 0 1 0 1 15 1 2 3 5 6 1 3 0 3 2 553
LYS 1020 1EJS.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 3 2 566
LYS 1044 1EJS.subC 8 0 1 1 1 0 0 1 1 6 2 0 1 1 0 2 3 0 4 2 566
LYS 1049 1EJS.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 3 2 566
LYS 1083 1EJS.subC 13 2 0 2 0 0 0 0 0 11 2 1 1 2 4 2 3 0 4 2 566
LYS 1089 1EJS.subC 9 2 0 2 0 1 0 1 1 6 3 0 0 0 3 1 2 0 5 2 566
LYS 1098 1EJS.subC 14 1 2 3 0 0 0 0 0 11 3 3 2 5 3 0 3 0 5 2 566
LYS 1124 1EJS.subC 10 0 1 1 0 0 0 0 0 9 1 1 0 1 5 1 2 0 5 2 566
LYS 1196 1EJS.subC 8 0 0 0 0 0 0 0 0 8 0 0 2 2 4 1 0 0 0 2 566
LYS 1382 1EJS.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 0 2 566
LYS 1401 1EJS.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 3 2 566
LYS 1406 1EJS.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 2 3 2 0 0 2 566
LYS 1428 1EJS.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 4 2 566
LYS 1443 1EJS.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 566
LYS 1449 1EJS.subC 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 1 2 566
LYS 1515 1EJS.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 4 2 566
LYS 1522 1EJS.subC 17 0 0 0 0 0 1 0 1 16 1 2 3 5 6 1 4 0 3 2 566
LYS 1020 1EJT.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 3 2 566
LYS 1044 1EJT.subC 8 0 1 1 1 0 0 1 1 6 2 0 1 1 0 2 3 0 4 2 566
LYS 1049 1EJT.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 3 2 566
LYS 1083 1EJT.subC 13 2 0 2 0 0 0 0 0 11 2 1 1 2 4 2 3 0 4 2 566
LYS 1089 1EJT.subC 9 2 0 2 0 1 0 1 1 6 3 0 0 0 3 1 2 0 5 2 566
LYS 1098 1EJT.subC 14 1 2 3 0 0 0 0 0 11 3 3 2 5 3 0 3 0 5 2 566
LYS 1124 1EJT.subC 10 0 1 1 0 0 0 0 0 9 1 1 0 1 5 1 2 0 4 2 566
LYS 1196 1EJT.subC 7 0 0 0 0 0 0 0 0 7 0 0 2 2 4 1 0 0 0 2 566
LYS 1382 1EJT.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 0 2 566
LYS 1401 1EJT.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 4 2 566
LYS 1406 1EJT.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 2 3 2 0 0 2 566
LYS 1428 1EJT.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 4 2 566
LYS 1443 1EJT.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 566
LYS 1449 1EJT.subC 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 1 2 566
LYS 1515 1EJT.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 4 2 566
LYS 1522 1EJT.subC 17 0 0 0 0 0 1 0 1 16 1 2 3 5 6 1 4 0 3 2 566
LYS 1020 1EJU.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 1 2 553
LYS 1044 1EJU.subC 8 0 1 1 1 0 0 1 1 6 2 0 1 1 0 2 3 0 4 2 553
LYS 1049 1EJU.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 3 2 553
LYS 1083 1EJU.subC 13 2 0 2 0 0 0 0 0 11 2 1 1 2 4 2 3 0 4 2 553
LYS 1089 1EJU.subC 9 2 0 2 0 1 0 1 1 6 3 0 0 0 3 1 2 0 5 2 553
LYS 1098 1EJU.subC 15 1 2 3 0 0 0 0 0 12 3 3 2 5 4 0 3 0 5 2 553
LYS 1124 1EJU.subC 10 0 1 1 0 0 0 0 0 9 1 1 0 1 5 1 2 0 5 2 553
LYS 1196 1EJU.subC 7 0 0 0 0 0 0 0 0 7 0 0 2 2 4 1 0 0 0 2 553
LYS 1382 1EJU.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 0 2 553
LYS 1401 1EJU.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 4 2 553
LYS 1406 1EJU.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 2 3 2 0 0 2 553
LYS 1428 1EJU.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 4 2 553
LYS 1443 1EJU.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 553
LYS 1449 1EJU.subC 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 1 2 553
LYS 1515 1EJU.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 4 2 553
LYS 1522 1EJU.subC 17 0 0 0 0 0 1 0 1 16 1 2 3 5 7 1 3 0 3 2 553
LYS 1020 1EJV.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 1 2 553
LYS 1044 1EJV.subC 8 0 1 1 1 0 0 1 1 6 2 0 1 1 0 2 3 0 4 2 553
LYS 1049 1EJV.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 3 2 553
LYS 1083 1EJV.subC 13 2 0 2 0 0 0 0 0 11 2 1 1 2 4 2 3 0 4 2 553
LYS 1089 1EJV.subC 9 2 0 2 0 1 0 1 1 6 3 0 0 0 3 1 2 0 5 2 553
LYS 1098 1EJV.subC 15 1 2 3 0 0 0 0 0 12 3 3 2 5 4 0 3 0 5 2 553
LYS 1124 1EJV.subC 10 0 1 1 0 0 0 0 0 9 1 1 0 1 5 1 2 0 5 2 553
LYS 1196 1EJV.subC 8 0 0 0 0 0 0 0 0 8 0 0 2 2 4 1 0 0 0 2 553
LYS 1382 1EJV.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 0 2 553
LYS 1401 1EJV.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 4 2 553
LYS 1406 1EJV.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 2 3 2 0 0 2 553
LYS 1428 1EJV.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 4 2 553
LYS 1443 1EJV.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 553
LYS 1449 1EJV.subC 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 1 2 553
LYS 1515 1EJV.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 4 2 553
LYS 1522 1EJV.subC 17 0 0 0 0 0 1 0 1 16 1 2 3 5 6 1 4 0 3 2 553
LYS 1020 1EJW.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 3 2 566
LYS 1044 1EJW.subC 8 0 1 1 1 0 0 1 1 6 2 0 1 1 0 2 3 0 4 2 566
LYS 1049 1EJW.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 3 2 566
LYS 1083 1EJW.subC 13 2 0 2 0 0 0 0 0 11 2 1 1 2 4 2 3 0 4 2 566
LYS 1089 1EJW.subC 9 2 0 2 0 1 0 1 1 6 3 0 0 0 3 1 2 0 5 2 566
LYS 1098 1EJW.subC 14 1 2 3 0 0 0 0 0 11 3 3 2 5 3 0 3 0 5 2 566
LYS 1124 1EJW.subC 10 0 1 1 0 0 0 0 0 9 1 1 0 1 5 1 2 0 5 2 566
LYS 1196 1EJW.subC 8 0 0 0 0 0 0 0 0 8 0 0 2 2 4 1 0 0 0 2 566
LYS 1382 1EJW.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 0 2 566
LYS 1401 1EJW.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 4 2 566
LYS 1406 1EJW.subC 13 1 0 1 0 1 1 1 2 10 3 3 1 4 2 3 1 0 0 2 566
LYS 1428 1EJW.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 4 2 566
LYS 1443 1EJW.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 566
LYS 1449 1EJW.subC 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 1 2 566
LYS 1515 1EJW.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 4 2 566
LYS 1522 1EJW.subC 17 0 0 0 0 0 1 0 1 16 1 2 3 5 6 1 4 0 3 2 566
LYS 1020 1EJX.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 3 2 556
LYS 1044 1EJX.subC 8 0 1 1 1 0 0 1 1 6 2 0 1 1 0 2 3 0 3 2 556
LYS 1049 1EJX.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 3 2 556
LYS 1083 1EJX.subC 14 2 0 2 0 0 0 0 0 12 2 1 1 2 4 2 4 0 4 2 556
LYS 1089 1EJX.subC 10 3 0 3 0 1 0 1 1 6 4 0 0 0 3 1 2 0 6 2 556
LYS 1098 1EJX.subC 14 1 2 3 0 0 0 0 0 11 3 3 2 5 3 0 3 0 5 2 556
LYS 1124 1EJX.subC 10 0 1 1 0 0 0 0 0 9 1 1 0 1 5 1 2 0 6 2 556
LYS 1196 1EJX.subC 8 0 0 0 0 0 0 0 0 8 0 0 2 2 4 1 0 0 0 2 556
LYS 1382 1EJX.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 2 2 556
LYS 1401 1EJX.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 4 2 556
LYS 1406 1EJX.subC 15 1 0 1 0 1 1 1 2 12 3 3 1 4 2 3 3 0 0 2 556
LYS 1428 1EJX.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 3 2 556
LYS 1443 1EJX.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 556
LYS 1449 1EJX.subC 15 1 0 1 0 1 0 1 1 13 2 1 1 2 5 0 6 0 1 2 556
LYS 1515 1EJX.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 5 2 556
LYS 1522 1EJX.subC 18 0 0 0 0 0 1 0 1 17 1 2 3 5 7 1 4 0 3 2 556
LYS 39 1EPV.subA 9 0 1 1 0 0 0 0 0 8 1 0 0 0 3 1 4 0 2 2 382
LYS 76 1EPV.subA 7 1 1 2 0 1 0 1 1 4 3 0 0 0 2 0 2 0 1 2 382
LYS 140 1EPV.subA 4 2 0 2 0 0 0 0 0 2 2 0 0 0 1 0 1 0 4 2 382
LYS 146 1EPV.subA 5 0 2 2 0 1 0 1 1 2 3 1 0 1 0 1 0 0 0 2 382
LYS 234 1EPV.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 3 0 1 0 1 2 382
LYS 242 1EPV.subA 6 1 1 2 0 0 0 0 0 4 2 0 1 1 2 0 1 0 1 2 382
LYS 255 1EPV.subA 5 0 0 0 0 0 1 0 1 4 1 1 0 1 0 0 3 0 2 2 382
LYS 256 1EPV.subA 5 0 2 2 0 0 0 0 0 3 2 0 1 1 0 1 1 0 1 2 382
LYS 262 1EPV.subA 6 0 1 1 0 0 0 0 0 5 1 1 0 1 2 1 1 0 1 2 382
LYS 303 1EPV.subA 6 0 1 1 0 0 1 0 1 4 2 0 1 1 1 0 2 0 0 2 382
LYS 328 1EPV.subA 9 0 1 1 0 1 1 1 2 6 3 3 1 4 0 1 1 0 1 2 382
LYS 372 1EPV.subA 3 1 0 1 0 1 1 1 2 0 3 0 0 0 0 0 0 0 1 2 382
LYS 77 1EYW.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 2 2 330
LYS 82 1EYW.subA 9 1 1 2 0 1 0 1 1 6 3 0 0 0 2 1 3 0 9 2 330
LYS 175 1EYW.subA 3 0 0 0 0 0 0 0 0 3 0 0 1 1 2 0 0 0 3 2 330
LYS 185 1EYW.subA 8 0 2 2 0 1 0 1 1 5 3 0 0 0 1 0 4 0 0 2 330
LYS 285 1EYW.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 3 2 330
LYS 294 1EYW.subA 4 0 0 0 0 1 0 1 1 3 1 0 2 2 0 0 1 0 3 2 330
LYS 339 1EYW.subA 6 0 1 1 1 0 0 1 1 4 2 0 0 0 1 1 2 0 0 2 330
LYS 77 1EZ2.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 0 2 329
LYS 82 1EZ2.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 2 1 2 0 7 2 329
LYS 175 1EZ2.subA 3 0 0 0 0 0 0 0 0 3 0 0 1 1 2 0 0 0 1 2 329
LYS 185 1EZ2.subA 9 0 2 2 0 1 0 1 1 6 3 1 0 1 1 0 4 0 2 2 329
LYS 285 1EZ2.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 1 2 329
LYS 294 1EZ2.subA 4 0 0 0 0 1 0 1 1 3 1 0 2 2 0 0 1 0 2 2 329
LYS 339 1EZ2.subA 6 0 1 1 1 0 0 1 1 4 2 0 0 0 1 1 2 0 0 2 329
LYS 39 1FTX.subA 9 0 1 1 0 0 0 0 0 8 1 0 0 0 3 1 4 0 2 2 380
LYS 76 1FTX.subA 7 1 1 2 0 1 0 1 1 4 3 0 0 0 2 0 2 0 0 2 380
LYS 140 1FTX.subA 5 2 0 2 0 1 0 1 1 2 3 0 0 0 1 0 1 0 5 2 380
LYS 146 1FTX.subA 5 0 2 2 0 1 0 1 1 2 3 1 0 1 0 1 0 0 0 2 380
LYS 234 1FTX.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 3 0 1 0 1 2 380
LYS 242 1FTX.subA 6 1 1 2 0 0 0 0 0 4 2 0 1 1 2 0 1 0 1 2 380
LYS 255 1FTX.subA 5 0 0 0 0 0 1 0 1 4 1 1 0 1 0 0 3 0 2 2 380
LYS 256 1FTX.subA 5 0 2 2 0 0 0 0 0 3 2 0 1 1 0 1 1 0 1 2 380
LYS 262 1FTX.subA 6 0 1 1 0 0 0 0 0 5 1 1 0 1 2 1 1 0 1 2 380
LYS 303 1FTX.subA 6 0 1 1 0 0 1 0 1 4 2 0 1 1 1 0 2 0 0 2 380
LYS 328 1FTX.subA 8 0 0 0 0 1 1 1 2 6 2 3 1 4 0 1 1 0 3 2 380
LYS 372 1FTX.subA 3 1 0 1 0 1 1 1 2 0 3 0 0 0 0 0 0 0 0 2 380
LYS 20 1FWA.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 1 2 566
LYS 44 1FWA.subC 8 0 1 1 1 0 0 1 1 6 2 0 1 1 0 2 3 0 6 2 566
LYS 49 1FWA.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 2 2 566
LYS 83 1FWA.subC 13 2 0 2 0 0 0 0 0 11 2 1 1 2 4 2 3 0 4 2 566
LYS 89 1FWA.subC 10 3 0 3 0 1 0 1 1 6 4 0 0 0 3 1 2 0 2 2 566
LYS 98 1FWA.subC 14 1 2 3 0 0 0 0 0 11 3 3 2 5 3 0 3 0 5 2 566
LYS 124 1FWA.subC 10 0 1 1 0 0 0 0 0 9 1 1 0 1 5 1 2 0 8 2 566
LYS 196 1FWA.subC 8 0 0 0 0 0 0 0 0 8 0 0 2 2 4 1 0 0 0 2 566
LYS 382 1FWA.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 0 2 566
LYS 401 1FWA.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 3 2 566
LYS 406 1FWA.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 2 3 2 0 0 2 566
LYS 428 1FWA.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 2 2 566
LYS 443 1FWA.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 566
LYS 449 1FWA.subC 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 0 2 566
LYS 515 1FWA.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 3 2 566
LYS 522 1FWA.subC 16 0 0 0 0 0 1 0 1 15 1 2 3 5 6 1 3 0 3 2 566
LYS 20 1FWB.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 1 2 566
LYS 44 1FWB.subC 8 0 1 1 1 0 0 1 1 6 2 0 1 1 0 2 3 0 6 2 566
LYS 49 1FWB.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 2 2 566
LYS 83 1FWB.subC 13 2 0 2 0 0 0 0 0 11 2 1 1 2 4 2 3 0 4 2 566
LYS 89 1FWB.subC 10 3 0 3 0 1 0 1 1 6 4 0 0 0 3 1 2 0 2 2 566
LYS 98 1FWB.subC 14 1 2 3 0 0 0 0 0 11 3 3 2 5 3 0 3 0 5 2 566
LYS 124 1FWB.subC 10 0 1 1 0 0 0 0 0 9 1 1 0 1 5 1 2 0 8 2 566
LYS 196 1FWB.subC 7 0 0 0 0 0 0 0 0 7 0 0 2 2 4 1 0 0 0 2 566
LYS 382 1FWB.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 0 2 566
LYS 401 1FWB.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 3 2 566
LYS 406 1FWB.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 2 3 2 0 0 2 566
LYS 428 1FWB.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 2 2 566
LYS 443 1FWB.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 566
LYS 449 1FWB.subC 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 0 2 566
LYS 515 1FWB.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 3 2 566
LYS 522 1FWB.subC 17 0 0 0 0 0 1 0 1 16 1 2 3 5 7 1 3 0 3 2 566
LYS 20 1FWC.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 1 2 566
LYS 44 1FWC.subC 8 0 1 1 1 0 0 1 1 6 2 0 1 1 0 2 3 0 6 2 566
LYS 49 1FWC.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 2 2 566
LYS 83 1FWC.subC 13 2 0 2 0 0 0 0 0 11 2 1 1 2 4 2 3 0 4 2 566
LYS 89 1FWC.subC 10 3 0 3 0 1 0 1 1 6 4 0 0 0 3 1 2 0 2 2 566
LYS 98 1FWC.subC 14 1 2 3 0 0 0 0 0 11 3 3 2 5 3 0 3 0 5 2 566
LYS 124 1FWC.subC 10 0 1 1 0 0 0 0 0 9 1 1 0 1 5 1 2 0 8 2 566
LYS 196 1FWC.subC 8 0 0 0 0 0 0 0 0 8 0 0 2 2 4 1 0 0 0 2 566
LYS 382 1FWC.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 0 2 566
LYS 401 1FWC.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 3 2 566
LYS 406 1FWC.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 2 3 2 0 0 2 566
LYS 428 1FWC.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 2 2 566
LYS 443 1FWC.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 566
LYS 449 1FWC.subC 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 0 2 566
LYS 515 1FWC.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 3 2 566
LYS 522 1FWC.subC 16 0 0 0 0 0 1 0 1 15 1 2 3 5 6 1 3 0 3 2 566
LYS 20 1FWD.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 1 2 566
LYS 44 1FWD.subC 8 0 1 1 1 0 0 1 1 6 2 0 1 1 0 2 3 0 6 2 566
LYS 49 1FWD.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 3 2 566
LYS 83 1FWD.subC 14 2 0 2 0 0 0 0 0 12 2 1 1 2 4 2 4 0 4 2 566
LYS 89 1FWD.subC 10 3 0 3 0 1 0 1 1 6 4 0 0 0 3 1 2 0 2 2 566
LYS 98 1FWD.subC 14 1 2 3 0 0 0 0 0 11 3 3 2 5 3 0 3 0 5 2 566
LYS 124 1FWD.subC 10 0 1 1 0 0 0 0 0 9 1 1 0 1 5 1 2 0 8 2 566
LYS 196 1FWD.subC 8 0 0 0 0 0 0 0 0 8 0 0 2 2 4 1 0 0 0 2 566
LYS 382 1FWD.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 0 2 566
LYS 401 1FWD.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 3 2 566
LYS 406 1FWD.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 2 3 2 0 0 2 566
LYS 428 1FWD.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 2 2 566
LYS 443 1FWD.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 566
LYS 449 1FWD.subC 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 0 2 566
LYS 515 1FWD.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 3 2 566
LYS 522 1FWD.subC 16 0 0 0 0 0 1 0 1 15 1 2 3 5 6 1 3 0 3 2 566
LYS 20 1FWE.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 1 2 552
LYS 44 1FWE.subC 8 0 1 1 1 0 0 1 1 6 2 0 1 1 0 2 3 0 6 2 552
LYS 49 1FWE.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 1 2 552
LYS 83 1FWE.subC 14 2 0 2 0 0 0 0 0 12 2 1 1 2 4 2 4 0 4 2 552
LYS 89 1FWE.subC 10 3 0 3 0 1 0 1 1 6 4 0 0 0 3 1 2 0 2 2 552
LYS 98 1FWE.subC 14 1 2 3 0 0 0 0 0 11 3 3 2 5 3 0 3 0 4 2 552
LYS 124 1FWE.subC 10 0 1 1 0 0 0 0 0 9 1 1 0 1 5 1 2 0 7 2 552
LYS 196 1FWE.subC 8 0 0 0 0 0 0 0 0 8 0 0 2 2 4 1 0 0 0 2 552
LYS 382 1FWE.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 0 2 552
LYS 401 1FWE.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 4 2 552
LYS 406 1FWE.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 2 3 2 0 0 2 552
LYS 428 1FWE.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 1 2 552
LYS 443 1FWE.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 552
LYS 449 1FWE.subC 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 0 2 552
LYS 515 1FWE.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 4 2 552
LYS 522 1FWE.subC 17 0 0 0 0 0 1 0 1 16 1 2 3 5 7 1 3 0 3 2 552
LYS 20 1FWF.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 1 2 551
LYS 44 1FWF.subC 8 0 1 1 1 0 0 1 1 6 2 0 1 1 0 2 3 0 6 2 551
LYS 49 1FWF.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 3 2 551
LYS 83 1FWF.subC 14 2 0 2 0 0 0 0 0 12 2 1 1 2 4 2 4 0 4 2 551
LYS 89 1FWF.subC 10 3 0 3 0 1 0 1 1 6 4 0 0 0 3 1 2 0 2 2 551
LYS 98 1FWF.subC 14 1 2 3 0 0 0 0 0 11 3 3 2 5 3 0 3 0 5 2 551
LYS 124 1FWF.subC 10 0 1 1 0 0 0 0 0 9 1 1 0 1 5 1 2 0 7 2 551
LYS 196 1FWF.subC 8 0 0 0 0 0 0 0 0 8 0 0 2 2 4 1 0 0 0 2 551
LYS 382 1FWF.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 0 2 551
LYS 401 1FWF.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 3 2 551
LYS 406 1FWF.subC 13 1 0 1 0 1 1 1 2 10 3 3 1 4 2 3 1 0 0 2 551
LYS 428 1FWF.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 2 2 551
LYS 443 1FWF.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 551
LYS 449 1FWF.subC 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 0 2 551
LYS 515 1FWF.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 3 2 551
LYS 522 1FWF.subC 17 0 0 0 0 0 1 0 1 16 1 2 3 5 7 1 3 0 3 2 551
LYS 20 1FWG.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 1 2 566
LYS 44 1FWG.subC 8 0 1 1 1 0 0 1 1 6 2 0 1 1 0 2 3 0 6 2 566
LYS 49 1FWG.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 2 2 566
LYS 83 1FWG.subC 14 2 0 2 0 0 0 0 0 12 2 1 1 2 4 2 4 0 4 2 566
LYS 89 1FWG.subC 10 3 0 3 0 1 0 1 1 6 4 0 0 0 3 1 2 0 2 2 566
LYS 98 1FWG.subC 14 1 2 3 0 0 0 0 0 11 3 3 2 5 3 0 3 0 5 2 566
LYS 124 1FWG.subC 10 0 1 1 0 0 0 0 0 9 1 1 0 1 5 1 2 0 7 2 566
LYS 196 1FWG.subC 8 0 0 0 0 0 0 0 0 8 0 0 2 2 4 1 0 0 0 2 566
LYS 382 1FWG.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 0 2 566
LYS 401 1FWG.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 3 2 566
LYS 406 1FWG.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 2 3 2 0 0 2 566
LYS 428 1FWG.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 2 2 566
LYS 443 1FWG.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 566
LYS 449 1FWG.subC 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 0 2 566
LYS 515 1FWG.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 3 2 566
LYS 522 1FWG.subC 17 0 0 0 0 0 1 0 1 16 1 2 3 5 7 1 3 0 3 2 566
LYS 20 1FWH.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 1 2 566
LYS 44 1FWH.subC 8 0 1 1 1 0 0 1 1 6 2 0 1 1 0 2 3 0 6 2 566
LYS 49 1FWH.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 3 2 566
LYS 83 1FWH.subC 13 2 0 2 0 0 0 0 0 11 2 1 1 2 4 2 3 0 4 2 566
LYS 89 1FWH.subC 10 3 0 3 0 1 0 1 1 6 4 0 0 0 3 1 2 0 2 2 566
LYS 98 1FWH.subC 14 1 2 3 0 0 0 0 0 11 3 3 2 5 3 0 3 0 5 2 566
LYS 124 1FWH.subC 9 0 1 1 0 0 0 0 0 8 1 1 0 1 4 1 2 0 8 2 566
LYS 196 1FWH.subC 8 0 0 0 0 0 0 0 0 8 0 0 2 2 4 1 0 0 0 2 566
LYS 382 1FWH.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 0 2 566
LYS 401 1FWH.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 3 2 566
LYS 406 1FWH.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 2 3 2 0 0 2 566
LYS 428 1FWH.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 2 2 566
LYS 443 1FWH.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 566
LYS 449 1FWH.subC 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 0 2 566
LYS 515 1FWH.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 3 2 566
LYS 522 1FWH.subC 16 0 0 0 0 0 1 0 1 15 1 2 3 5 6 1 3 0 3 2 566
LYS 20 1FWI.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 3 2 551
LYS 44 1FWI.subC 8 0 1 1 1 0 0 1 1 6 2 0 1 1 0 2 3 0 3 2 551
LYS 49 1FWI.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 2 2 551
LYS 83 1FWI.subC 13 2 0 2 0 0 0 0 0 11 2 1 1 2 4 2 3 0 4 2 551
LYS 89 1FWI.subC 10 3 0 3 0 1 0 1 1 6 4 0 0 0 3 1 2 0 3 2 551
LYS 98 1FWI.subC 15 1 2 3 0 0 0 0 0 12 3 3 2 5 4 0 3 0 5 2 551
LYS 124 1FWI.subC 10 0 1 1 0 0 0 0 0 9 1 1 0 1 5 1 2 0 8 2 551
LYS 196 1FWI.subC 7 0 0 0 0 0 0 0 0 7 0 0 2 2 4 1 0 0 0 2 551
LYS 382 1FWI.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 1 2 551
LYS 401 1FWI.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 2 2 551
LYS 406 1FWI.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 2 3 2 0 0 2 551
LYS 428 1FWI.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 3 2 551
LYS 443 1FWI.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 551
LYS 449 1FWI.subC 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 1 2 551
LYS 515 1FWI.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 3 2 551
LYS 522 1FWI.subC 17 0 0 0 0 0 1 0 1 16 1 2 3 5 7 1 3 0 3 2 551
LYS 20 1FWJ.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 1 2 566
LYS 44 1FWJ.subC 9 0 2 2 1 0 0 1 1 6 3 0 1 1 0 2 3 0 2 2 566
LYS 49 1FWJ.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 1 2 566
LYS 83 1FWJ.subC 13 2 0 2 0 0 0 0 0 11 2 1 1 2 4 2 3 0 4 2 566
LYS 89 1FWJ.subC 10 3 0 3 0 1 0 1 1 6 4 0 0 0 3 1 2 0 1 2 566
LYS 98 1FWJ.subC 15 1 2 3 0 0 0 0 0 12 3 3 2 5 4 0 3 0 2 2 566
LYS 124 1FWJ.subC 10 0 1 1 0 0 0 0 0 9 1 1 0 1 5 1 2 0 1 2 566
LYS 196 1FWJ.subC 7 0 0 0 0 0 0 0 0 7 0 0 2 2 4 1 0 0 0 2 566
LYS 382 1FWJ.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 0 2 566
LYS 401 1FWJ.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 1 2 566
LYS 406 1FWJ.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 2 3 2 0 0 2 566
LYS 428 1FWJ.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 1 2 566
LYS 443 1FWJ.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 566
LYS 449 1FWJ.subC 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 0 2 566
LYS 515 1FWJ.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 3 2 566
LYS 522 1FWJ.subC 17 0 0 0 0 0 1 0 1 16 1 2 3 5 7 1 3 0 3 2 566
LYS 8 1GK8.subA 2 0 0 0 0 0 0 0 0 2 0 1 0 1 1 0 0 0 0 2 469
LYS 14 1GK8.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 7 2 469
LYS 18 1GK8.subA 6 1 0 1 0 0 0 0 0 5 1 2 0 2 1 0 2 0 12 2 469
LYS 81 1GK8.subA 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 6 2 469
LYS 128 1GK8.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 0 2 469
LYS 146 1GK8.subA 5 0 0 0 0 0 0 0 0 4 0 1 0 1 2 0 1 0 2 2 469
LYS 161 1GK8.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 1 0 2 0 2 2 469
LYS 164 1GK8.subA 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 4 2 469
LYS 175 1GK8.subA 8 1 0 1 1 0 0 1 1 6 2 1 0 1 3 0 2 1 0 2 469
LYS 177 1GK8.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 2 2 469
LYS 183 1GK8.subA 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 3 2 469
LYS 227 1GK8.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 2 0 0 3 2 469
LYS 236 1GK8.subA 12 1 1 2 0 1 0 1 1 9 3 1 1 2 3 0 4 0 2 2 469
LYS 252 1GK8.subA 5 0 2 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 11 2 469
LYS 258 1GK8.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 2 1 2 0 6 2 469
LYS 316 1GK8.subA 13 1 0 1 0 2 0 2 2 10 3 0 1 1 2 2 5 0 1 2 469
LYS 334 1GK8.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 2 2 469
LYS 356 1GK8.subA 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 6 2 469
LYS 450 1GK8.subA 4 0 0 0 0 1 0 1 1 3 1 1 0 1 1 1 0 0 4 2 469
LYS 463 1GK8.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 2 0 0 7 2 469
LYS 466 1GK8.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 7 2 469
LYS 474 1GK8.subA 6 2 0 2 0 1 0 1 1 3 3 1 0 1 0 0 2 0 9 2 469
LYS 6 1GKP.subA 8 2 0 2 0 0 0 0 0 6 2 0 1 1 1 0 4 0 4 2 458
LYS 18 1GKP.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 2 1 0 0 11 2 458
LYS 49 1GKP.subA 9 1 0 1 0 0 0 0 0 8 1 2 0 2 3 2 1 0 6 2 458
LYS 72 1GKP.subA 9 2 0 2 0 0 0 0 0 7 2 0 3 3 2 2 0 0 9 2 458
LYS 80 1GKP.subA 8 0 1 1 0 0 1 0 1 6 2 2 0 2 1 2 1 0 4 2 458
LYS 110 1GKP.subA 8 1 0 1 0 1 0 1 1 6 2 1 1 2 0 3 1 0 8 2 458
LYS 112 1GKP.subA 6 0 0 0 0 0 1 0 1 5 1 1 0 1 1 2 1 0 5 2 458
LYS 129 1GKP.subA 10 3 1 4 0 1 0 1 1 5 5 3 1 4 0 1 0 0 10 2 458
LYS 133 1GKP.subA 5 2 1 3 0 0 0 0 0 2 3 1 1 2 0 0 0 0 7 2 458
LYS 156 1GKP.subA 5 0 1 1 0 1 0 1 1 3 2 0 1 1 1 1 0 0 9 2 458
LYS 174 1GKP.subA 9 0 1 1 0 1 0 1 1 7 2 1 0 1 3 2 1 0 8 2 458
LYS 196 1GKP.subA 6 0 1 1 0 1 0 1 1 4 2 1 1 2 0 0 2 0 8 2 458
LYS 202 1GKP.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 1 1 1 0 4 2 458
LYS 243 1GKP.subA 5 1 0 1 0 0 0 0 0 4 1 0 1 1 2 0 1 0 10 2 458
LYS 253 1GKP.subA 12 1 0 1 0 0 0 0 0 11 1 0 0 0 4 2 5 0 2 2 458
LYS 272 1GKP.subA 7 2 2 4 0 1 0 1 1 2 5 1 0 1 0 0 1 0 11 2 458
LYS 284 1GKP.subA 7 0 1 1 0 0 0 0 0 6 1 1 0 1 1 2 2 0 6 2 458
LYS 294 1GKP.subA 4 2 0 2 0 1 0 1 1 1 3 0 1 1 0 0 0 0 8 2 458
LYS 298 1GKP.subA 4 1 0 1 0 1 0 1 1 2 2 0 1 1 0 0 1 0 11 2 458
LYS 324 1GKP.subA 14 1 0 1 0 0 0 0 0 13 1 1 2 3 2 2 6 0 1 2 458
LYS 328 1GKP.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 1 0 3 0 6 2 458
LYS 370 1GKP.subA 9 2 0 2 1 1 0 2 2 5 4 2 0 2 1 1 1 0 11 2 458
LYS 373 1GKP.subA 6 1 0 1 1 0 0 1 1 4 2 1 0 1 1 1 1 0 6 2 458
LYS 381 1GKP.subA 12 3 0 3 0 2 0 2 2 7 5 1 0 1 2 1 3 0 6 2 458
LYS 406 1GKP.subA 4 0 1 1 0 0 0 0 0 3 1 2 0 2 0 0 1 0 2 2 458
LYS 434 1GKP.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 1 1 2 0 6 2 458
LYS 446 1GKP.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 1 1 0 0 5 2 458
LYS 450 1GKP.subA 6 1 0 1 0 1 0 1 1 4 2 0 0 0 2 1 1 0 13 2 458
LYS 6 1GKQ.subA 7 2 0 2 0 0 0 0 0 5 2 0 1 1 0 0 4 0 0 2 458
LYS 18 1GKQ.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 2 1 0 0 0 2 458
LYS 49 1GKQ.subA 9 1 0 1 0 0 0 0 0 8 1 2 0 2 3 2 1 0 1 2 458
LYS 72 1GKQ.subA 9 2 0 2 0 0 0 0 0 7 2 0 3 3 2 2 0 0 1 2 458
LYS 80 1GKQ.subA 8 0 1 1 0 0 1 0 1 6 2 2 0 2 1 2 1 0 0 2 458
LYS 110 1GKQ.subA 9 1 1 2 0 1 0 1 1 6 3 1 1 2 0 3 1 0 3 2 458
LYS 112 1GKQ.subA 6 0 0 0 0 0 1 0 1 5 1 1 0 1 1 2 1 0 1 2 458
LYS 129 1GKQ.subA 10 3 1 4 0 1 0 1 1 5 5 3 1 4 0 1 0 0 3 2 458
LYS 133 1GKQ.subA 5 2 1 3 0 0 0 0 0 2 3 1 1 2 0 0 0 0 1 2 458
LYS 156 1GKQ.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 1 1 0 0 3 2 458
LYS 174 1GKQ.subA 8 0 1 1 0 1 0 1 1 6 2 1 0 1 3 1 1 0 4 2 458
LYS 196 1GKQ.subA 6 0 1 1 0 1 0 1 1 4 2 1 1 2 0 0 2 0 1 2 458
LYS 202 1GKQ.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 1 1 1 0 0 2 458
LYS 243 1GKQ.subA 5 1 0 1 0 0 0 0 0 4 1 0 1 1 2 0 1 0 2 2 458
LYS 253 1GKQ.subA 12 1 0 1 0 0 0 0 0 11 1 0 0 0 4 2 5 0 3 2 458
LYS 272 1GKQ.subA 7 2 2 4 0 1 0 1 1 2 5 1 0 1 0 0 1 0 6 2 458
LYS 284 1GKQ.subA 7 0 1 1 0 0 0 0 0 6 1 1 0 1 1 2 2 0 1 2 458
LYS 294 1GKQ.subA 4 2 0 2 0 1 0 1 1 1 3 0 1 1 0 0 0 0 3 2 458
LYS 298 1GKQ.subA 4 1 0 1 0 1 0 1 1 2 2 0 1 1 0 0 1 0 2 2 458
LYS 324 1GKQ.subA 13 1 0 1 0 0 0 0 0 12 1 1 2 3 2 2 5 0 1 2 458
LYS 328 1GKQ.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 1 2 458
LYS 370 1GKQ.subA 9 2 0 2 1 0 0 1 1 6 3 2 0 2 2 1 1 0 3 2 458
LYS 373 1GKQ.subA 6 1 0 1 1 0 0 1 1 4 2 1 0 1 1 1 1 0 1 2 458
LYS 381 1GKQ.subA 11 3 0 3 0 2 0 2 2 6 5 1 0 1 2 0 3 0 2 2 458
LYS 406 1GKQ.subA 4 0 1 1 0 0 0 0 0 3 1 2 0 2 0 0 1 0 0 2 458
LYS 434 1GKQ.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 1 1 2 0 2 2 458
LYS 446 1GKQ.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 2 1 0 0 0 2 458
LYS 450 1GKQ.subA 6 1 0 1 0 1 0 1 1 4 2 0 0 0 2 1 1 0 0 2 458
LYS 7 1GKR.subA 8 2 0 2 0 0 0 0 0 6 2 2 1 3 0 0 3 0 0 2 451
LYS 25 1GKR.subA 6 2 0 2 0 0 0 0 0 4 2 0 0 0 1 1 2 0 1 2 451
LYS 28 1GKR.subA 7 2 0 2 0 0 0 0 0 5 2 1 0 1 3 0 1 0 3 2 451
LYS 50 1GKR.subA 7 1 0 1 0 0 0 0 0 6 1 0 0 0 3 1 2 0 1 2 451
LYS 69 1GKR.subA 9 2 0 2 0 2 0 2 2 5 4 0 2 2 1 0 2 0 2 2 451
LYS 109 1GKR.subA 14 0 2 2 0 0 0 0 0 12 2 2 0 2 3 2 5 0 3 2 451
LYS 110 1GKR.subA 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 3 2 0 0 2 451
LYS 111 1GKR.subA 4 0 1 1 0 0 0 0 0 3 1 0 2 2 0 0 1 0 2 2 451
LYS 118 1GKR.subA 11 0 3 3 0 1 0 1 1 7 4 1 1 2 1 0 4 0 1 2 451
LYS 137 1GKR.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 1 0 3 0 0 2 451
LYS 195 1GKR.subA 4 0 0 0 0 0 0 0 0 4 0 0 3 3 1 0 0 0 0 2 451
LYS 198 1GKR.subA 5 0 0 0 0 0 0 0 0 5 0 0 1 1 2 1 1 0 0 2 451
LYS 203 1GKR.subA 6 2 0 2 0 0 0 0 0 4 2 0 0 0 2 1 1 0 0 2 451
LYS 229 1GKR.subA 6 0 1 1 0 0 0 0 0 5 1 0 2 2 0 0 3 0 5 2 451
LYS 283 1GKR.subA 12 0 0 0 1 0 0 1 1 11 1 0 1 1 4 2 4 0 1 2 451
LYS 321 1GKR.subA 14 1 1 2 1 0 1 1 2 10 4 0 1 1 5 2 2 0 1 2 451
LYS 326 1GKR.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 1 0 0 0 2 451
LYS 330 1GKR.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 2 1 0 0 0 2 451
LYS 352 1GKR.subA 6 0 0 0 0 1 0 1 1 5 1 1 2 3 2 0 0 0 1 2 451
LYS 367 1GKR.subA 7 1 2 3 0 0 0 0 0 4 3 1 0 1 1 0 2 0 2 2 451
LYS 370 1GKR.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 1 1 1 0 1 2 451
LYS 378 1GKR.subA 12 3 0 3 0 1 0 1 1 8 4 1 1 2 2 1 3 0 2 2 451
LYS 399 1GKR.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 2 0 2 0 0 2 451
LYS 410 1GKR.subA 9 0 0 0 0 1 1 1 2 7 2 2 1 3 1 2 1 0 1 2 451
LYS 436 1GKR.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 0 2 451
LYS 6 1H01.subA 5 0 1 1 1 0 0 1 1 3 2 0 1 1 0 1 1 0 1 2 287
LYS 9 1H01.subA 7 0 2 2 0 0 0 0 0 5 2 0 0 0 1 1 3 0 2 2 287
LYS 20 1H01.subA 9 0 1 1 0 0 1 0 1 7 2 0 0 0 1 2 4 0 2 2 287
LYS 24 1H01.subA 7 0 1 1 0 1 0 1 1 5 2 0 3 3 0 1 1 0 0 2 287
LYS 34 1H01.subA 8 0 0 0 2 1 0 3 3 5 3 0 0 0 0 2 2 0 1 2 287
LYS 56 1H01.subA 9 1 1 2 0 0 0 0 0 7 2 1 0 1 0 0 6 0 3 2 287
LYS 65 1H01.subA 6 0 1 1 0 0 0 0 0 5 1 0 1 1 0 0 4 0 5 2 287
LYS 75 1H01.subA 8 0 1 1 1 1 0 2 2 5 3 1 1 2 0 1 2 0 0 2 287
LYS 88 1H01.subA 7 2 1 3 0 0 0 0 0 4 3 0 1 1 1 0 2 0 3 2 287
LYS 89 1H01.subA 3 1 0 1 0 0 0 0 0 2 1 0 1 1 0 1 0 0 3 2 287
LYS 105 1H01.subA 8 1 0 1 0 0 0 0 0 7 1 1 0 1 1 2 3 0 6 2 287
LYS 129 1H01.subA 9 1 0 1 0 0 0 0 0 8 1 2 2 4 1 1 2 0 6 2 287
LYS 142 1H01.subA 11 0 1 1 0 0 0 0 0 10 1 0 2 2 2 0 6 0 4 2 287
LYS 178 1H01.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 0 2 287
LYS 237 1H01.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 2 2 0 0 2 2 287
LYS 242 1H01.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 2 2 0 0 6 2 287
LYS 250 1H01.subA 5 1 0 1 0 0 0 0 0 4 1 1 1 2 0 1 1 0 2 2 287
LYS 273 1H01.subA 7 1 0 1 0 1 1 1 2 4 3 0 1 1 1 1 1 0 4 2 287
LYS 278 1H01.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 3 0 0 0 5 2 287
LYS 291 1H01.subA 7 0 1 1 0 0 0 0 0 6 1 2 1 3 3 0 0 0 8 2 287
LYS 9 1HL8.subA 5 1 0 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 0 2 426
LYS 22 1HL8.subA 3 0 0 0 1 0 0 1 1 2 1 0 0 0 1 1 0 0 1 2 426
LYS 26 1HL8.subA 5 1 0 1 1 0 0 1 1 3 2 0 0 0 2 1 0 0 1 2 426
LYS 28 1HL8.subA 13 1 0 1 1 0 0 1 1 11 2 2 1 3 3 2 3 0 3 2 426
LYS 75 1HL8.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 0 2 426
LYS 85 1HL8.subA 4 0 1 1 0 0 0 0 0 3 1 1 0 1 0 1 1 0 0 2 426
LYS 95 1HL8.subA 6 1 2 3 0 1 0 1 1 2 4 0 0 0 0 2 0 0 0 2 426
LYS 104 1HL8.subA 4 1 1 2 0 0 0 0 0 2 2 1 0 1 0 1 0 0 0 2 426
LYS 115 1HL8.subA 5 1 0 1 0 0 0 0 0 4 1 0 0 0 2 1 1 0 0 2 426
LYS 116 1HL8.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 0 2 426
LYS 120 1HL8.subA 6 0 0 0 0 0 0 0 0 6 0 0 0 0 3 3 0 0 2 2 426
LYS 127 1HL8.subA 13 1 0 1 0 0 1 0 1 11 2 2 1 3 2 5 1 0 2 2 426
LYS 138 1HL8.subA 4 0 1 1 0 0 0 0 0 3 1 1 0 1 1 1 0 0 3 2 426
LYS 146 1HL8.subA 4 1 0 1 0 1 0 1 1 2 2 0 1 1 0 0 1 0 1 2 426
LYS 150 1HL8.subA 6 0 0 0 0 2 0 2 2 4 2 1 0 1 2 0 1 0 0 2 426
LYS 159 1HL8.subA 5 1 1 2 0 0 0 0 0 3 2 0 0 0 3 0 0 0 3 2 426
LYS 207 1HL8.subA 5 1 0 1 0 0 0 0 0 4 1 0 1 1 0 3 0 0 2 2 426
LYS 230 1HL8.subA 6 0 1 1 0 0 0 0 0 5 1 1 1 2 2 1 0 0 2 2 426
LYS 232 1HL8.subA 8 0 2 2 0 0 0 0 0 6 2 0 0 0 3 2 1 0 5 2 426
LYS 236 1HL8.subA 8 2 1 3 0 0 0 0 0 5 3 0 0 0 1 2 2 0 2 2 426
LYS 245 1HL8.subA 7 1 0 1 0 0 1 0 1 5 2 0 1 1 0 2 2 0 1 2 426
LYS 263 1HL8.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 0 4 0 0 0 2 426
LYS 279 1HL8.subA 11 2 0 2 1 0 0 1 1 8 3 0 1 1 2 4 1 0 3 2 426
LYS 317 1HL8.subA 10 1 0 1 0 1 0 1 1 8 2 1 0 1 4 1 2 0 0 2 426
LYS 328 1HL8.subA 7 1 0 1 0 0 0 0 0 6 1 1 0 1 4 1 0 0 1 2 426
LYS 350 1HL8.subA 8 1 1 2 0 1 0 1 1 5 3 2 0 2 0 2 1 0 0 2 426
LYS 367 1HL8.subA 6 0 2 2 0 0 0 0 0 4 2 2 0 2 2 0 0 0 0 2 426
LYS 379 1HL8.subA 10 1 0 1 0 2 0 2 2 7 3 1 0 1 4 2 0 0 0 2 426
LYS 395 1HL8.subA 5 0 2 2 0 0 0 0 0 3 2 1 0 1 0 0 2 0 1 2 426
LYS 420 1HL8.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 0 1 1 0 1 2 426
LYS 424 1HL8.subA 3 1 0 1 0 0 0 0 0 2 1 0 1 1 1 0 0 0 0 2 426
LYS 432 1HL8.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 1 0 1 0 0 2 426
LYS 433 1HL8.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 0 1 0 0 2 426
LYS 9 1HL9.subA 5 1 0 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 0 2 421
LYS 22 1HL9.subA 4 0 0 0 1 0 0 1 1 3 1 0 0 0 1 1 1 0 2 2 421
LYS 26 1HL9.subA 5 1 0 1 1 0 0 1 1 3 2 0 0 0 2 1 0 0 4 2 421
LYS 28 1HL9.subA 13 1 0 1 1 0 0 1 1 11 2 2 1 3 3 2 3 0 5 2 421
LYS 75 1HL9.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 0 1 1 0 0 2 421
LYS 85 1HL9.subA 4 0 1 1 0 0 0 0 0 3 1 1 0 1 0 1 1 0 0 2 421
LYS 95 1HL9.subA 7 1 2 3 0 1 0 1 1 3 4 0 0 0 0 3 0 0 1 2 421
LYS 104 1HL9.subA 4 1 1 2 0 0 0 0 0 2 2 1 0 1 0 1 0 0 0 2 421
LYS 115 1HL9.subA 5 1 0 1 0 0 0 0 0 4 1 0 0 0 2 1 1 0 1 2 421
LYS 116 1HL9.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 0 2 421
LYS 120 1HL9.subA 6 0 0 0 0 0 0 0 0 6 0 0 0 0 3 3 0 0 3 2 421
LYS 127 1HL9.subA 14 1 0 1 0 0 1 0 1 12 2 2 1 3 3 5 1 0 2 2 421
LYS 138 1HL9.subA 3 0 1 1 0 0 0 0 0 2 1 1 0 1 0 1 0 0 1 2 421
LYS 146 1HL9.subA 4 1 0 1 0 1 0 1 1 2 2 0 1 1 0 0 1 0 1 2 421
LYS 150 1HL9.subA 6 0 0 0 0 2 0 2 2 4 2 1 0 1 2 0 1 0 0 2 421
LYS 159 1HL9.subA 5 1 1 2 0 0 0 0 0 3 2 0 0 0 3 0 0 0 0 2 421
LYS 207 1HL9.subA 5 1 0 1 0 0 0 0 0 4 1 0 1 1 0 3 0 0 4 2 421
LYS 230 1HL9.subA 6 0 1 1 0 0 0 0 0 5 1 1 1 2 2 1 0 0 3 2 421
LYS 232 1HL9.subA 9 0 2 2 0 0 0 0 0 7 2 0 0 0 3 2 2 0 2 2 421
LYS 236 1HL9.subA 7 1 1 2 0 0 0 0 0 5 2 0 0 0 1 2 2 0 1 2 421
LYS 245 1HL9.subA 7 1 0 1 0 0 1 0 1 5 2 0 1 1 0 2 2 0 1 2 421
LYS 263 1HL9.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 0 4 0 0 1 2 421
LYS 279 1HL9.subA 11 2 0 2 1 0 0 1 1 8 3 0 1 1 2 4 1 0 5 2 421
LYS 317 1HL9.subA 11 1 0 1 0 1 0 1 1 9 2 1 0 1 5 1 2 0 1 2 421
LYS 328 1HL9.subA 7 1 0 1 0 0 0 0 0 6 1 1 0 1 4 1 0 0 0 2 421
LYS 350 1HL9.subA 8 1 1 2 0 1 0 1 1 5 3 2 0 2 0 2 1 0 3 2 421
LYS 367 1HL9.subA 6 0 2 2 0 0 0 0 0 4 2 2 0 2 2 0 0 0 1 2 421
LYS 379 1HL9.subA 10 1 0 1 0 2 0 2 2 7 3 1 0 1 4 2 0 0 0 2 421
LYS 395 1HL9.subA 4 0 1 1 0 0 0 0 0 3 1 1 0 1 0 0 2 0 0 2 421
LYS 420 1HL9.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 0 1 1 0 1 2 421
LYS 424 1HL9.subA 4 1 1 2 0 0 0 0 0 2 2 0 1 1 1 0 0 0 1 2 421
LYS 432 1HL9.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 1 0 1 0 0 2 421
LYS 433 1HL9.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 0 1 0 0 2 421
LYS 2 1IE7.subC 2 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 2 0 1 2 570
LYS 33 1IE7.subC 4 1 1 2 0 0 0 0 0 2 2 0 0 0 0 2 0 0 13 2 570
LYS 48 1IE7.subC 7 0 0 0 0 0 0 0 0 7 0 0 1 1 3 1 2 0 7 2 570
LYS 84 1IE7.subC 12 1 0 1 0 0 0 0 0 11 1 1 1 2 4 3 2 0 4 2 570
LYS 90 1IE7.subC 11 2 1 3 0 1 0 1 1 7 4 1 1 2 0 1 4 0 9 2 570
LYS 99 1IE7.subC 14 3 1 4 0 0 0 0 0 10 4 3 0 3 3 1 3 0 6 2 570
LYS 127 1IE7.subC 11 1 2 3 0 1 0 1 1 7 4 0 0 0 3 1 3 0 5 2 570
LYS 169 1IE7.subC 7 0 1 1 0 0 0 0 0 6 1 1 0 1 4 0 1 0 7 2 570
LYS 182 1IE7.subC 5 0 1 1 1 0 0 1 1 3 2 0 1 1 0 1 1 0 3 2 570
LYS 185 1IE7.subC 6 0 1 1 1 0 0 1 1 4 2 1 0 1 0 0 3 0 3 2 570
LYS 199 1IE7.subC 8 0 0 0 0 0 1 0 1 7 1 0 1 1 3 1 1 0 0 2 570
LYS 326 1IE7.subC 6 1 0 1 0 0 1 0 1 4 2 0 2 2 0 0 2 0 2 2 570
LYS 383 1IE7.subC 12 2 1 3 0 0 0 0 0 9 3 1 2 3 3 1 2 0 2 2 570
LYS 385 1IE7.subC 10 3 0 3 0 1 0 1 1 6 4 0 0 0 3 1 2 0 11 2 570
LYS 386 1IE7.subC 5 2 1 3 0 0 0 0 0 2 3 0 1 1 0 0 1 0 9 2 570
LYS 395 1IE7.subC 7 0 3 3 1 0 0 1 1 3 4 1 2 3 0 0 0 0 10 2 570
LYS 404 1IE7.subC 12 1 1 2 1 1 0 2 2 8 4 0 1 1 1 3 3 0 6 2 570
LYS 409 1IE7.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 1 2 4 0 0 2 570
LYS 431 1IE7.subC 12 0 3 3 0 0 1 0 1 8 4 1 0 1 3 1 3 0 6 2 570
LYS 441 1IE7.subC 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 1 0 0 4 2 570
LYS 446 1IE7.subC 7 0 0 0 0 0 0 0 0 7 0 0 1 1 2 3 1 0 3 2 570
LYS 452 1IE7.subC 14 1 0 1 0 0 0 0 0 13 1 2 2 4 3 1 5 0 3 2 570
LYS 497 1IE7.subC 8 1 0 1 1 0 0 1 1 6 2 3 1 4 1 0 1 0 11 2 570
LYS 507 1IE7.subC 5 0 1 1 0 0 0 0 0 4 1 0 0 0 2 0 2 0 10 2 570
LYS 511 1IE7.subC 5 0 1 1 0 1 1 1 2 2 3 0 0 0 1 0 1 0 9 2 570
LYS 518 1IE7.subC 7 1 0 1 1 0 1 1 2 4 3 1 1 2 1 0 1 0 9 2 570
LYS 525 1IE7.subC 16 1 0 1 0 0 1 0 1 14 2 2 1 3 6 0 5 0 2 2 570
LYS 526 1IE7.subC 2 1 0 1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 4 2 570
LYS 529 1IE7.subC 8 1 0 1 0 0 0 0 0 7 1 1 1 2 0 2 3 0 11 2 570
LYS 547 1IE7.subC 7 1 2 3 0 0 0 0 0 4 3 0 0 0 1 0 3 0 3 2 570
LYS 559 1IE7.subC 3 0 1 1 0 0 0 0 0 2 1 0 0 0 0 0 2 0 5 2 570
LYS 14 1IR1.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 1 0 0 6 2 464
LYS 18 1IR1.subA 6 1 0 1 0 0 0 0 0 5 1 2 0 2 1 0 2 0 8 2 464
LYS 21 1IR1.subA 8 1 2 3 0 0 0 0 0 5 3 0 0 0 2 2 1 0 3 2 464
LYS 81 1IR1.subA 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 3 2 464
LYS 128 1IR1.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 1 2 464
LYS 146 1IR1.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 1 0 3 0 3 2 464
LYS 161 1IR1.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 1 0 2 0 4 2 464
LYS 164 1IR1.subA 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 3 2 464
LYS 175 1IR1.subA 8 1 0 1 1 0 0 1 1 6 2 1 0 1 3 0 2 1 1 2 464
LYS 177 1IR1.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 3 2 464
LYS 183 1IR1.subA 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 3 2 464
LYS 227 1IR1.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 2 0 0 5 2 464
LYS 236 1IR1.subA 12 1 1 2 0 1 0 1 1 9 3 1 1 2 3 0 4 0 2 2 464
LYS 252 1IR1.subA 5 1 1 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 7 2 464
LYS 305 1IR1.subA 8 1 0 1 0 2 0 2 2 5 3 1 2 3 0 0 2 0 7 2 464
LYS 316 1IR1.subA 11 1 0 1 0 2 0 2 2 8 3 0 1 1 2 1 4 0 1 2 464
LYS 334 1IR1.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 3 2 464
LYS 356 1IR1.subA 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 9 2 464
LYS 450 1IR1.subA 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 3 2 464
LYS 463 1IR1.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 2 0 0 7 2 464
LYS 466 1IR1.subA 6 0 2 2 0 0 1 0 1 3 3 0 0 0 0 1 2 0 7 2 464
LYS 8 1IR2.subA 1 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 2 468
LYS 14 1IR2.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 6 2 468
LYS 18 1IR2.subA 5 1 0 1 0 0 0 0 0 4 1 2 0 2 0 0 2 0 6 2 468
LYS 81 1IR2.subA 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 6 2 468
LYS 128 1IR2.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 0 2 468
LYS 146 1IR2.subA 5 0 0 0 0 0 0 0 0 4 0 1 0 1 2 0 1 0 2 2 468
LYS 161 1IR2.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 1 0 2 0 3 2 468
LYS 164 1IR2.subA 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 3 2 468
LYS 175 1IR2.subA 8 1 0 1 1 0 0 1 1 6 2 1 0 1 3 0 2 1 1 2 468
LYS 177 1IR2.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 2 2 468
LYS 183 1IR2.subA 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 3 2 468
LYS 227 1IR2.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 2 0 0 3 2 468
LYS 236 1IR2.subA 12 1 1 2 0 1 0 1 1 9 3 1 1 2 3 0 4 0 2 2 468
LYS 252 1IR2.subA 5 0 2 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 10 2 468
LYS 258 1IR2.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 2 1 2 0 7 2 468
LYS 316 1IR2.subA 13 1 0 1 0 2 0 2 2 10 3 0 1 1 2 2 5 0 1 2 468
LYS 334 1IR2.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 3 2 468
LYS 356 1IR2.subA 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 5 2 468
LYS 450 1IR2.subA 5 0 0 0 0 1 0 1 1 4 1 1 0 1 2 1 0 0 5 2 468
LYS 463 1IR2.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 1 2 1 0 4 2 468
LYS 466 1IR2.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 7 2 468
LYS 474 1IR2.subA 6 2 0 2 0 1 0 1 1 3 3 1 0 1 0 0 2 0 4 2 468
LYS 8 1J79.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 0 1 4 0 3 2 343
LYS 26 1J79.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 9 2 343
LYS 131 1J79.subA 5 0 3 3 0 1 0 1 1 1 4 0 0 0 0 0 1 0 3 2 343
LYS 172 1J79.subA 9 0 0 0 0 1 0 1 1 8 1 0 0 0 3 2 3 0 6 2 343
LYS 181 1J79.subA 6 2 0 2 0 0 1 0 1 3 3 1 0 1 1 0 1 0 2 2 343
LYS 226 1J79.subA 11 0 2 2 0 1 1 1 2 7 4 2 0 2 1 1 3 0 6 2 343
LYS 259 1J79.subA 12 0 1 1 0 2 2 2 4 7 5 0 0 0 6 1 0 0 3 2 343
LYS 346 1J79.subA 4 0 1 1 0 0 0 0 0 3 1 1 0 1 0 0 2 0 0 2 343
LYS 50 1JBV.subA 12 0 1 1 0 0 0 0 0 11 1 3 1 4 4 0 2 2 6 2 404
LYS 172 1JBV.subA 3 0 0 0 0 0 1 0 1 2 1 0 1 1 0 0 1 0 0 2 404
LYS 183 1JBV.subA 5 0 0 0 0 0 1 0 1 4 1 1 0 1 2 0 1 0 5 2 404
LYS 190 1JBV.subA 10 0 0 0 0 1 1 1 2 8 2 1 1 2 2 0 4 0 7 2 404
LYS 211 1JBV.subA 9 0 0 0 0 0 0 0 0 9 0 2 0 2 3 0 4 0 9 2 404
LYS 230 1JBV.subA 4 0 0 0 0 1 0 1 1 3 1 1 0 1 2 0 0 0 6 2 404
LYS 232 1JBV.subA 3 0 0 0 0 1 0 1 1 2 1 0 0 0 1 0 1 0 2 2 404
LYS 273 1JBV.subA 9 2 0 2 0 0 0 0 0 7 2 0 1 1 1 1 4 0 5 2 404
LYS 277 1JBV.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 1 0 1 0 5 2 404
LYS 303 1JBV.subA 8 1 1 2 0 0 0 0 0 6 2 1 0 1 0 1 4 0 2 2 404
LYS 329 1JBV.subA 7 0 0 0 0 0 0 0 0 7 0 1 1 2 1 2 2 0 2 2 404
LYS 388 1JBV.subA 10 1 1 2 0 1 0 1 1 7 3 1 0 1 2 1 3 0 1 2 404
LYS 18 1JBW.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 1 0 2 0 1 2 410
LYS 50 1JBW.subA 12 0 1 1 0 0 0 0 0 11 1 3 1 4 4 0 2 2 5 2 410
LYS 172 1JBW.subA 2 0 0 0 0 0 0 0 0 2 0 0 1 1 0 0 1 0 1 2 410
LYS 183 1JBW.subA 5 0 0 0 0 0 1 0 1 4 1 1 0 1 2 0 1 0 4 2 410
LYS 190 1JBW.subA 8 0 0 0 0 1 0 1 1 7 1 0 1 1 2 0 4 0 8 2 410
LYS 211 1JBW.subA 9 0 0 0 0 0 0 0 0 9 0 2 0 2 3 0 4 0 7 2 410
LYS 230 1JBW.subA 3 0 0 0 0 0 0 0 0 3 0 1 0 1 2 0 0 0 4 2 410
LYS 232 1JBW.subA 3 0 0 0 0 1 0 1 1 2 1 0 0 0 1 0 1 0 4 2 410
LYS 273 1JBW.subA 9 2 0 2 0 0 0 0 0 7 2 0 1 1 1 1 4 0 7 2 410
LYS 277 1JBW.subA 3 0 0 0 0 0 0 0 0 3 0 0 1 1 1 0 1 0 7 2 410
LYS 303 1JBW.subA 8 1 1 2 0 0 0 0 0 6 2 1 0 1 0 1 4 0 3 2 410
LYS 329 1JBW.subA 7 0 0 0 0 0 0 0 0 7 0 1 1 2 2 1 2 0 1 2 410
LYS 346 1JBW.subA 2 2 0 2 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 2 410
LYS 388 1JBW.subA 10 1 1 2 0 1 0 1 1 7 3 1 0 1 2 1 3 0 1 2 410
LYS 77 1JGM.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 8 2 333
LYS 82 1JGM.subA 9 1 1 2 0 1 0 1 1 6 3 0 0 0 3 1 2 0 14 2 333
LYS 175 1JGM.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 3 2 333
LYS 185 1JGM.subA 9 0 2 2 0 1 0 1 1 6 3 1 0 1 1 0 4 0 2 2 333
LYS 285 1JGM.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 3 2 333
LYS 294 1JGM.subA 4 0 0 0 0 1 0 1 1 3 1 0 2 2 0 0 1 0 3 2 333
LYS 339 1JGM.subA 6 0 1 1 1 0 0 1 1 4 2 0 0 0 1 1 2 0 2 2 333
LYS 3 1K1D.subA 10 0 2 2 0 0 0 0 0 8 2 1 0 1 2 1 4 0 0 2 460
LYS 6 1K1D.subA 8 1 0 1 0 0 1 0 1 6 2 0 1 1 1 0 4 0 0 2 460
LYS 24 1K1D.subA 6 1 1 2 0 0 0 0 0 4 2 1 0 1 1 0 2 0 0 2 460
LYS 27 1K1D.subA 7 2 0 2 0 0 0 0 0 5 2 1 0 1 3 0 1 0 0 2 460
LYS 38 1K1D.subA 2 0 1 1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 0 2 460
LYS 46 1K1D.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 0 2 460
LYS 71 1K1D.subA 8 2 0 2 0 0 0 0 0 6 2 2 1 3 1 1 1 0 0 2 460
LYS 97 1K1D.subA 5 0 1 1 0 0 0 0 0 4 1 1 2 3 1 0 0 0 0 2 460
LYS 102 1K1D.subA 5 0 2 2 1 0 0 1 1 2 3 0 0 0 0 0 2 0 0 2 460
LYS 103 1K1D.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 0 2 460
LYS 111 1K1D.subA 8 1 0 1 0 0 0 0 0 7 1 1 1 2 2 2 1 0 0 2 460
LYS 115 1K1D.subA 9 2 1 3 0 0 0 0 0 6 3 0 1 1 4 1 0 0 0 2 460
LYS 139 1K1D.subA 5 0 2 2 1 0 0 1 1 2 3 0 0 0 1 0 1 0 0 2 460
LYS 156 1K1D.subA 7 2 0 2 0 0 0 0 0 5 2 0 2 2 1 2 0 0 0 2 460
LYS 174 1K1D.subA 8 0 1 1 0 0 0 0 0 7 1 1 1 2 3 0 2 0 0 2 460
LYS 195 1K1D.subA 3 1 0 1 0 0 0 0 0 2 1 1 0 1 0 1 0 0 0 2 460
LYS 196 1K1D.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 1 1 0 0 2 460
LYS 248 1K1D.subA 10 0 3 3 0 0 0 0 0 7 3 1 1 2 4 0 1 0 0 2 460
LYS 255 1K1D.subA 6 0 2 2 0 0 0 0 0 4 2 0 1 1 2 0 1 0 0 2 460
LYS 277 1K1D.subA 6 0 2 2 0 0 0 0 0 4 2 1 0 1 1 2 0 0 0 2 460
LYS 284 1K1D.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 2 3 0 0 0 2 460
LYS 294 1K1D.subA 5 1 1 2 0 0 0 0 0 3 2 0 1 1 0 1 1 0 0 2 460
LYS 305 1K1D.subA 8 0 0 0 0 1 0 1 1 7 1 0 3 3 1 1 2 0 0 2 460
LYS 322 1K1D.subA 5 2 1 3 0 0 0 0 0 2 3 0 0 0 1 1 0 0 0 2 460
LYS 325 1K1D.subA 13 1 1 2 0 0 0 0 0 11 2 0 2 2 3 4 2 0 0 2 460
LYS 334 1K1D.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 1 1 2 0 0 2 460
LYS 355 1K1D.subA 7 0 0 0 1 0 0 1 1 6 1 2 1 3 1 1 1 0 0 2 460
LYS 356 1K1D.subA 6 0 1 1 1 1 0 2 2 3 3 1 0 1 2 0 0 0 0 2 460
LYS 374 1K1D.subA 6 0 0 0 0 1 0 1 1 5 1 1 0 1 1 1 2 0 0 2 460
LYS 381 1K1D.subA 6 1 0 1 1 1 0 2 2 3 3 0 0 0 1 0 2 0 0 2 460
LYS 382 1K1D.subA 12 3 0 3 1 1 0 2 2 7 5 1 0 1 3 0 3 0 0 2 460
LYS 422 1K1D.subA 7 0 1 1 0 1 0 1 1 5 2 1 0 1 1 0 3 0 0 2 460
LYS 441 1K1D.subA 6 1 1 2 0 0 0 0 0 4 2 0 1 1 1 0 2 0 0 2 460
LYS 446 1K1D.subA 4 0 0 0 0 0 0 0 0 4 0 0 0 0 2 1 1 0 0 2 460
LYS 454 1K1D.subA 7 0 3 3 0 1 0 1 1 3 4 0 0 0 1 1 1 0 0 2 460
LYS 457 1K1D.subA 8 0 1 1 0 0 0 0 0 7 1 1 0 1 3 1 2 0 0 2 460
LYS 28 1K38.subA 5 1 0 1 0 1 0 1 1 3 2 1 0 1 0 2 0 0 3 2 239
LYS 36 1K38.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 5 2 239
LYS 60 1K38.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 1 0 1 0 4 2 239
LYS 61 1K38.subA 4 0 0 0 0 2 0 2 2 2 2 0 0 0 0 1 1 0 5 2 239
LYS 125 1K38.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 1 0 1 0 0 2 239
LYS 131 1K38.subA 8 2 0 2 0 1 0 1 1 5 3 0 0 0 3 1 1 0 1 2 239
LYS 137 1K38.subA 3 0 0 0 0 2 0 2 2 1 2 0 0 0 0 0 1 0 0 2 239
LYS 138 1K38.subA 6 0 0 0 0 1 0 1 1 5 1 0 0 0 1 2 2 0 1 2 239
LYS 172 1K38.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 1 1 3 0 7 2 239
LYS 189 1K38.subA 11 1 0 1 0 3 0 3 3 7 4 0 2 2 0 2 3 0 2 2 239
LYS 205 1K38.subA 9 0 0 0 0 0 0 0 0 9 0 4 0 4 1 0 3 0 3 2 239
LYS 243 1K38.subA 11 2 0 2 0 1 0 1 1 8 3 0 2 2 1 2 3 0 4 2 239
LYS 30 1K4E.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 1 0 0 1 2 245
LYS 45 1K4E.subA 8 0 1 1 0 0 0 0 0 7 1 4 0 4 1 1 1 0 3 2 245
LYS 49 1K4E.subA 4 0 1 1 0 0 0 0 0 3 1 3 0 3 0 0 0 0 5 2 245
LYS 61 1K4E.subA 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 7 2 245
LYS 84 1K4E.subA 5 0 0 0 0 0 0 0 0 5 0 0 2 2 1 0 2 0 2 2 245
LYS 91 1K4E.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 5 2 245
LYS 95 1K4E.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 3 2 245
LYS 100 1K4E.subA 3 0 1 1 0 0 0 0 0 1 1 0 1 1 0 0 0 0 7 2 245
LYS 134 1K4E.subA 5 0 0 0 1 1 0 2 2 3 2 0 1 1 0 1 1 0 6 2 245
LYS 137 1K4E.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 9 2 245
LYS 138 1K4E.subA 4 0 0 0 1 0 0 1 1 3 1 1 0 1 0 2 0 0 5 2 245
LYS 152 1K4E.subA 7 1 1 2 0 0 0 0 0 5 2 0 1 1 2 1 1 0 8 2 245
LYS 177 1K4E.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 0 0 2 0 1 2 245
LYS 182 1K4E.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 0 1 0 3 2 245
LYS 189 1K4E.subA 10 0 1 1 0 0 1 0 1 8 2 0 2 2 0 2 4 0 1 2 245
LYS 205 1K4E.subA 10 0 0 0 0 0 0 0 0 10 0 5 0 5 0 0 4 0 3 2 245
LYS 228 1K4E.subA 9 0 4 4 0 0 0 0 0 5 4 0 0 0 1 2 2 0 0 2 245
LYS 246 1K4E.subA 6 0 1 1 0 0 0 0 0 5 1 1 1 2 0 0 3 0 3 2 245
LYS 251 1K4E.subA 6 0 0 0 0 1 0 1 1 5 1 3 0 3 1 1 0 0 9 2 245
LYS 256 1K4E.subA 6 0 2 2 0 0 0 0 0 4 2 2 0 2 0 0 2 0 2 2 245
LYS 30 1K4F.subA 5 0 1 1 0 0 0 0 0 4 1 1 1 2 1 1 0 0 1 2 244
LYS 45 1K4F.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 1 1 1 0 5 2 244
LYS 49 1K4F.subA 3 0 0 0 0 0 0 0 0 3 0 3 0 3 0 0 0 0 2 2 244
LYS 61 1K4F.subA 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 12 2 244
LYS 84 1K4F.subA 4 0 0 0 0 0 0 0 0 4 0 0 2 2 1 0 1 0 5 2 244
LYS 91 1K4F.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 7 2 244
LYS 95 1K4F.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 2 2 244
LYS 100 1K4F.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 5 2 244
LYS 134 1K4F.subA 5 0 0 0 1 1 0 2 2 3 2 0 1 1 0 1 1 0 10 2 244
LYS 137 1K4F.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 8 2 244
LYS 138 1K4F.subA 4 0 0 0 1 0 0 1 1 3 1 1 0 1 0 2 0 0 10 2 244
LYS 152 1K4F.subA 7 1 1 2 0 0 0 0 0 5 2 0 1 1 2 1 1 0 8 2 244
LYS 177 1K4F.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 0 0 2 0 6 2 244
LYS 182 1K4F.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 0 1 0 5 2 244
LYS 189 1K4F.subA 10 0 1 1 0 0 1 0 1 8 2 0 2 2 0 2 4 0 3 2 244
LYS 205 1K4F.subA 11 0 0 0 0 0 0 0 0 11 0 5 0 5 1 0 4 0 4 2 244
LYS 228 1K4F.subA 9 0 4 4 0 0 0 0 0 5 4 0 0 0 1 2 2 0 1 2 244
LYS 246 1K4F.subA 7 0 1 1 0 0 0 0 0 6 1 1 1 2 0 0 4 0 6 2 244
LYS 251 1K4F.subA 6 0 0 0 0 1 0 1 1 5 1 3 0 3 1 1 0 0 12 2 244
LYS 256 1K4F.subA 6 0 2 2 0 0 0 0 0 4 2 2 0 2 0 0 2 0 2 2 244
LYS 30 1K54.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 1 0 0 0 2 244
LYS 45 1K54.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 1 1 1 0 1 2 244
LYS 49 1K54.subA 4 0 1 1 0 0 0 0 0 3 1 3 0 3 0 0 0 0 3 2 244
LYS 61 1K54.subA 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 3 2 244
LYS 84 1K54.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 1 2 244
LYS 91 1K54.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 3 2 244
LYS 95 1K54.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 2 2 244
LYS 100 1K54.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 2 2 244
LYS 134 1K54.subA 5 0 0 0 1 1 0 2 2 3 2 0 1 1 0 1 1 0 8 2 244
LYS 137 1K54.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 4 2 244
LYS 138 1K54.subA 4 0 0 0 1 0 0 1 1 3 1 1 0 1 0 2 0 0 2 2 244
LYS 152 1K54.subA 6 1 1 2 0 0 0 0 0 4 2 0 0 0 2 1 1 0 1 2 244
LYS 177 1K54.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 0 0 2 0 1 2 244
LYS 182 1K54.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 0 1 0 2 2 244
LYS 189 1K54.subA 10 0 1 1 0 0 1 0 1 8 2 0 2 2 0 2 4 0 2 2 244
LYS 205 1K54.subA 11 0 0 0 0 0 0 0 0 11 0 5 0 5 1 0 4 0 3 2 244
LYS 228 1K54.subA 11 0 4 4 0 0 0 0 0 7 4 0 0 0 2 2 3 0 1 2 244
LYS 246 1K54.subA 7 0 1 1 0 0 0 0 0 6 1 1 1 2 0 0 4 0 1 2 244
LYS 251 1K54.subA 6 0 0 0 0 1 0 1 1 5 1 3 0 3 1 1 0 0 5 2 244
LYS 256 1K54.subA 5 0 1 1 0 0 0 0 0 4 1 2 0 2 0 0 2 0 2 2 244
LYS 30 1K55.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 1 0 0 0 2 245
LYS 45 1K55.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 1 1 1 0 0 2 245
LYS 49 1K55.subA 3 0 0 0 0 0 0 0 0 3 0 3 0 3 0 0 0 0 1 2 245
LYS 61 1K55.subA 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 4 2 245
LYS 84 1K55.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 2 2 245
LYS 91 1K55.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 5 2 245
LYS 95 1K55.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 5 2 245
LYS 100 1K55.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 5 2 245
LYS 134 1K55.subA 5 0 0 0 1 1 0 2 2 3 2 0 1 1 0 1 1 0 8 2 245
LYS 137 1K55.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 4 2 245
LYS 138 1K55.subA 4 0 0 0 1 0 0 1 1 3 1 1 0 1 0 2 0 0 3 2 245
LYS 152 1K55.subA 7 1 1 2 0 0 0 0 0 5 2 0 1 1 2 1 1 0 4 2 245
LYS 177 1K55.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 0 0 2 0 1 2 245
LYS 182 1K55.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 0 1 0 2 2 245
LYS 189 1K55.subA 10 0 1 1 0 0 1 0 1 8 2 0 2 2 0 2 4 0 2 2 245
LYS 205 1K55.subA 11 0 0 0 0 0 0 0 0 11 0 5 0 5 1 0 4 0 3 2 245
LYS 228 1K55.subA 11 0 4 4 0 0 0 0 0 7 4 0 0 0 2 2 3 0 1 2 245
LYS 246 1K55.subA 7 0 1 1 0 0 0 0 0 6 1 1 1 2 0 0 4 0 1 2 245
LYS 251 1K55.subA 6 0 0 0 0 1 0 1 1 5 1 3 0 3 1 1 0 0 6 2 245
LYS 256 1K55.subA 5 0 1 1 0 0 0 0 0 4 1 2 0 2 0 0 2 0 4 2 245
LYS 30 1K56.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 1 0 0 0 2 243
LYS 45 1K56.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 1 1 1 0 0 2 243
LYS 49 1K56.subA 3 0 0 0 0 0 0 0 0 3 0 3 0 3 0 0 0 0 0 2 243
LYS 61 1K56.subA 6 1 1 2 0 1 0 1 1 3 3 1 0 1 1 1 0 0 3 2 243
LYS 84 1K56.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 2 2 243
LYS 91 1K56.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 2 2 243
LYS 95 1K56.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 3 2 243
LYS 100 1K56.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 3 2 243
LYS 134 1K56.subA 5 0 0 0 1 1 0 2 2 3 2 0 1 1 0 1 1 0 5 2 243
LYS 137 1K56.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 4 2 243
LYS 138 1K56.subA 4 0 0 0 1 0 0 1 1 3 1 1 0 1 0 2 0 0 3 2 243
LYS 152 1K56.subA 7 1 1 2 0 0 0 0 0 5 2 0 1 1 2 1 1 0 2 2 243
LYS 177 1K56.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 0 0 2 0 1 2 243
LYS 182 1K56.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 0 1 0 3 2 243
LYS 189 1K56.subA 10 0 1 1 0 0 1 0 1 8 2 0 2 2 0 2 4 0 1 2 243
LYS 205 1K56.subA 11 0 0 0 0 0 0 0 0 11 0 5 0 5 1 0 4 0 3 2 243
LYS 228 1K56.subA 10 0 4 4 0 0 0 0 0 6 4 0 0 0 1 2 3 0 1 2 243
LYS 246 1K56.subA 7 0 1 1 0 0 0 0 0 6 1 1 1 2 0 0 4 0 0 2 243
LYS 251 1K56.subA 6 0 0 0 0 1 0 1 1 5 1 3 0 3 1 1 0 0 8 2 243
LYS 256 1K56.subA 5 0 1 1 0 0 0 0 0 4 1 2 0 2 0 0 2 0 1 2 243
LYS 30 1K57.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 1 0 0 0 2 243
LYS 45 1K57.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 1 1 1 0 0 2 243
LYS 49 1K57.subA 3 0 0 0 0 0 0 0 0 3 0 3 0 3 0 0 0 0 0 2 243
LYS 61 1K57.subA 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 2 2 243
LYS 84 1K57.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 1 2 243
LYS 91 1K57.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 1 2 243
LYS 95 1K57.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 3 2 243
LYS 100 1K57.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 4 2 243
LYS 134 1K57.subA 6 0 0 0 1 1 0 2 2 4 2 0 1 1 0 1 2 0 3 2 243
LYS 137 1K57.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 5 2 243
LYS 138 1K57.subA 4 0 0 0 1 0 0 1 1 3 1 1 0 1 0 2 0 0 3 2 243
LYS 152 1K57.subA 7 1 1 2 0 0 0 0 0 5 2 0 1 1 2 1 1 0 1 2 243
LYS 177 1K57.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 0 0 2 0 0 2 243
LYS 182 1K57.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 0 1 0 3 2 243
LYS 189 1K57.subA 10 0 1 1 0 0 1 0 1 8 2 0 2 2 0 2 4 0 2 2 243
LYS 205 1K57.subA 11 0 0 0 0 0 0 0 0 11 0 5 0 5 1 0 4 0 2 2 243
LYS 228 1K57.subA 10 0 4 4 0 0 0 0 0 6 4 0 0 0 1 2 3 0 0 2 243
LYS 246 1K57.subA 9 0 1 1 0 0 0 0 0 8 1 1 1 2 1 0 5 0 2 2 243
LYS 251 1K57.subA 6 0 0 0 0 1 0 1 1 5 1 3 0 3 1 1 0 0 5 2 243
LYS 256 1K57.subA 5 0 1 1 0 0 0 0 0 4 1 2 0 2 0 0 2 0 2 2 243
LYS 30 1K6S.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 1 0 0 1 2 244
LYS 45 1K6S.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 1 1 1 0 3 2 244
LYS 49 1K6S.subA 3 0 0 0 0 0 0 0 0 3 0 3 0 3 0 0 0 0 0 2 244
LYS 61 1K6S.subA 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 4 2 244
LYS 84 1K6S.subA 5 0 0 0 0 0 0 0 0 5 0 0 2 2 1 0 2 0 1 2 244
LYS 91 1K6S.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 6 2 244
LYS 95 1K6S.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 5 2 244
LYS 100 1K6S.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 6 2 244
LYS 134 1K6S.subA 5 0 0 0 1 1 0 2 2 3 2 0 1 1 0 1 1 0 3 2 244
LYS 137 1K6S.subA 4 0 0 0 0 0 0 0 0 4 0 0 2 2 0 1 1 0 5 2 244
LYS 138 1K6S.subA 4 0 0 0 1 0 0 1 1 3 1 1 0 1 0 2 0 0 3 2 244
LYS 152 1K6S.subA 7 1 1 2 0 0 0 0 0 5 2 0 1 1 2 1 1 0 4 2 244
LYS 177 1K6S.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 0 0 2 0 0 2 244
LYS 182 1K6S.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 0 1 0 3 2 244
LYS 189 1K6S.subA 10 0 1 1 0 0 1 0 1 8 2 0 2 2 0 2 4 0 0 2 244
LYS 205 1K6S.subA 11 0 0 0 0 0 0 0 0 11 0 5 0 5 1 0 4 0 4 2 244
LYS 228 1K6S.subA 10 0 4 4 0 0 0 0 0 6 4 0 0 0 1 2 3 0 0 2 244
LYS 246 1K6S.subA 7 0 1 1 0 0 0 0 0 6 1 1 1 2 0 0 4 0 4 2 244
LYS 251 1K6S.subA 6 0 0 0 0 1 0 1 1 5 1 3 0 3 1 1 0 0 5 2 244
LYS 256 1K6S.subA 6 0 2 2 0 0 0 0 0 4 2 2 0 2 0 0 2 0 1 2 244
LYS 20 1KRB.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 1 2 566
LYS 44 1KRB.subC 8 0 1 1 1 0 0 1 1 6 2 0 1 1 0 2 3 0 1 2 566
LYS 49 1KRB.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 1 2 566
LYS 83 1KRB.subC 14 2 0 2 0 0 0 0 0 12 2 1 1 2 4 2 4 0 3 2 566
LYS 89 1KRB.subC 10 3 0 3 0 1 0 1 1 6 4 0 0 0 3 1 2 0 1 2 566
LYS 98 1KRB.subC 14 1 2 3 0 0 0 0 0 11 3 3 2 5 3 0 3 0 1 2 566
LYS 124 1KRB.subC 9 0 1 1 0 0 0 0 0 8 1 1 0 1 4 1 2 0 1 2 566
LYS 196 1KRB.subC 7 0 0 0 0 0 0 0 0 7 0 0 2 2 4 1 0 0 0 2 566
LYS 382 1KRB.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 0 2 566
LYS 401 1KRB.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 1 2 566
LYS 406 1KRB.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 2 3 2 0 0 2 566
LYS 428 1KRB.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 1 2 566
LYS 443 1KRB.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 566
LYS 449 1KRB.subC 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 0 2 566
LYS 515 1KRB.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 2 2 566
LYS 522 1KRB.subC 17 0 0 0 0 0 1 0 1 16 1 2 3 5 7 1 3 0 3 2 566
LYS 39 1L6F.subA 9 0 1 1 0 0 0 0 0 8 1 0 0 0 3 1 4 0 2 2 382
LYS 76 1L6F.subA 7 1 1 2 0 1 0 1 1 4 3 0 0 0 2 0 2 0 0 2 382
LYS 140 1L6F.subA 4 2 0 2 0 0 0 0 0 2 2 0 0 0 1 0 1 0 3 2 382
LYS 146 1L6F.subA 7 0 3 3 0 1 0 1 1 3 4 1 0 1 0 1 1 0 0 2 382
LYS 234 1L6F.subA 8 1 0 1 0 0 0 0 0 7 1 1 0 1 4 0 2 0 1 2 382
LYS 242 1L6F.subA 7 1 1 2 0 0 0 0 0 5 2 0 2 2 2 0 1 0 1 2 382
LYS 255 1L6F.subA 5 0 0 0 0 0 1 0 1 4 1 1 0 1 0 0 3 0 1 2 382
LYS 256 1L6F.subA 5 0 2 2 0 0 0 0 0 3 2 0 1 1 0 1 1 0 0 2 382
LYS 262 1L6F.subA 6 0 1 1 0 0 0 0 0 5 1 1 0 1 2 1 1 0 0 2 382
LYS 303 1L6F.subA 6 0 1 1 0 0 1 0 1 4 2 0 1 1 1 0 2 0 1 2 382
LYS 328 1L6F.subA 9 0 0 0 0 1 1 1 2 7 2 3 1 4 1 1 1 0 3 2 382
LYS 372 1L6F.subA 3 1 0 1 0 1 1 1 2 0 3 0 0 0 0 0 0 0 2 2 382
LYS 39 1L6G.subA 9 0 1 1 0 0 0 0 0 8 1 0 0 0 3 1 4 0 2 2 382
LYS 76 1L6G.subA 7 1 1 2 0 1 0 1 1 4 3 0 0 0 2 0 2 0 0 2 382
LYS 140 1L6G.subA 4 2 0 2 0 0 0 0 0 2 2 0 0 0 1 0 1 0 2 2 382
LYS 146 1L6G.subA 7 0 2 2 0 1 0 1 1 4 3 1 0 1 1 1 1 0 0 2 382
LYS 234 1L6G.subA 8 1 0 1 0 0 0 0 0 7 1 1 0 1 4 0 2 0 2 2 382
LYS 242 1L6G.subA 7 1 1 2 0 0 0 0 0 5 2 0 2 2 2 0 1 0 1 2 382
LYS 255 1L6G.subA 5 0 0 0 0 0 1 0 1 4 1 1 0 1 0 0 3 0 2 2 382
LYS 256 1L6G.subA 5 0 2 2 0 0 0 0 0 3 2 0 1 1 0 1 1 0 1 2 382
LYS 262 1L6G.subA 6 0 1 1 0 0 0 0 0 5 1 1 0 1 2 1 1 0 1 2 382
LYS 303 1L6G.subA 6 0 1 1 0 0 1 0 1 4 2 0 1 1 1 0 2 0 1 2 382
LYS 328 1L6G.subA 9 0 0 0 0 1 1 1 2 7 2 3 1 4 1 1 1 0 2 2 382
LYS 372 1L6G.subA 3 1 0 1 0 1 1 1 2 0 3 0 0 0 0 0 0 0 0 2 382
LYS 56 1M6K.subA 8 0 2 2 0 0 0 0 0 6 2 0 2 2 2 1 1 0 4 2 250
LYS 58 1M6K.subA 11 0 1 1 0 0 0 0 0 10 1 2 1 3 5 0 2 0 4 2 250
LYS 87 1M6K.subA 3 1 0 1 0 0 0 0 0 2 1 1 1 2 0 0 0 0 4 2 250
LYS 91 1M6K.subA 4 0 0 0 0 0 0 0 0 4 0 0 1 1 0 2 1 0 8 2 250
LYS 94 1M6K.subA 3 1 0 1 0 0 0 0 0 2 1 1 0 1 0 1 0 0 5 2 250
LYS 97 1M6K.subA 9 1 0 1 0 0 0 0 0 8 1 1 1 2 2 3 1 0 7 2 250
LYS 109 1M6K.subA 6 0 1 1 0 0 0 0 0 5 1 2 2 4 1 0 0 0 8 2 250
LYS 126 1M6K.subA 5 0 1 1 0 0 0 0 0 4 1 0 1 1 0 0 3 0 3 2 250
LYS 131 1M6K.subA 6 0 0 0 0 0 0 0 0 6 0 0 1 1 2 1 2 0 7 2 250
LYS 133 1M6K.subA 5 0 0 0 1 0 0 1 1 4 1 0 2 2 0 0 2 0 7 2 250
LYS 137 1M6K.subA 5 1 0 1 1 0 0 1 1 3 2 0 2 2 0 0 1 0 7 2 250
LYS 166 1M6K.subA 8 0 0 0 0 0 0 0 0 8 0 2 2 4 1 0 3 0 9 2 250
LYS 178 1M6K.subA 8 1 0 1 0 1 0 1 1 6 2 0 2 2 0 2 2 0 8 2 250
LYS 187 1M6K.subA 7 1 1 2 0 0 0 0 0 5 2 1 1 2 1 1 1 0 7 2 250
LYS 208 1M6K.subA 8 1 0 1 0 0 0 0 0 7 1 2 1 3 0 2 2 0 4 2 250
LYS 212 1M6K.subA 11 0 1 1 0 0 0 0 0 10 1 4 0 4 2 0 3 0 5 2 250
LYS 236 1M6K.subA 4 0 0 0 0 0 0 0 0 4 0 3 1 4 0 0 0 0 5 2 250
LYS 240 1M6K.subA 7 0 0 0 0 0 1 0 1 6 1 0 0 0 2 1 3 0 4 2 250
LYS 260 1M6K.subA 11 0 1 1 0 0 0 0 0 10 1 3 1 4 1 1 4 0 7 2 250
LYS 262 1M6K.subA 7 0 0 0 0 0 0 0 0 7 0 2 0 2 1 1 3 0 9 2 250
LYS 263 1M6K.subA 3 0 0 0 0 0 0 0 0 3 0 1 1 2 0 0 1 0 12 2 250
LYS 6 1NFG.subA 7 2 0 2 0 0 0 0 0 5 2 1 1 2 0 0 3 0 0 2 457
LYS 24 1NFG.subA 5 2 0 2 0 0 0 0 0 3 2 1 0 1 0 0 2 0 2 2 457
LYS 27 1NFG.subA 7 2 0 2 0 0 0 0 0 5 2 2 0 2 2 0 1 0 3 2 457
LYS 106 1NFG.subA 6 1 1 2 0 0 0 0 0 4 2 0 0 0 2 1 1 0 0 2 457
LYS 114 1NFG.subA 8 0 0 0 0 0 0 0 0 8 0 1 0 1 6 0 1 0 3 2 457
LYS 166 1NFG.subA 6 0 0 0 0 0 0 0 0 6 0 2 0 2 1 0 3 0 3 2 457
LYS 170 1NFG.subA 8 1 2 3 0 0 0 0 0 5 3 1 0 1 1 1 2 0 1 2 457
LYS 173 1NFG.subA 4 1 1 2 0 0 0 0 0 2 2 1 0 1 0 0 1 0 0 2 457
LYS 194 1NFG.subA 5 1 1 2 0 0 0 0 0 3 2 0 0 0 0 2 1 0 3 2 457
LYS 200 1NFG.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 1 2 0 0 0 2 457
LYS 251 1NFG.subA 9 0 1 1 0 0 0 0 0 8 1 1 0 1 2 0 5 0 2 2 457
LYS 270 1NFG.subA 7 0 1 1 0 1 0 1 1 5 2 1 0 1 2 0 2 0 1 2 457
LYS 282 1NFG.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 2 3 0 0 2 2 457
LYS 292 1NFG.subA 3 0 0 0 0 0 1 0 1 2 1 1 0 1 1 0 0 0 1 2 457
LYS 293 1NFG.subA 5 2 2 4 0 0 0 0 0 1 4 0 0 0 1 0 0 0 2 2 457
LYS 320 1NFG.subA 4 2 0 2 0 0 0 0 0 2 2 0 0 0 1 1 0 0 0 2 457
LYS 323 1NFG.subA 13 1 0 1 0 0 1 0 1 11 2 0 1 1 3 4 3 0 0 2 457
LYS 371 1NFG.subA 5 0 0 0 0 1 0 1 1 4 1 1 0 1 1 1 1 0 5 2 457
LYS 379 1NFG.subA 12 3 0 3 0 1 0 1 1 8 4 1 1 2 2 0 4 0 2 2 457
LYS 419 1NFG.subA 5 0 0 0 0 0 1 0 1 4 1 0 0 0 1 0 3 0 0 2 457
LYS 421 1NFG.subA 8 0 2 2 0 0 0 0 0 6 2 0 2 2 3 0 1 0 0 2 457
LYS 425 1NFG.subA 9 1 1 2 0 1 0 1 1 6 3 1 0 1 3 0 2 0 0 2 457
LYS 432 1NFG.subA 9 1 1 2 0 1 0 1 1 6 3 0 0 0 2 0 4 0 3 2 457
LYS 448 1NFG.subA 6 2 0 2 0 1 0 1 1 3 3 1 0 1 1 1 0 0 0 2 457
LYS 451 1NFG.subA 6 1 0 1 0 1 0 1 1 4 2 0 0 0 1 1 2 0 0 2 457
LYS 454 1NFG.subA 7 1 0 1 0 1 0 1 1 5 2 2 0 2 2 1 0 0 2 2 457
LYS 456 1NFG.subA 11 0 0 0 0 2 1 2 3 8 3 2 1 3 1 2 2 0 0 2 457
LYS 39 1NIU.subA 9 0 1 1 0 0 0 0 0 8 1 0 0 0 3 1 4 0 2 2 382
LYS 76 1NIU.subA 7 1 1 2 0 1 0 1 1 4 3 0 0 0 2 0 2 0 2 2 382
LYS 140 1NIU.subA 5 2 0 2 0 1 0 1 1 2 3 0 0 0 1 0 1 0 4 2 382
LYS 146 1NIU.subA 5 0 2 2 0 1 0 1 1 2 3 1 0 1 0 1 0 0 0 2 382
LYS 234 1NIU.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 3 0 1 0 1 2 382
LYS 242 1NIU.subA 6 1 1 2 0 0 0 0 0 4 2 0 1 1 2 0 1 0 1 2 382
LYS 255 1NIU.subA 5 0 0 0 0 0 1 0 1 4 1 1 0 1 0 0 3 0 2 2 382
LYS 256 1NIU.subA 5 0 2 2 0 0 0 0 0 3 2 0 1 1 0 1 1 0 1 2 382
LYS 262 1NIU.subA 6 0 1 1 0 0 0 0 0 5 1 1 0 1 2 1 1 0 1 2 382
LYS 303 1NIU.subA 6 0 1 1 0 0 1 0 1 4 2 0 1 1 1 0 2 0 2 2 382
LYS 328 1NIU.subA 9 0 1 1 0 1 1 1 2 6 3 3 1 4 0 1 1 0 1 2 382
LYS 372 1NIU.subA 3 1 0 1 0 1 1 1 2 0 3 0 0 0 0 0 0 0 1 2 382
LYS 6 1OIR.subA 7 0 0 0 0 0 0 0 0 7 0 0 1 1 0 3 3 0 0 2 287
LYS 9 1OIR.subA 7 0 2 2 0 0 0 0 0 5 2 0 0 0 1 1 3 0 1 2 287
LYS 20 1OIR.subA 8 0 1 1 0 0 0 0 0 7 1 0 0 0 1 2 4 0 1 2 287
LYS 24 1OIR.subA 8 0 1 1 0 1 0 1 1 6 2 0 3 3 0 1 2 0 0 2 287
LYS 34 1OIR.subA 6 0 0 0 1 0 0 1 1 5 1 0 0 0 0 2 2 0 1 2 287
LYS 56 1OIR.subA 9 1 1 2 0 0 0 0 0 7 2 1 0 1 0 0 6 0 1 2 287
LYS 65 1OIR.subA 6 0 1 1 0 0 0 0 0 5 1 0 1 1 0 0 4 0 5 2 287
LYS 75 1OIR.subA 9 0 1 1 1 1 0 2 2 6 3 1 1 2 0 1 3 0 0 2 287
LYS 88 1OIR.subA 5 2 0 2 0 0 0 0 0 3 2 0 1 1 1 0 1 0 0 2 287
LYS 89 1OIR.subA 3 1 0 1 0 0 0 0 0 2 1 0 1 1 0 1 0 0 0 2 287
LYS 105 1OIR.subA 8 1 0 1 0 0 0 0 0 7 1 1 0 1 1 2 3 0 2 2 287
LYS 129 1OIR.subA 9 1 0 1 0 0 0 0 0 8 1 2 2 4 1 1 2 0 5 2 287
LYS 142 1OIR.subA 10 0 1 1 0 0 0 0 0 9 1 0 1 1 2 0 6 0 3 2 287
LYS 178 1OIR.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 0 2 287
LYS 237 1OIR.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 2 2 0 0 0 2 287
LYS 242 1OIR.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 2 2 0 0 2 2 287
LYS 250 1OIR.subA 5 1 0 1 0 0 0 0 0 4 1 1 1 2 0 1 1 0 1 2 287
LYS 273 1OIR.subA 4 1 0 1 0 1 1 1 2 1 3 0 1 1 0 0 0 0 2 2 287
LYS 278 1OIR.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 3 0 0 0 3 2 287
LYS 291 1OIR.subA 7 0 1 1 0 0 0 0 0 6 1 2 1 3 3 0 0 0 4 2 287
LYS 34 1ONW.subA 8 1 0 1 0 0 0 0 0 7 1 0 2 2 3 0 2 0 3 2 376
LYS 119 1ONW.subA 11 0 1 1 0 0 0 0 0 10 1 2 0 2 3 0 5 0 5 2 376
LYS 150 1ONW.subA 5 1 1 2 0 0 0 0 0 3 2 2 0 2 0 0 1 0 2 2 376
LYS 194 1ONW.subA 15 2 0 2 0 1 0 1 1 12 3 0 0 0 7 0 5 0 4 2 376
LYS 206 1ONW.subA 1 0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 0 0 4 2 376
LYS 207 1ONW.subA 5 0 0 0 0 0 0 0 0 5 0 2 1 3 2 0 0 0 1 2 376
LYS 225 1ONW.subA 9 0 0 0 0 1 0 1 1 8 1 2 0 2 3 0 3 0 10 2 376
LYS 247 1ONW.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 1 2 1 0 3 2 376
LYS 319 1ONW.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 0 0 2 0 6 2 376
LYS 345 1ONW.subA 12 2 0 2 0 1 0 1 1 9 3 3 0 3 3 1 2 0 5 2 376
LYS 373 1ONW.subA 8 0 1 1 0 1 0 1 1 6 2 0 0 0 3 0 3 0 5 2 376
LYS 377 1ONW.subA 6 1 1 2 0 0 0 0 0 4 2 0 1 1 1 0 2 0 3 2 376
LYS 380 1ONW.subA 4 1 1 2 0 0 0 0 0 2 2 0 0 0 2 0 0 0 1 2 376
LYS 384 1ONW.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 2 0 1 0 3 2 376
LYS 34 1ONX.subA 8 1 0 1 0 0 0 0 0 7 1 0 2 2 3 0 2 0 1 2 389
LYS 119 1ONX.subA 9 0 1 1 0 0 0 0 0 8 1 2 0 2 1 0 5 0 3 2 389
LYS 150 1ONX.subA 5 1 1 2 0 0 0 0 0 3 2 2 0 2 0 0 1 0 1 2 389
LYS 194 1ONX.subA 15 2 0 2 0 1 0 1 1 12 3 0 0 0 7 0 5 0 1 2 389
LYS 206 1ONX.subA 1 0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 0 0 1 2 389
LYS 207 1ONX.subA 5 0 0 0 0 0 0 0 0 5 0 2 1 3 2 0 0 0 1 2 389
LYS 225 1ONX.subA 9 0 0 0 0 1 0 1 1 8 1 1 1 2 3 0 3 0 2 2 389
LYS 247 1ONX.subA 7 0 1 1 0 1 0 1 1 5 2 0 0 0 2 2 1 0 1 2 389
LYS 319 1ONX.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 0 0 2 0 0 2 389
LYS 345 1ONX.subA 12 2 0 2 0 1 0 1 1 9 3 3 0 3 3 1 2 0 4 2 389
LYS 373 1ONX.subA 7 0 1 1 0 1 0 1 1 5 2 0 0 0 2 0 3 0 1 2 389
LYS 377 1ONX.subA 6 1 1 2 0 0 0 0 0 4 2 0 1 1 1 0 2 0 1 2 389
LYS 380 1ONX.subA 4 1 1 2 0 0 0 0 0 2 2 0 0 0 2 0 0 0 0 2 389
LYS 384 1ONX.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 0 2 389
LYS 77 1P6B.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 0 2 330
LYS 82 1P6B.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 2 1 2 0 7 2 330
LYS 175 1P6B.subA 3 0 0 0 0 0 0 0 0 3 0 0 1 1 2 0 0 0 1 2 330
LYS 185 1P6B.subA 8 0 2 2 0 1 0 1 1 5 3 0 0 0 1 0 4 0 2 2 330
LYS 285 1P6B.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 2 2 330
LYS 294 1P6B.subA 4 0 0 0 0 1 0 1 1 3 1 0 2 2 0 0 1 0 2 2 330
LYS 339 1P6B.subA 7 0 1 1 1 1 0 2 2 4 3 0 0 0 1 1 2 0 0 2 330
LYS 77 1P6C.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 2 2 331
LYS 82 1P6C.subA 10 1 1 2 0 1 0 1 1 7 3 0 0 0 4 1 2 0 6 2 331
LYS 175 1P6C.subA 3 0 0 0 0 0 0 0 0 3 0 0 1 1 2 0 0 0 0 2 331
LYS 185 1P6C.subA 8 0 2 2 0 1 0 1 1 5 3 0 0 0 1 0 4 0 2 2 331
LYS 285 1P6C.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 0 2 331
LYS 294 1P6C.subA 4 0 0 0 0 1 0 1 1 3 1 0 2 2 0 0 1 0 2 2 331
LYS 339 1P6C.subA 7 0 1 1 1 1 0 2 2 4 3 0 0 0 1 1 2 0 0 2 331
LYS 34 1PO9.subA 8 1 0 1 0 0 0 0 0 7 1 0 2 2 3 0 2 0 2 2 373
LYS 119 1PO9.subA 11 0 1 1 0 0 0 0 0 10 1 2 0 2 3 0 5 0 4 2 373
LYS 150 1PO9.subA 5 1 1 2 0 0 0 0 0 3 2 2 0 2 0 0 1 0 1 2 373
LYS 194 1PO9.subA 15 2 0 2 0 1 0 1 1 12 3 0 0 0 7 0 5 0 7 2 373
LYS 206 1PO9.subA 1 0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 0 0 4 2 373
LYS 207 1PO9.subA 5 0 0 0 0 0 0 0 0 5 0 2 1 3 2 0 0 0 3 2 373
LYS 225 1PO9.subA 9 0 0 0 0 1 0 1 1 8 1 2 0 2 3 0 3 0 7 2 373
LYS 247 1PO9.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 1 2 1 0 2 2 373
LYS 319 1PO9.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 0 0 2 0 2 2 373
LYS 345 1PO9.subA 12 2 0 2 0 1 0 1 1 9 3 3 0 3 3 1 2 0 4 2 373
LYS 373 1PO9.subA 5 0 0 0 0 1 0 1 1 4 1 0 0 0 2 0 2 0 1 2 373
LYS 377 1PO9.subA 5 1 0 1 0 0 0 0 0 4 1 0 1 1 1 0 2 0 0 2 373
LYS 380 1PO9.subA 4 1 1 2 0 0 0 0 0 2 2 0 0 0 2 0 0 0 0 2 373
LYS 384 1PO9.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 1 2 373
LYS 34 1POJ.subA 8 1 0 1 0 0 0 0 0 7 1 0 2 2 3 0 2 0 0 2 388
LYS 119 1POJ.subA 10 0 1 1 0 0 0 0 0 9 1 2 0 2 2 0 5 0 1 2 388
LYS 150 1POJ.subA 5 1 1 2 0 0 0 0 0 3 2 2 0 2 0 0 1 0 1 2 388
LYS 194 1POJ.subA 15 2 0 2 0 1 0 1 1 12 3 0 0 0 7 0 5 0 0 2 388
LYS 206 1POJ.subA 1 0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 0 0 0 2 388
LYS 207 1POJ.subA 5 0 0 0 0 0 0 0 0 5 0 2 1 3 2 0 0 0 0 2 388
LYS 225 1POJ.subA 9 0 0 0 0 1 0 1 1 8 1 1 1 2 3 0 3 0 0 2 388
LYS 247 1POJ.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 1 2 1 0 0 2 388
LYS 319 1POJ.subA 3 1 0 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 0 2 388
LYS 345 1POJ.subA 11 2 0 2 0 0 0 0 0 9 2 3 0 3 3 1 2 0 0 2 388
LYS 373 1POJ.subA 4 0 0 0 0 0 0 0 0 4 0 0 0 0 2 0 2 0 0 2 388
LYS 377 1POJ.subA 6 1 1 2 0 0 0 0 0 4 2 0 1 1 1 0 2 0 0 2 388
LYS 380 1POJ.subA 4 1 1 2 0 0 0 0 0 2 2 0 0 0 2 0 0 0 0 2 388
LYS 384 1POJ.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 1 2 388
LYS 34 1POK.subA 8 1 0 1 0 0 0 0 0 7 1 0 2 2 3 0 2 0 0 2 377
LYS 119 1POK.subA 11 0 1 1 0 0 0 0 0 10 1 2 0 2 3 0 5 0 1 2 377
LYS 150 1POK.subA 5 1 1 2 0 0 0 0 0 3 2 2 0 2 0 0 1 0 1 2 377
LYS 194 1POK.subA 15 2 0 2 0 1 0 1 1 12 3 0 0 0 7 0 5 0 1 2 377
LYS 206 1POK.subA 1 0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 0 0 0 2 377
LYS 207 1POK.subA 5 0 0 0 0 0 0 0 0 5 0 2 1 3 2 0 0 0 1 2 377
LYS 225 1POK.subA 9 0 0 0 0 1 0 1 1 8 1 2 0 2 3 0 3 0 2 2 377
LYS 247 1POK.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 1 2 1 0 0 2 377
LYS 319 1POK.subA 5 1 0 1 0 1 0 1 1 3 2 0 1 1 0 0 2 0 1 2 377
LYS 345 1POK.subA 12 2 0 2 0 1 0 1 1 9 3 3 0 3 3 1 2 0 1 2 377
LYS 373 1POK.subA 5 0 0 0 0 1 0 1 1 4 1 0 0 0 2 0 2 0 0 2 377
LYS 377 1POK.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 0 2 377
LYS 380 1POK.subA 4 1 1 2 0 0 0 0 0 2 2 0 0 0 2 0 0 0 0 2 377
LYS 384 1POK.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 0 2 377
LYS 9 1PU6.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 1 1 2 0 1 2 217
LYS 12 1PU6.subA 6 0 0 0 0 0 0 0 0 6 0 1 0 1 0 1 3 0 2 2 217
LYS 18 1PU6.subA 3 1 0 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 1 2 217
LYS 30 1PU6.subA 7 1 1 2 0 0 0 0 0 5 2 0 1 1 2 1 1 0 3 2 217
LYS 44 1PU6.subA 4 0 1 1 0 0 0 0 0 3 1 1 0 1 1 1 0 0 1 2 217
LYS 50 1PU6.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 0 2 217
LYS 56 1PU6.subA 7 0 2 2 0 0 0 0 0 5 2 0 2 2 0 0 3 0 1 2 217
LYS 70 1PU6.subA 6 0 1 1 0 0 0 0 0 5 1 0 1 1 0 1 3 0 4 2 217
LYS 71 1PU6.subA 5 0 0 0 0 0 0 0 0 5 0 0 1 1 0 1 3 0 4 2 217
LYS 79 1PU6.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 0 0 3 0 6 2 217
LYS 93 1PU6.subA 14 0 0 0 0 1 0 1 1 13 1 3 3 6 2 2 3 0 1 2 217
LYS 95 1PU6.subA 6 1 0 1 0 1 0 1 1 4 2 0 2 2 1 1 0 0 3 2 217
LYS 106 1PU6.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 1 0 1 0 0 2 217
LYS 115 1PU6.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 2 2 1 0 5 2 217
LYS 127 1PU6.subA 2 0 0 0 0 0 0 0 0 2 0 0 1 1 1 0 0 0 5 2 217
LYS 131 1PU6.subA 7 2 1 3 0 1 0 1 1 3 4 0 0 0 1 1 1 0 5 2 217
LYS 144 1PU6.subA 9 0 2 2 0 0 0 0 0 7 2 0 1 1 3 0 3 0 2 2 217
LYS 151 1PU6.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 0 1 1 0 2 2 217
LYS 158 1PU6.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 0 1 3 0 1 2 217
LYS 159 1PU6.subA 4 0 0 0 1 0 0 1 1 3 1 0 0 0 0 1 2 0 1 2 217
LYS 176 1PU6.subA 6 0 2 2 0 0 1 0 1 3 3 0 1 1 1 1 0 0 2 2 217
LYS 211 1PU6.subA 7 0 2 2 0 0 0 0 0 5 2 1 1 2 0 2 1 0 3 2 217
LYS 213 1PU6.subA 6 0 0 0 1 0 0 1 1 5 1 1 1 2 0 1 2 0 2 2 217
LYS 217 1PU6.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 3 0 0 2 217
LYS 9 1PU7.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 1 1 2 0 3 2 215
LYS 12 1PU7.subA 6 0 0 0 0 0 0 0 0 6 0 1 0 1 0 1 3 0 2 2 215
LYS 18 1PU7.subA 3 1 0 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 2 2 215
LYS 30 1PU7.subA 8 1 1 2 0 0 0 0 0 6 2 0 2 2 2 1 1 0 6 2 215
LYS 44 1PU7.subA 4 0 1 1 0 0 0 0 0 3 1 1 0 1 1 1 0 0 1 2 215
LYS 50 1PU7.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 1 2 215
LYS 56 1PU7.subA 7 0 2 2 0 0 0 0 0 5 2 0 2 2 0 0 3 0 3 2 215
LYS 70 1PU7.subA 6 0 1 1 0 0 0 0 0 5 1 0 1 1 0 1 3 0 9 2 215
LYS 71 1PU7.subA 6 0 1 1 0 0 0 0 0 5 1 0 1 1 0 1 3 0 4 2 215
LYS 79 1PU7.subA 6 0 1 1 0 0 0 0 0 5 1 1 0 1 1 0 3 0 2 2 215
LYS 93 1PU7.subA 14 0 0 0 0 1 0 1 1 13 1 3 3 6 2 2 3 0 3 2 215
LYS 95 1PU7.subA 6 1 0 1 0 1 0 1 1 4 2 0 2 2 1 1 0 0 8 2 215
LYS 106 1PU7.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 1 0 1 0 1 2 215
LYS 115 1PU7.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 2 2 1 0 6 2 215
LYS 127 1PU7.subA 2 0 0 0 0 0 0 0 0 2 0 0 1 1 1 0 0 0 6 2 215
LYS 131 1PU7.subA 7 2 1 3 0 1 0 1 1 3 4 0 0 0 1 1 1 0 6 2 215
LYS 144 1PU7.subA 9 0 2 2 0 0 0 0 0 7 2 0 1 1 3 0 3 0 4 2 215
LYS 151 1PU7.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 3 2 215
LYS 158 1PU7.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 0 1 3 0 2 2 215
LYS 159 1PU7.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 0 1 3 0 4 2 215
LYS 176 1PU7.subA 6 0 2 2 0 0 1 0 1 3 3 0 1 1 1 1 0 0 2 2 215
LYS 211 1PU7.subA 6 0 2 2 0 0 0 0 0 4 2 1 1 2 0 1 1 0 1 2 215
LYS 213 1PU7.subA 3 0 0 0 0 0 0 0 0 3 0 1 1 2 0 0 1 0 0 2 215
LYS 9 1PU8.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 1 1 2 0 1 2 216
LYS 12 1PU8.subA 6 0 0 0 0 0 0 0 0 6 0 1 0 1 0 1 3 0 1 2 216
LYS 18 1PU8.subA 3 1 0 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 0 2 216
LYS 30 1PU8.subA 7 1 1 2 0 0 0 0 0 5 2 0 1 1 2 1 1 0 2 2 216
LYS 44 1PU8.subA 4 0 1 1 0 0 0 0 0 3 1 1 0 1 1 1 0 0 0 2 216
LYS 50 1PU8.subA 6 0 2 2 0 0 0 0 0 4 2 1 0 1 2 0 1 0 0 2 216
LYS 56 1PU8.subA 7 0 2 2 0 0 0 0 0 5 2 0 2 2 0 0 3 0 2 2 216
LYS 70 1PU8.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 0 1 3 0 1 2 216
LYS 71 1PU8.subA 5 0 0 0 0 0 0 0 0 5 0 0 1 1 0 1 3 0 0 2 216
LYS 79 1PU8.subA 6 0 1 1 0 0 0 0 0 5 1 1 0 1 1 0 3 0 1 2 216
LYS 93 1PU8.subA 14 0 0 0 0 1 0 1 1 13 1 3 3 6 2 2 3 0 2 2 216
LYS 95 1PU8.subA 6 1 0 1 0 1 0 1 1 4 2 0 2 2 1 1 0 0 2 2 216
LYS 106 1PU8.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 1 0 1 0 1 2 216
LYS 115 1PU8.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 2 2 1 0 3 2 216
LYS 127 1PU8.subA 2 0 0 0 0 0 0 0 0 2 0 0 1 1 1 0 0 0 3 2 216
LYS 131 1PU8.subA 7 2 1 3 0 1 0 1 1 3 4 0 0 0 1 1 1 0 2 2 216
LYS 144 1PU8.subA 9 0 2 2 0 0 0 0 0 7 2 0 1 1 3 0 3 0 1 2 216
LYS 151 1PU8.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 0 1 1 0 3 2 216
LYS 158 1PU8.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 0 1 3 0 1 2 216
LYS 159 1PU8.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 0 1 2 0 0 2 216
LYS 176 1PU8.subA 7 0 2 2 0 0 1 0 1 4 3 0 1 1 1 2 0 0 1 2 216
LYS 211 1PU8.subA 7 0 2 2 0 0 0 0 0 5 2 1 1 2 0 2 1 0 0 2 216
LYS 213 1PU8.subA 3 0 0 0 0 0 0 0 0 3 0 1 1 2 0 0 1 0 0 2 216
LYS 77 1QW7.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 6 2 336
LYS 82 1QW7.subA 9 1 1 2 0 1 0 1 1 6 3 0 0 0 3 1 2 0 10 2 336
LYS 175 1QW7.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 5 2 336
LYS 185 1QW7.subA 8 0 2 2 0 1 0 1 1 5 3 0 0 0 1 0 4 0 2 2 336
LYS 285 1QW7.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 3 2 336
LYS 294 1QW7.subA 5 0 0 0 0 1 0 1 1 4 1 0 2 2 0 1 1 0 5 2 336
LYS 339 1QW7.subA 7 0 1 1 1 0 0 1 1 5 2 0 0 0 1 1 3 0 1 2 336
LYS 33 1RCQ.subA 9 0 1 1 0 0 0 0 0 8 1 0 0 0 3 1 4 0 2 2 357
LYS 149 1RCQ.subA 8 0 1 1 0 1 0 1 1 6 2 1 0 1 2 0 3 0 2 2 357
LYS 152 1RCQ.subA 3 0 0 0 0 0 0 0 0 3 0 0 1 1 1 0 1 0 6 2 357
LYS 201 1RCQ.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 0 1 0 7 2 357
LYS 238 1RCQ.subA 6 0 0 0 0 1 0 1 1 5 1 2 0 2 2 0 1 0 3 2 357
LYS 291 1RCQ.subA 5 1 0 1 0 0 1 0 1 3 2 0 0 0 2 0 1 0 5 2 357
LYS 348 1RCQ.subA 5 1 0 1 0 1 0 1 1 3 2 0 1 1 1 0 1 0 5 2 357
LYS 81 1RQB.subA 4 0 0 0 0 2 0 2 2 2 2 1 0 1 0 0 1 0 5 2 472
LYS 114 1RQB.subA 10 2 2 4 0 1 0 1 1 5 5 1 1 2 0 2 1 0 3 2 472
LYS 143 1RQB.subA 6 1 0 1 0 0 0 0 0 4 1 0 0 0 3 0 1 0 7 2 472
LYS 144 1RQB.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 3 0 0 0 9 2 472
LYS 147 1RQB.subA 8 1 0 1 0 0 1 0 1 5 2 0 0 0 3 1 1 0 3 2 472
LYS 168 1RQB.subA 5 0 1 1 0 0 0 0 0 4 1 0 1 1 1 0 2 0 6 2 472
LYS 191 1RQB.subA 8 0 0 0 0 0 0 0 0 8 0 1 1 2 2 1 3 0 3 2 472
LYS 200 1RQB.subA 7 2 1 3 0 0 0 0 0 4 3 0 0 0 2 1 1 0 10 2 472
LYS 203 1RQB.subA 12 2 0 2 0 0 0 0 0 10 2 1 2 3 3 0 4 0 4 2 472
LYS 209 1RQB.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 1 1 0 0 5 2 472
LYS 229 1RQB.subA 9 0 1 1 0 0 0 0 0 7 1 1 0 1 3 1 2 0 3 2 472
LYS 277 1RQB.subA 5 2 0 2 0 1 1 1 2 1 4 0 0 0 0 0 1 0 4 2 472
LYS 283 1RQB.subA 5 1 0 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 7 2 472
LYS 288 1RQB.subA 11 2 0 2 0 1 0 1 1 8 3 1 0 1 3 1 3 0 3 2 472
LYS 290 1RQB.subA 7 0 2 2 1 0 0 1 1 4 3 0 0 0 1 1 2 0 3 2 472
LYS 291 1RQB.subA 5 0 1 1 1 0 0 1 1 3 2 0 0 0 0 2 1 0 2 2 472
LYS 295 1RQB.subA 3 0 1 1 0 0 0 0 0 2 1 2 0 2 0 0 0 0 0 2 472
LYS 304 1RQB.subA 3 0 0 0 0 0 0 0 0 3 0 2 0 2 0 1 0 0 1 2 472
LYS 327 1RQB.subA 8 1 2 3 0 0 0 0 0 3 3 0 0 0 1 0 2 0 0 2 472
LYS 340 1RQB.subA 5 0 0 0 0 2 0 2 2 3 2 0 0 0 2 1 0 0 4 2 472
LYS 369 1RQB.subA 5 0 1 1 0 1 0 1 1 2 2 0 0 0 1 1 0 0 1 2 472
LYS 393 1RQB.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 1 2 472
LYS 396 1RQB.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 1 2 472
LYS 404 1RQB.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 1 1 1 0 0 2 472
LYS 405 1RQB.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 1 1 1 0 2 2 472
LYS 421 1RQB.subA 4 0 3 3 0 0 0 0 0 1 3 0 1 1 0 0 0 0 3 2 472
LYS 424 1RQB.subA 3 0 2 2 0 0 0 0 0 1 2 1 0 1 0 0 0 0 1 2 472
LYS 430 1RQB.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 3 2 472
LYS 471 1RQB.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 2 0 1 0 0 2 472
LYS 81 1RQH.subA 4 0 0 0 0 2 0 2 2 2 2 1 0 1 0 0 1 0 7 2 471
LYS 114 1RQH.subA 10 2 2 4 0 1 0 1 1 5 5 1 1 2 0 2 1 0 1 2 471
LYS 143 1RQH.subA 6 1 0 1 0 0 0 0 0 5 1 0 0 0 3 0 2 0 4 2 471
LYS 144 1RQH.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 3 0 0 0 8 2 471
LYS 147 1RQH.subA 8 1 0 1 0 0 1 0 1 6 2 0 0 0 3 1 2 0 6 2 471
LYS 168 1RQH.subA 5 0 1 1 0 0 0 0 0 4 1 0 1 1 1 0 2 0 5 2 471
LYS 184 1RQH.subA 12 1 0 1 0 0 1 0 1 10 2 1 0 1 2 1 5 1 2 2 471
LYS 191 1RQH.subA 8 0 0 0 0 0 0 0 0 8 0 1 1 2 2 1 3 0 2 2 471
LYS 200 1RQH.subA 7 2 1 3 0 0 0 0 0 4 3 0 0 0 2 1 1 0 7 2 471
LYS 203 1RQH.subA 12 2 0 2 0 0 0 0 0 10 2 1 2 3 3 0 4 0 4 2 471
LYS 209 1RQH.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 1 1 0 0 5 2 471
LYS 229 1RQH.subA 10 0 1 1 0 0 0 0 0 9 1 1 1 2 3 1 3 0 2 2 471
LYS 277 1RQH.subA 5 2 0 2 0 1 1 1 2 1 4 0 0 0 0 0 1 0 6 2 471
LYS 283 1RQH.subA 5 1 0 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 8 2 471
LYS 288 1RQH.subA 10 2 0 2 0 1 0 1 1 7 3 1 0 1 2 1 3 0 4 2 471
LYS 290 1RQH.subA 7 0 2 2 1 0 0 1 1 4 3 0 0 0 1 1 2 0 5 2 471
LYS 291 1RQH.subA 5 0 1 1 1 0 0 1 1 3 2 0 0 0 0 2 1 0 3 2 471
LYS 295 1RQH.subA 3 0 1 1 0 0 0 0 0 2 1 2 0 2 0 0 0 0 0 2 471
LYS 304 1RQH.subA 3 0 0 0 0 0 0 0 0 3 0 2 0 2 0 1 0 0 4 2 471
LYS 327 1RQH.subA 8 1 2 3 0 0 0 0 0 5 3 0 0 0 2 0 3 0 0 2 471
LYS 340 1RQH.subA 5 0 0 0 0 2 0 2 2 3 2 0 0 0 2 1 0 0 5 2 471
LYS 369 1RQH.subA 5 0 1 1 0 1 0 1 1 3 2 0 0 0 1 1 1 0 0 2 471
LYS 393 1RQH.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 0 2 471
LYS 396 1RQH.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 0 2 471
LYS 404 1RQH.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 1 1 1 0 1 2 471
LYS 405 1RQH.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 1 0 0 2 2 471
LYS 421 1RQH.subA 4 0 3 3 0 0 0 0 0 1 3 0 1 1 0 0 0 0 3 2 471
LYS 424 1RQH.subA 3 0 2 2 0 0 0 0 0 1 2 1 0 1 0 0 0 0 1 2 471
LYS 430 1RQH.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 1 2 471
LYS 471 1RQH.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 2 0 1 0 0 2 471
LYS 81 1RR2.subA 4 0 0 0 0 2 0 2 2 2 2 1 0 1 0 0 1 0 3 2 472
LYS 114 1RR2.subA 10 2 2 4 0 1 0 1 1 5 5 1 1 2 0 2 1 0 0 2 472
LYS 143 1RR2.subA 7 1 0 1 0 0 0 0 0 6 1 0 0 0 3 0 3 0 4 2 472
LYS 144 1RR2.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 3 0 0 0 4 2 472
LYS 147 1RR2.subA 8 1 0 1 0 0 1 0 1 6 2 0 0 0 3 1 2 0 3 2 472
LYS 168 1RR2.subA 6 0 1 1 0 0 0 0 0 5 1 0 1 1 1 1 2 0 2 2 472
LYS 184 1RR2.subA 13 1 0 1 0 1 1 1 2 10 3 1 0 1 2 1 5 0 2 2 472
LYS 191 1RR2.subA 8 0 0 0 0 0 0 0 0 8 0 1 1 2 2 1 3 0 2 2 472
LYS 200 1RR2.subA 7 2 1 3 0 0 0 0 0 4 3 0 0 0 2 1 1 0 6 2 472
LYS 203 1RR2.subA 12 2 0 2 0 0 0 0 0 10 2 1 2 3 3 0 4 0 2 2 472
LYS 209 1RR2.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 1 1 0 0 3 2 472
LYS 229 1RR2.subA 10 0 1 1 0 0 0 0 0 9 1 1 1 2 3 1 3 0 1 2 472
LYS 277 1RR2.subA 5 2 0 2 0 1 1 1 2 1 4 0 0 0 0 0 1 0 4 2 472
LYS 283 1RR2.subA 5 1 0 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 5 2 472
LYS 288 1RR2.subA 10 2 0 2 0 1 0 1 1 7 3 1 0 1 2 1 3 0 4 2 472
LYS 290 1RR2.subA 7 0 2 2 1 0 0 1 1 4 3 0 0 0 1 1 2 0 2 2 472
LYS 291 1RR2.subA 5 0 1 1 1 0 0 1 1 3 2 0 0 0 0 2 1 0 2 2 472
LYS 295 1RR2.subA 3 0 1 1 0 0 0 0 0 2 1 2 0 2 0 0 0 0 0 2 472
LYS 304 1RR2.subA 3 0 0 0 0 0 0 0 0 3 0 2 0 2 0 1 0 0 2 2 472
LYS 327 1RR2.subA 8 1 2 3 0 0 0 0 0 5 3 0 0 0 1 0 4 0 0 2 472
LYS 340 1RR2.subA 6 0 0 0 0 2 0 2 2 4 2 0 0 0 2 1 1 0 3 2 472
LYS 369 1RR2.subA 5 0 1 1 0 1 0 1 1 3 2 0 0 0 1 1 1 0 0 2 472
LYS 393 1RR2.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 0 2 472
LYS 396 1RR2.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 0 2 0 0 2 472
LYS 404 1RR2.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 1 1 1 0 0 2 472
LYS 405 1RR2.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 2 1 0 0 0 2 472
LYS 421 1RR2.subA 5 0 3 3 0 0 0 0 0 2 3 0 1 1 0 1 0 0 3 2 472
LYS 424 1RR2.subA 3 0 2 2 0 0 0 0 0 1 2 1 0 1 0 0 0 0 0 2 472
LYS 430 1RR2.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 1 2 472
LYS 471 1RR2.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 2 0 1 0 0 2 472
LYS 14 1RXO.subB 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 0 2 455
LYS 18 1RXO.subB 4 1 0 1 0 0 0 0 0 3 1 1 0 1 1 0 1 0 0 2 455
LYS 21 1RXO.subB 7 1 2 3 0 0 0 0 0 4 3 0 0 0 2 1 1 0 1 2 455
LYS 81 1RXO.subB 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 1 2 455
LYS 128 1RXO.subB 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 1 2 455
LYS 146 1RXO.subB 5 0 0 0 0 0 0 0 0 5 0 1 0 1 1 0 3 0 4 2 455
LYS 161 1RXO.subB 6 1 1 2 0 0 1 0 1 3 3 0 0 0 1 0 2 0 3 2 455
LYS 164 1RXO.subB 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 2 2 455
LYS 175 1RXO.subB 7 1 0 1 0 0 0 0 0 6 1 1 0 1 3 0 2 0 3 2 455
LYS 177 1RXO.subB 7 1 1 2 0 0 0 0 0 5 2 0 1 1 1 0 3 0 10 2 455
LYS 183 1RXO.subB 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 4 2 455
LYS 227 1RXO.subB 5 0 1 1 0 0 0 0 0 4 1 0 0 0 2 2 0 0 4 2 455
LYS 236 1RXO.subB 12 1 1 2 0 1 0 1 1 9 3 1 1 2 3 0 4 0 2 2 455
LYS 252 1RXO.subB 5 1 1 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 3 2 455
LYS 305 1RXO.subB 7 0 2 2 0 2 0 2 2 3 4 0 3 3 0 0 0 0 1 2 455
LYS 316 1RXO.subB 11 1 0 1 0 2 0 2 2 8 3 0 1 1 2 1 4 0 1 2 455
LYS 334 1RXO.subB 8 0 1 1 0 1 0 1 1 6 2 2 0 2 1 0 3 0 1 2 455
LYS 356 1RXO.subB 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 2 2 455
LYS 450 1RXO.subB 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 0 2 455
LYS 463 1RXO.subB 4 0 0 0 0 0 1 0 1 3 1 0 0 0 0 2 1 0 1 2 455
LYS 81 1S3H.subA 4 0 0 0 0 2 0 2 2 2 2 1 0 1 0 0 1 0 0 2 472
LYS 114 1S3H.subA 10 2 2 4 0 1 0 1 1 5 5 1 1 2 0 2 1 0 0 2 472
LYS 143 1S3H.subA 6 1 0 1 0 0 0 0 0 5 1 0 0 0 3 0 2 0 0 2 472
LYS 144 1S3H.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 3 0 0 0 0 2 472
LYS 147 1S3H.subA 8 1 0 1 0 0 1 0 1 6 2 0 0 0 3 1 2 0 0 2 472
LYS 168 1S3H.subA 6 0 1 1 0 0 0 0 0 5 1 0 1 1 1 1 2 0 0 2 472
LYS 191 1S3H.subA 7 0 0 0 0 0 0 0 0 7 0 0 1 1 2 1 3 0 0 2 472
LYS 200 1S3H.subA 7 2 1 3 0 0 0 0 0 4 3 0 0 0 2 1 1 0 0 2 472
LYS 203 1S3H.subA 12 2 0 2 0 0 0 0 0 10 2 1 2 3 3 0 4 0 0 2 472
LYS 209 1S3H.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 1 1 0 0 0 2 472
LYS 229 1S3H.subA 10 0 1 1 0 0 0 0 0 9 1 1 1 2 3 1 3 0 0 2 472
LYS 277 1S3H.subA 5 2 0 2 0 1 1 1 2 1 4 0 0 0 0 0 1 0 0 2 472
LYS 283 1S3H.subA 5 1 0 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 0 2 472
LYS 288 1S3H.subA 10 2 0 2 0 1 0 1 1 7 3 1 0 1 2 1 3 0 0 2 472
LYS 290 1S3H.subA 7 0 2 2 1 0 0 1 1 4 3 0 0 0 1 1 2 0 0 2 472
LYS 291 1S3H.subA 5 0 1 1 1 0 0 1 1 3 2 0 0 0 0 2 1 0 0 2 472
LYS 295 1S3H.subA 3 0 1 1 0 0 0 0 0 2 1 2 0 2 0 0 0 0 0 2 472
LYS 304 1S3H.subA 3 0 0 0 0 0 0 0 0 3 0 2 0 2 0 1 0 0 0 2 472
LYS 327 1S3H.subA 8 1 2 3 0 0 0 0 0 5 3 0 0 0 1 0 4 0 0 2 472
LYS 340 1S3H.subA 5 0 0 0 0 2 0 2 2 3 2 0 0 0 2 1 0 0 0 2 472
LYS 369 1S3H.subA 5 0 1 1 0 1 0 1 1 3 2 0 0 0 1 1 1 0 0 2 472
LYS 393 1S3H.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 0 2 472
LYS 396 1S3H.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 0 2 472
LYS 404 1S3H.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 1 1 1 0 0 2 472
LYS 405 1S3H.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 0 2 472
LYS 421 1S3H.subA 4 0 3 3 0 0 0 0 0 1 3 0 1 1 0 0 0 0 0 2 472
LYS 424 1S3H.subA 3 0 2 2 0 0 0 0 0 1 2 1 0 1 0 0 0 0 0 2 472
LYS 430 1S3H.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 0 2 472
LYS 471 1S3H.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 2 0 1 0 0 2 472
LYS 2 1S3T.subC 2 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 2 0 0 2 570
LYS 33 1S3T.subC 4 1 1 2 0 0 0 0 0 2 2 0 0 0 0 2 0 0 5 2 570
LYS 48 1S3T.subC 7 0 0 0 0 0 0 0 0 7 0 0 1 1 3 1 2 0 4 2 570
LYS 84 1S3T.subC 13 1 0 1 0 0 0 0 0 12 1 1 1 2 4 4 2 0 4 2 570
LYS 90 1S3T.subC 11 2 1 3 0 1 0 1 1 7 4 1 1 2 0 1 4 0 4 2 570
LYS 99 1S3T.subC 15 3 1 4 0 0 0 0 0 11 4 3 0 3 4 1 3 0 5 2 570
LYS 127 1S3T.subC 11 1 2 3 0 1 0 1 1 7 4 0 0 0 3 1 3 0 1 2 570
LYS 169 1S3T.subC 7 0 1 1 0 0 0 0 0 6 1 1 0 1 4 0 1 0 5 2 570
LYS 182 1S3T.subC 5 0 1 1 1 0 0 1 1 3 2 0 1 1 0 1 1 0 1 2 570
LYS 185 1S3T.subC 6 0 1 1 1 0 0 1 1 4 2 1 0 1 0 0 3 0 3 2 570
LYS 199 1S3T.subC 7 0 0 0 0 0 1 0 1 6 1 0 1 1 3 1 1 0 0 2 570
LYS 326 1S3T.subC 6 1 0 1 0 0 1 0 1 4 2 0 2 2 0 0 2 0 2 2 570
LYS 383 1S3T.subC 12 2 1 3 0 0 0 0 0 9 3 1 2 3 3 1 2 0 1 2 570
LYS 385 1S3T.subC 10 3 0 3 0 1 0 1 1 6 4 0 0 0 3 1 2 0 2 2 570
LYS 386 1S3T.subC 2 1 0 1 0 0 0 0 0 1 1 0 1 1 0 0 0 0 1 2 570
LYS 395 1S3T.subC 3 0 2 2 0 0 0 0 0 1 2 0 1 1 0 0 0 0 0 2 570
LYS 404 1S3T.subC 11 1 1 2 0 1 0 1 1 8 3 0 1 1 1 3 3 0 2 2 570
LYS 409 1S3T.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 1 2 4 0 0 2 570
LYS 431 1S3T.subC 12 0 3 3 0 0 1 0 1 8 4 1 0 1 3 1 3 0 2 2 570
LYS 441 1S3T.subC 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 1 0 0 4 2 570
LYS 446 1S3T.subC 7 0 0 0 0 0 0 0 0 7 0 0 1 1 2 3 1 0 2 2 570
LYS 452 1S3T.subC 12 1 0 1 0 0 0 0 0 11 1 2 2 4 2 1 4 0 1 2 570
LYS 497 1S3T.subC 8 1 0 1 1 0 0 1 1 6 2 3 1 4 1 0 1 0 4 2 570
LYS 507 1S3T.subC 5 0 1 1 0 0 0 0 0 4 1 0 0 0 2 0 2 0 3 2 570
LYS 511 1S3T.subC 5 0 1 1 0 1 1 1 2 2 3 0 0 0 1 0 1 0 2 2 570
LYS 518 1S3T.subC 7 1 0 1 1 0 1 1 2 4 3 1 1 2 1 0 1 0 3 2 570
LYS 525 1S3T.subC 16 1 0 1 0 0 1 0 1 14 2 2 1 3 6 0 5 0 2 2 570
LYS 526 1S3T.subC 2 1 0 1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 2 2 570
LYS 529 1S3T.subC 8 1 0 1 0 0 0 0 0 7 1 1 1 2 0 2 3 0 4 2 570
LYS 547 1S3T.subC 6 1 1 2 0 0 0 0 0 4 2 0 0 0 1 0 3 0 2 2 570
LYS 559 1S3T.subC 3 0 1 1 0 0 0 0 0 2 1 0 0 0 0 0 2 0 2 2 570
LYS 81 1U5J.subA 4 0 0 0 0 2 0 2 2 2 2 1 0 1 0 0 1 0 0 2 471
LYS 114 1U5J.subA 10 2 2 4 0 1 0 1 1 5 5 1 1 2 0 2 1 0 0 2 471
LYS 143 1U5J.subA 6 1 0 1 0 0 0 0 0 5 1 0 0 0 3 0 2 0 0 2 471
LYS 144 1U5J.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 3 0 0 0 0 2 471
LYS 147 1U5J.subA 8 1 0 1 0 0 1 0 1 6 2 0 0 0 3 1 2 0 1 2 471
LYS 168 1U5J.subA 5 0 1 1 0 0 0 0 0 4 1 0 1 1 1 0 2 0 1 2 471
LYS 191 1U5J.subA 8 0 0 0 0 0 0 0 0 8 0 1 1 2 2 1 3 0 0 2 471
LYS 200 1U5J.subA 7 2 1 3 0 0 0 0 0 4 3 0 0 0 2 1 1 0 2 2 471
LYS 203 1U5J.subA 12 2 0 2 0 0 0 0 0 10 2 1 2 3 3 0 4 0 0 2 471
LYS 209 1U5J.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 1 1 0 0 0 2 471
LYS 229 1U5J.subA 9 0 1 1 0 0 0 0 0 8 1 1 1 2 2 1 3 0 0 2 471
LYS 277 1U5J.subA 5 2 0 2 0 1 1 1 2 1 4 0 0 0 0 0 1 0 0 2 471
LYS 283 1U5J.subA 6 1 0 1 0 1 1 1 2 3 3 1 0 1 1 1 0 0 0 2 471
LYS 288 1U5J.subA 9 2 0 2 0 0 0 0 0 7 2 1 0 1 2 1 3 0 0 2 471
LYS 290 1U5J.subA 7 0 2 2 1 0 0 1 1 4 3 0 0 0 1 1 2 0 0 2 471
LYS 291 1U5J.subA 4 0 1 1 1 0 0 1 1 2 2 0 0 0 0 1 1 0 0 2 471
LYS 295 1U5J.subA 3 0 1 1 0 0 0 0 0 2 1 2 0 2 0 0 0 0 0 2 471
LYS 304 1U5J.subA 3 0 0 0 0 0 0 0 0 3 0 2 0 2 0 1 0 0 0 2 471
LYS 327 1U5J.subA 6 1 1 2 0 0 0 0 0 4 2 0 0 0 1 0 3 0 0 2 471
LYS 340 1U5J.subA 5 0 0 0 0 2 0 2 2 3 2 0 0 0 2 1 0 0 0 2 471
LYS 369 1U5J.subA 5 0 1 1 0 1 0 1 1 3 2 0 0 0 1 1 1 0 0 2 471
LYS 393 1U5J.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 0 2 471
LYS 396 1U5J.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 1 0 3 0 0 2 471
LYS 404 1U5J.subA 6 1 0 1 1 0 0 1 1 4 2 1 0 1 1 1 1 0 0 2 471
LYS 405 1U5J.subA 5 0 0 0 1 0 0 1 1 4 1 0 0 0 1 1 2 0 0 2 471
LYS 421 1U5J.subA 5 0 3 3 0 0 0 0 0 2 3 0 1 1 0 1 0 0 0 2 471
LYS 424 1U5J.subA 4 0 2 2 0 0 0 0 0 2 2 1 1 2 0 0 0 0 0 2 471
LYS 430 1U5J.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 0 2 471
LYS 471 1U5J.subA 5 1 0 1 0 0 0 0 0 4 1 0 1 1 2 0 1 0 0 2 471
LYS 6 1UAG.subA 9 2 1 3 0 0 0 0 0 6 3 0 2 2 1 1 2 0 2 2 428
LYS 45 1UAG.subA 5 1 0 1 0 0 0 0 0 4 1 0 0 0 2 0 2 0 0 2 428
LYS 115 1UAG.subA 13 0 1 1 0 0 0 0 0 12 1 4 2 6 2 0 3 0 3 2 428
LYS 127 1UAG.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 3 0 2 0 5 2 428
LYS 206 1UAG.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 0 2 0 3 2 428
LYS 251 1UAG.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 1 1 0 0 2 428
LYS 254 1UAG.subA 6 0 1 1 0 1 0 1 1 4 2 0 1 1 0 1 2 0 2 2 428
LYS 259 1UAG.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 3 2 428
LYS 262 1UAG.subA 7 0 1 1 0 0 0 0 0 6 1 3 0 3 0 1 2 0 1 2 428
LYS 291 1UAG.subA 8 0 0 0 0 0 0 0 0 8 0 4 0 4 3 0 1 0 1 2 428
LYS 319 1UAG.subA 9 1 0 1 0 1 1 1 2 6 3 2 0 2 2 1 1 0 5 2 428
LYS 348 1UAG.subA 12 3 0 3 0 0 1 0 1 8 4 2 1 3 4 0 1 0 7 2 428
LYS 420 1UAG.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 1 0 0 1 2 428
LYS 434 1UAG.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 0 2 428
LYS 2 1UBP.subC 2 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 2 0 1 2 570
LYS 33 1UBP.subC 4 1 1 2 0 0 0 0 0 2 2 0 0 0 0 2 0 0 15 2 570
LYS 48 1UBP.subC 7 0 0 0 0 0 0 0 0 7 0 0 1 1 3 1 2 0 10 2 570
LYS 84 1UBP.subC 12 1 0 1 0 0 0 0 0 11 1 1 1 2 4 3 2 0 4 2 570
LYS 90 1UBP.subC 11 2 1 3 0 1 0 1 1 7 4 1 1 2 0 1 4 0 10 2 570
LYS 99 1UBP.subC 15 3 1 4 0 0 0 0 0 11 4 3 0 3 4 1 3 0 6 2 570
LYS 127 1UBP.subC 11 1 2 3 0 1 0 1 1 7 4 0 0 0 3 1 3 0 5 2 570
LYS 169 1UBP.subC 7 0 1 1 0 0 0 0 0 6 1 1 0 1 4 0 1 0 10 2 570
LYS 182 1UBP.subC 5 0 1 1 1 0 0 1 1 3 2 0 1 1 0 1 1 0 4 2 570
LYS 185 1UBP.subC 6 0 1 1 1 0 0 1 1 4 2 1 0 1 0 0 3 0 5 2 570
LYS 199 1UBP.subC 8 0 0 0 0 0 1 0 1 7 1 0 1 1 3 1 1 0 0 2 570
LYS 326 1UBP.subC 6 1 0 1 0 0 1 0 1 4 2 0 2 2 0 0 2 0 6 2 570
LYS 383 1UBP.subC 12 2 1 3 0 0 0 0 0 9 3 1 2 3 3 1 2 0 2 2 570
LYS 385 1UBP.subC 10 3 0 3 0 1 0 1 1 6 4 0 0 0 3 1 2 0 8 2 570
LYS 386 1UBP.subC 5 2 1 3 0 0 0 0 0 2 3 0 1 1 0 0 1 0 8 2 570
LYS 395 1UBP.subC 4 0 3 3 0 0 0 0 0 1 3 0 1 1 0 0 0 0 5 2 570
LYS 404 1UBP.subC 11 1 1 2 0 1 0 1 1 8 3 0 1 1 1 3 3 0 4 2 570
LYS 409 1UBP.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 1 2 4 0 0 2 570
LYS 431 1UBP.subC 12 0 3 3 0 0 1 0 1 8 4 1 0 1 3 1 3 0 7 2 570
LYS 441 1UBP.subC 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 1 0 0 4 2 570
LYS 446 1UBP.subC 7 0 0 0 0 0 0 0 0 7 0 0 1 1 2 3 1 0 3 2 570
LYS 452 1UBP.subC 13 1 0 1 0 0 0 0 0 12 1 2 2 4 3 1 4 0 3 2 570
LYS 497 1UBP.subC 8 1 0 1 1 0 0 1 1 6 2 3 1 4 1 0 1 0 12 2 570
LYS 507 1UBP.subC 5 0 1 1 0 0 0 0 0 4 1 0 0 0 2 0 2 0 13 2 570
LYS 511 1UBP.subC 5 0 1 1 0 1 1 1 2 2 3 0 0 0 1 0 1 0 10 2 570
LYS 518 1UBP.subC 7 1 0 1 1 0 1 1 2 4 3 1 1 2 1 0 1 0 12 2 570
LYS 525 1UBP.subC 16 1 0 1 0 0 1 0 1 14 2 2 1 3 6 0 5 0 4 2 570
LYS 526 1UBP.subC 2 1 0 1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 8 2 570
LYS 529 1UBP.subC 8 1 0 1 0 0 0 0 0 7 1 1 1 2 0 2 3 0 9 2 570
LYS 547 1UBP.subC 7 1 2 3 0 0 0 0 0 4 3 0 0 0 1 0 3 0 5 2 570
LYS 559 1UBP.subC 3 0 1 1 0 0 0 0 0 2 1 0 0 0 0 0 2 0 7 2 570
LYS 14 1UPM.subB 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 1 0 0 5 2 467
LYS 18 1UPM.subB 5 1 0 1 0 0 0 0 0 4 1 2 0 2 0 0 2 0 7 2 467
LYS 21 1UPM.subB 8 1 2 3 0 0 0 0 0 5 3 0 0 0 2 2 1 0 1 2 467
LYS 81 1UPM.subB 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 2 2 467
LYS 128 1UPM.subB 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 1 0 0 1 2 467
LYS 146 1UPM.subB 5 0 0 0 0 0 0 0 0 5 0 1 0 1 1 0 3 0 4 2 467
LYS 161 1UPM.subB 6 1 1 2 0 0 1 0 1 3 3 0 0 0 1 0 2 0 4 2 467
LYS 164 1UPM.subB 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 3 2 467
LYS 175 1UPM.subB 8 1 0 1 1 0 0 1 1 6 2 1 0 1 3 0 2 0 0 2 467
LYS 177 1UPM.subB 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 0 1 2 467
LYS 183 1UPM.subB 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 2 2 467
LYS 227 1UPM.subB 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 2 0 0 3 2 467
LYS 236 1UPM.subB 12 1 1 2 0 1 0 1 1 9 3 1 1 2 3 0 4 0 1 2 467
LYS 252 1UPM.subB 5 1 1 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 4 2 467
LYS 305 1UPM.subB 8 1 0 1 0 2 0 2 2 5 3 1 2 3 0 0 2 0 4 2 467
LYS 316 1UPM.subB 11 1 0 1 0 2 0 2 2 8 3 0 1 1 2 1 4 0 1 2 467
LYS 334 1UPM.subB 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 3 2 467
LYS 356 1UPM.subB 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 4 2 467
LYS 450 1UPM.subB 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 0 2 467
LYS 463 1UPM.subB 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 2 0 0 3 2 467
LYS 466 1UPM.subB 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 3 2 467
LYS 14 1UPP.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 1 0 0 0 2 467
LYS 18 1UPP.subA 5 1 0 1 0 0 0 0 0 4 1 2 0 2 0 0 2 0 3 2 467
LYS 21 1UPP.subA 7 1 2 3 0 0 0 0 0 4 3 0 0 0 1 2 1 0 0 2 467
LYS 81 1UPP.subA 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 0 2 467
LYS 128 1UPP.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 1 0 0 1 2 467
LYS 146 1UPP.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 1 0 3 0 2 2 467
LYS 161 1UPP.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 1 0 2 0 0 2 467
LYS 164 1UPP.subA 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 0 2 467
LYS 175 1UPP.subA 8 1 0 1 1 0 0 1 1 6 2 1 0 1 3 0 2 0 1 2 467
LYS 177 1UPP.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 0 2 2 467
LYS 183 1UPP.subA 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 0 2 467
LYS 227 1UPP.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 2 0 0 2 2 467
LYS 236 1UPP.subA 12 1 1 2 0 1 0 1 1 9 3 1 1 2 3 0 4 0 1 2 467
LYS 252 1UPP.subA 5 1 1 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 5 2 467
LYS 305 1UPP.subA 9 1 1 2 0 2 0 2 2 5 4 1 2 3 0 0 2 0 1 2 467
LYS 316 1UPP.subA 12 1 0 1 0 2 0 2 2 9 3 0 1 1 2 2 4 0 0 2 467
LYS 334 1UPP.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 2 2 467
LYS 356 1UPP.subA 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 2 2 467
LYS 450 1UPP.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 0 1 0 0 1 2 467
LYS 463 1UPP.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 2 0 0 0 2 467
LYS 466 1UPP.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 2 2 467
LYS 14 1UW9.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 5 2 465
LYS 18 1UW9.subA 5 1 0 1 0 0 0 0 0 4 1 2 0 2 0 0 2 0 5 2 465
LYS 81 1UW9.subA 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 3 2 465
LYS 128 1UW9.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 0 2 465
LYS 146 1UW9.subA 5 0 0 0 0 0 0 0 0 4 0 1 0 1 2 0 1 0 1 2 465
LYS 161 1UW9.subA 7 1 1 2 0 0 1 0 1 4 3 0 0 0 1 1 2 0 3 2 465
LYS 164 1UW9.subA 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 3 2 465
LYS 175 1UW9.subA 8 1 0 1 1 0 0 1 1 6 2 1 0 1 3 0 2 1 0 2 465
LYS 177 1UW9.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 1 2 465
LYS 183 1UW9.subA 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 2 2 465
LYS 227 1UW9.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 2 0 0 3 2 465
LYS 236 1UW9.subA 12 1 1 2 0 1 0 1 1 9 3 1 1 2 3 0 4 0 2 2 465
LYS 252 1UW9.subA 5 0 2 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 6 2 465
LYS 258 1UW9.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 2 1 2 0 4 2 465
LYS 316 1UW9.subA 13 1 0 1 0 2 0 2 2 10 3 0 1 1 2 2 5 0 1 2 465
LYS 334 1UW9.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 3 2 465
LYS 356 1UW9.subA 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 3 2 465
LYS 450 1UW9.subA 4 0 0 0 0 1 0 1 1 3 1 1 0 1 1 1 0 0 1 2 465
LYS 463 1UW9.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 2 0 0 4 2 465
LYS 466 1UW9.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 1 2 465
LYS 474 1UW9.subA 6 2 0 2 0 1 0 1 1 3 3 1 0 1 0 0 2 0 3 2 465
LYS 14 1UWA.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 5 2 465
LYS 18 1UWA.subA 5 1 0 1 0 0 0 0 0 4 1 2 0 2 0 0 2 0 5 2 465
LYS 81 1UWA.subA 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 2 2 465
LYS 128 1UWA.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 1 0 0 1 2 465
LYS 146 1UWA.subA 5 0 0 0 0 0 0 0 0 4 0 1 0 1 2 0 1 0 3 2 465
LYS 161 1UWA.subA 7 1 1 2 0 0 1 0 1 4 3 0 0 0 1 1 2 0 2 2 465
LYS 164 1UWA.subA 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 4 2 465
LYS 175 1UWA.subA 8 1 0 1 1 0 0 1 1 6 2 1 0 1 3 0 2 1 0 2 465
LYS 177 1UWA.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 1 2 465
LYS 183 1UWA.subA 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 2 2 465
LYS 227 1UWA.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 2 0 0 2 2 465
LYS 236 1UWA.subA 12 1 1 2 0 1 0 1 1 9 3 1 1 2 3 0 4 0 3 2 465
LYS 252 1UWA.subA 5 0 2 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 5 2 465
LYS 258 1UWA.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 2 1 2 0 4 2 465
LYS 316 1UWA.subA 12 1 0 1 0 2 0 2 2 9 3 0 1 1 2 2 4 0 1 2 465
LYS 334 1UWA.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 2 2 465
LYS 356 1UWA.subA 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 5 2 465
LYS 450 1UWA.subA 4 0 0 0 0 1 0 1 1 3 1 1 0 1 1 1 0 0 2 2 465
LYS 463 1UWA.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 2 0 0 2 2 465
LYS 466 1UWA.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 4 2 465
LYS 474 1UWA.subA 6 2 0 2 0 1 0 1 1 3 3 1 0 1 0 0 2 0 1 2 465
LYS 14 1UZD.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 1 0 0 1 2 465
LYS 18 1UZD.subA 5 1 0 1 0 0 0 0 0 4 1 2 0 2 0 0 2 0 2 2 465
LYS 81 1UZD.subA 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 1 2 465
LYS 128 1UZD.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 1 0 0 0 2 465
LYS 146 1UZD.subA 6 0 0 0 0 0 0 0 0 5 0 1 0 1 2 1 1 0 0 2 465
LYS 161 1UZD.subA 8 1 1 2 0 0 1 0 1 5 3 0 0 0 1 0 4 0 2 2 465
LYS 164 1UZD.subA 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 2 2 465
LYS 175 1UZD.subA 8 1 0 1 1 0 0 1 1 6 2 1 0 1 3 0 2 1 0 2 465
LYS 177 1UZD.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 0 2 465
LYS 183 1UZD.subA 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 2 2 465
LYS 227 1UZD.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 2 0 0 1 2 465
LYS 236 1UZD.subA 13 1 1 2 0 2 0 2 2 9 4 1 1 2 3 0 4 0 3 2 465
LYS 252 1UZD.subA 5 0 2 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 2 2 465
LYS 258 1UZD.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 2 1 2 0 3 2 465
LYS 316 1UZD.subA 12 1 0 1 0 2 0 2 2 9 3 0 1 1 2 2 4 0 1 2 465
LYS 334 1UZD.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 1 2 465
LYS 356 1UZD.subA 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 3 2 465
LYS 450 1UZD.subA 4 0 0 0 0 1 0 1 1 3 1 1 0 1 1 1 0 0 1 2 465
LYS 463 1UZD.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 2 0 0 2 2 465
LYS 466 1UZD.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 2 2 465
LYS 474 1UZD.subA 6 2 0 2 0 1 0 1 1 3 3 1 0 1 0 0 2 0 1 2 465
LYS 14 1UZH.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 3 2 465
LYS 18 1UZH.subA 5 1 0 1 0 0 0 0 0 4 1 2 0 2 0 0 2 0 6 2 465
LYS 81 1UZH.subA 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 3 2 465
LYS 128 1UZH.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 1 2 465
LYS 146 1UZH.subA 5 0 0 0 0 0 0 0 0 4 0 1 0 1 2 0 1 0 2 2 465
LYS 161 1UZH.subA 8 1 1 2 0 0 1 0 1 5 3 0 0 0 1 0 4 0 5 2 465
LYS 164 1UZH.subA 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 3 2 465
LYS 175 1UZH.subA 8 1 0 1 1 0 0 1 1 6 2 1 0 1 3 0 2 1 1 2 465
LYS 177 1UZH.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 2 2 465
LYS 183 1UZH.subA 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 4 2 465
LYS 227 1UZH.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 2 0 0 2 2 465
LYS 236 1UZH.subA 13 1 1 2 0 2 0 2 2 9 4 1 1 2 3 0 4 0 3 2 465
LYS 252 1UZH.subA 5 0 2 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 5 2 465
LYS 258 1UZH.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 2 1 2 0 4 2 465
LYS 316 1UZH.subA 13 1 0 1 0 2 0 2 2 10 3 0 1 1 2 2 5 0 1 2 465
LYS 334 1UZH.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 4 2 465
LYS 356 1UZH.subA 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 3 2 465
LYS 450 1UZH.subA 4 0 0 0 0 1 0 1 1 3 1 1 0 1 1 1 0 0 0 2 465
LYS 463 1UZH.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 2 0 0 1 2 465
LYS 466 1UZH.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 2 2 465
LYS 474 1UZH.subA 6 2 0 2 0 1 0 1 1 3 3 1 0 1 0 0 2 0 5 2 465
LYS 38 1VFH.subA 9 0 1 1 0 0 0 0 0 8 1 2 0 2 2 2 2 0 1 2 382
LYS 196 1VFH.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 1 1 0 0 0 2 382
LYS 260 1VFH.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 0 0 4 0 0 2 382
LYS 308 1VFH.subA 6 1 0 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 0 2 382
LYS 38 1VFS.subA 9 0 1 1 0 0 0 0 0 8 1 2 0 2 2 2 2 0 3 2 383
LYS 196 1VFS.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 1 1 0 0 0 2 383
LYS 260 1VFS.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 0 0 4 0 0 2 383
LYS 308 1VFS.subA 6 1 0 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 0 2 383
LYS 38 1VFT.subA 9 0 1 1 0 0 0 0 0 8 1 2 0 2 2 2 2 0 3 2 382
LYS 196 1VFT.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 1 1 0 0 0 2 382
LYS 260 1VFT.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 0 0 4 0 0 2 382
LYS 308 1VFT.subA 6 1 0 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 0 2 382
LYS 26 1W78.subA 6 1 0 1 0 0 1 0 1 4 2 2 0 2 0 1 1 0 3 2 414
LYS 46 1W78.subA 6 1 0 1 1 0 0 1 1 4 2 0 0 0 3 0 1 0 12 2 414
LYS 60 1W78.subA 12 0 1 1 0 0 0 0 0 11 1 2 1 3 5 0 2 2 5 2 414
LYS 77 1W78.subA 6 1 0 1 0 0 0 0 0 5 1 0 2 2 0 1 2 0 4 2 414
LYS 136 1W78.subA 8 0 0 0 1 0 0 1 1 7 1 0 1 1 2 2 2 0 8 2 414
LYS 196 1W78.subA 7 1 1 2 0 1 0 1 1 4 3 1 0 1 2 0 1 0 12 2 414
LYS 216 1W78.subA 11 0 1 1 0 1 0 1 1 9 2 1 0 1 5 1 2 0 7 2 414
LYS 318 1W78.subA 6 0 0 0 0 0 0 0 0 6 0 2 0 2 2 0 2 0 3 2 414
LYS 322 1W78.subA 3 0 0 0 0 0 0 0 0 3 0 0 1 1 1 0 1 0 3 2 414
LYS 336 1W78.subA 8 3 0 3 0 0 2 0 2 3 5 0 0 0 0 0 3 1 8 2 414
LYS 346 1W78.subA 8 0 0 0 0 0 1 0 1 7 1 1 0 1 1 2 3 0 6 2 414
LYS 376 1W78.subA 7 2 0 2 0 0 0 0 0 5 2 1 0 1 1 3 0 0 3 2 414
LYS 393 1W78.subA 5 2 1 3 0 0 0 0 0 2 3 0 0 0 2 0 0 0 5 2 414
LYS 46 1W7K.subA 5 1 0 1 0 0 0 0 0 4 1 0 0 0 3 0 1 0 10 2 411
LYS 60 1W7K.subA 12 0 1 1 0 0 0 0 0 11 1 2 1 3 5 0 2 2 5 2 411
LYS 77 1W7K.subA 6 1 0 1 0 0 0 0 0 5 1 0 2 2 0 1 2 0 2 2 411
LYS 136 1W7K.subA 7 0 0 0 0 0 0 0 0 7 0 0 1 1 2 2 2 0 5 2 411
LYS 196 1W7K.subA 6 1 1 2 0 1 0 1 1 3 3 1 0 1 2 0 0 0 10 2 411
LYS 216 1W7K.subA 11 0 1 1 0 1 0 1 1 9 2 1 0 1 5 1 2 0 7 2 411
LYS 318 1W7K.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 2 0 2 0 3 2 411
LYS 322 1W7K.subA 2 0 0 0 0 0 0 0 0 2 0 0 1 1 1 0 0 0 2 2 411
LYS 336 1W7K.subA 9 3 0 3 0 0 2 0 2 4 5 0 1 1 0 0 3 1 9 2 411
LYS 346 1W7K.subA 7 0 0 0 0 0 1 0 1 6 1 1 0 1 1 1 3 0 10 2 411
LYS 376 1W7K.subA 7 2 0 2 0 0 0 0 0 5 2 1 0 1 1 3 0 0 3 2 411
LYS 393 1W7K.subA 5 2 1 3 0 0 0 0 0 2 3 0 0 0 2 0 0 0 5 2 411
LYS 14 1WDD.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 5 2 465
LYS 18 1WDD.subA 5 1 0 1 0 0 0 0 0 4 1 2 0 2 0 0 2 0 7 2 465
LYS 21 1WDD.subA 8 1 2 3 0 0 0 0 0 5 3 0 0 0 2 2 1 0 5 2 465
LYS 32 1WDD.subA 6 2 1 3 0 0 0 0 0 3 3 2 0 2 0 1 0 0 10 2 465
LYS 81 1WDD.subA 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 3 2 465
LYS 128 1WDD.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 2 2 465
LYS 146 1WDD.subA 5 0 0 0 0 0 0 0 0 5 0 3 0 3 1 0 1 0 3 2 465
LYS 161 1WDD.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 1 0 2 0 3 2 465
LYS 164 1WDD.subA 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 3 2 465
LYS 175 1WDD.subA 8 1 0 1 1 0 0 1 1 6 2 1 0 1 3 0 2 1 1 2 465
LYS 177 1WDD.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 3 2 465
LYS 183 1WDD.subA 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 0 2 465
LYS 227 1WDD.subA 6 0 2 2 0 0 0 0 0 4 2 1 0 1 1 2 0 0 6 2 465
LYS 236 1WDD.subA 12 1 1 2 0 1 0 1 1 9 3 2 1 3 2 0 4 0 1 2 465
LYS 252 1WDD.subA 5 0 2 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 9 2 465
LYS 305 1WDD.subA 8 1 1 2 0 2 0 2 2 4 4 0 2 2 0 0 2 0 8 2 465
LYS 316 1WDD.subA 13 1 0 1 0 2 0 2 2 10 3 0 1 1 2 2 5 0 1 2 465
LYS 334 1WDD.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 2 2 465
LYS 356 1WDD.subA 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 6 2 465
LYS 450 1WDD.subA 5 0 0 0 0 1 0 1 1 4 1 1 0 1 2 1 0 0 1 2 465
LYS 463 1WDD.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 2 2 0 0 2 2 465
LYS 466 1WDD.subA 6 0 1 1 0 0 1 0 1 4 2 0 0 0 1 1 2 0 5 2 465
LYS 474 1WDD.subA 5 1 1 2 0 0 0 0 0 3 2 0 0 0 1 0 2 0 1 2 465
LYS 8 1XGE.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 0 1 4 0 6 2 343
LYS 26 1XGE.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 13 2 343
LYS 131 1XGE.subA 5 0 3 3 0 1 0 1 1 1 4 0 0 0 0 0 1 0 4 2 343
LYS 172 1XGE.subA 8 0 0 0 0 1 0 1 1 7 1 0 0 0 3 2 2 0 7 2 343
LYS 181 1XGE.subA 6 2 0 2 0 0 1 0 1 3 3 1 0 1 1 0 1 0 7 2 343
LYS 226 1XGE.subA 11 0 2 2 0 1 1 1 2 7 4 2 0 2 1 1 3 0 6 2 343
LYS 259 1XGE.subA 12 0 1 1 0 2 2 2 4 7 5 0 0 0 6 1 0 0 3 2 343
LYS 346 1XGE.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 0 0 2 0 6 2 343
LYS 39 1XQK.subA 9 0 1 1 0 0 0 0 0 8 1 0 0 0 3 1 4 0 3 2 382
LYS 76 1XQK.subA 8 1 1 2 0 1 0 1 1 5 3 0 1 1 2 0 2 0 0 2 382
LYS 140 1XQK.subA 4 2 0 2 0 0 0 0 0 2 2 0 0 0 1 0 1 0 4 2 382
LYS 146 1XQK.subA 7 0 3 3 0 1 0 1 1 3 4 1 0 1 0 1 1 0 3 2 382
LYS 234 1XQK.subA 7 1 0 1 0 0 0 0 0 6 1 1 0 1 3 0 2 0 7 2 382
LYS 242 1XQK.subA 6 1 1 2 0 0 0 0 0 4 2 0 1 1 2 0 1 0 6 2 382
LYS 255 1XQK.subA 5 0 0 0 0 0 1 0 1 4 1 1 0 1 0 0 3 0 4 2 382
LYS 256 1XQK.subA 5 0 2 2 0 0 0 0 0 3 2 0 1 1 0 1 1 0 2 2 382
LYS 262 1XQK.subA 6 0 1 1 0 0 0 0 0 5 1 1 0 1 2 1 1 0 0 2 382
LYS 303 1XQK.subA 8 0 1 1 0 0 1 0 1 6 2 0 1 1 1 0 4 0 2 2 382
LYS 328 1XQK.subA 9 0 0 0 0 1 1 1 2 7 2 3 1 4 1 1 1 0 3 2 382
LYS 372 1XQK.subA 3 1 0 1 0 1 1 1 2 0 3 0 0 0 0 0 0 0 4 2 382
LYS 39 1XQL.subA 9 0 1 1 0 0 0 0 0 8 1 0 0 0 3 1 4 0 3 2 382
LYS 76 1XQL.subA 7 1 1 2 0 1 0 1 1 4 3 0 0 0 2 0 2 0 1 2 382
LYS 140 1XQL.subA 4 2 0 2 0 0 0 0 0 2 2 0 0 0 1 0 1 0 4 2 382
LYS 146 1XQL.subA 7 0 3 3 0 1 0 1 1 3 4 1 0 1 0 1 1 0 1 2 382
LYS 234 1XQL.subA 7 1 0 1 0 0 0 0 0 6 1 1 0 1 4 0 1 0 2 2 382
LYS 242 1XQL.subA 6 1 1 2 0 0 0 0 0 4 2 0 1 1 2 0 1 0 2 2 382
LYS 255 1XQL.subA 5 0 0 0 0 0 1 0 1 4 1 1 0 1 0 0 3 0 3 2 382
LYS 256 1XQL.subA 5 0 2 2 0 0 0 0 0 3 2 0 1 1 0 1 1 0 1 2 382
LYS 262 1XQL.subA 6 0 1 1 0 0 0 0 0 5 1 1 0 1 2 1 1 0 1 2 382
LYS 303 1XQL.subA 6 0 1 1 0 0 1 0 1 4 2 0 1 1 1 0 2 0 3 2 382
LYS 328 1XQL.subA 10 0 1 1 0 1 1 1 2 7 3 3 1 4 1 1 1 0 3 2 382
LYS 372 1XQL.subA 3 1 0 1 0 1 1 1 2 0 3 0 0 0 0 0 0 0 2 2 382
LYS 34 1YBQ.subB 8 1 0 1 0 0 0 0 0 7 1 0 2 2 3 0 2 0 0 2 389
LYS 119 1YBQ.subB 11 0 1 1 0 0 0 0 0 10 1 2 0 2 3 0 5 0 5 2 389
LYS 150 1YBQ.subB 5 1 1 2 0 0 0 0 0 3 2 2 0 2 0 0 1 0 1 2 389
LYS 194 1YBQ.subB 15 2 0 2 0 1 0 1 1 12 3 0 0 0 7 0 5 0 3 2 389
LYS 206 1YBQ.subB 1 0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 0 0 2 2 389
LYS 207 1YBQ.subB 5 0 0 0 0 0 0 0 0 5 0 2 1 3 2 0 0 0 2 2 389
LYS 225 1YBQ.subB 10 0 0 0 0 1 0 1 1 9 1 2 1 3 3 0 3 0 4 2 389
LYS 247 1YBQ.subB 7 0 1 1 0 1 0 1 1 5 2 0 0 0 2 2 1 0 0 2 389
LYS 319 1YBQ.subB 4 1 0 1 0 0 0 0 0 3 1 0 1 1 0 0 2 0 0 2 389
LYS 345 1YBQ.subB 12 2 0 2 0 1 0 1 1 9 3 3 0 3 3 1 2 0 3 2 389
LYS 373 1YBQ.subB 6 0 0 0 0 1 0 1 1 5 1 0 0 0 2 0 3 0 2 2 389
LYS 377 1YBQ.subB 5 1 0 1 0 0 0 0 0 4 1 0 1 1 1 0 2 0 1 2 389
LYS 380 1YBQ.subB 4 1 1 2 0 0 0 0 0 2 2 0 0 0 2 0 0 0 1 2 389
LYS 384 1YBQ.subB 6 0 1 1 0 0 0 0 0 5 1 2 0 2 2 0 1 0 0 2 389
LYS 34 2AQO.subA 8 1 0 1 0 0 0 0 0 7 1 0 2 2 3 0 2 0 2 2 378
LYS 119 2AQO.subA 11 0 1 1 0 0 0 0 0 10 1 2 0 2 3 0 5 0 4 2 378
LYS 150 2AQO.subA 5 1 1 2 0 0 0 0 0 3 2 2 0 2 0 0 1 0 1 2 378
LYS 194 2AQO.subA 15 2 0 2 0 1 0 1 1 12 3 0 0 0 7 0 5 0 3 2 378
LYS 206 2AQO.subA 1 0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 0 0 3 2 378
LYS 207 2AQO.subA 5 0 0 0 0 0 0 0 0 5 0 2 1 3 2 0 0 0 4 2 378
LYS 225 2AQO.subA 9 0 0 0 0 1 0 1 1 8 1 1 1 2 3 0 3 0 3 2 378
LYS 247 2AQO.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 1 2 1 0 1 2 378
LYS 319 2AQO.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 0 0 2 0 0 2 378
LYS 345 2AQO.subA 12 2 0 2 0 1 0 1 1 9 3 3 0 3 3 1 2 0 6 2 378
LYS 373 2AQO.subA 7 0 1 1 0 1 0 1 1 5 2 0 0 0 2 0 3 0 2 2 378
LYS 377 2AQO.subA 6 1 1 2 0 0 0 0 0 4 2 0 1 1 1 0 2 0 1 2 378
LYS 380 2AQO.subA 4 1 1 2 0 0 0 0 0 2 2 0 0 0 2 0 0 0 0 2 378
LYS 384 2AQO.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 2 2 378
LYS 34 2AQV.subA 8 1 0 1 0 0 0 0 0 7 1 0 2 2 3 0 2 0 1 2 375
LYS 119 2AQV.subA 11 0 1 1 0 0 0 0 0 10 1 2 0 2 3 0 5 0 5 2 375
LYS 150 2AQV.subA 5 1 1 2 0 0 0 0 0 3 2 2 0 2 0 0 1 0 2 2 375
LYS 194 2AQV.subA 15 2 0 2 0 1 0 1 1 12 3 0 0 0 7 0 5 0 2 2 375
LYS 206 2AQV.subA 1 0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 0 0 1 2 375
LYS 207 2AQV.subA 5 0 0 0 0 0 0 0 0 5 0 2 1 3 2 0 0 0 0 2 375
LYS 225 2AQV.subA 9 0 0 0 0 1 0 1 1 8 1 1 1 2 3 0 3 0 2 2 375
LYS 247 2AQV.subA 7 0 1 1 0 1 0 1 1 5 2 0 0 0 2 2 1 0 2 2 375
LYS 319 2AQV.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 0 0 2 0 0 2 375
LYS 345 2AQV.subA 12 2 0 2 0 1 0 1 1 9 3 3 0 3 3 1 2 0 4 2 375
LYS 373 2AQV.subA 7 0 1 1 0 1 0 1 1 5 2 0 0 0 2 0 3 0 2 2 375
LYS 377 2AQV.subA 6 1 1 2 0 0 0 0 0 4 2 0 1 1 1 0 2 0 1 2 375
LYS 380 2AQV.subA 4 1 1 2 0 0 0 0 0 2 2 0 0 0 2 0 0 0 0 2 375
LYS 384 2AQV.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 0 2 375
LYS 8 2E25.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 0 1 3 0 0 2 343
LYS 26 2E25.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 2 0 2 0 0 2 343
LYS 131 2E25.subA 4 0 2 2 0 1 0 1 1 1 3 0 0 0 0 0 1 0 0 2 343
LYS 172 2E25.subA 9 0 0 0 0 1 0 1 1 8 1 0 0 0 3 2 3 0 0 2 343
LYS 181 2E25.subA 6 2 0 2 0 0 1 0 1 3 3 1 0 1 1 0 1 0 0 2 343
LYS 226 2E25.subA 11 0 2 2 0 1 1 1 2 7 4 2 0 2 1 1 3 0 0 2 343
LYS 259 2E25.subA 12 0 1 1 0 2 2 2 4 7 5 0 0 0 6 1 0 0 1 2 343
LYS 346 2E25.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 0 0 2 0 0 2 343
LYS 8 2EG6.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 0 1 4 0 5 2 343
LYS 26 2EG6.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 12 2 343
LYS 131 2EG6.subA 4 0 2 2 0 1 0 1 1 1 3 0 0 0 0 0 1 0 4 2 343
LYS 172 2EG6.subA 8 0 0 0 0 1 0 1 1 7 1 0 0 0 3 1 3 0 8 2 343
LYS 181 2EG6.subA 6 2 0 2 0 0 1 0 1 3 3 1 0 1 1 0 1 0 6 2 343
LYS 226 2EG6.subA 11 0 2 2 0 1 1 1 2 7 4 2 0 2 1 1 3 0 6 2 343
LYS 259 2EG6.subA 12 0 1 1 0 2 2 2 4 7 5 0 0 0 6 1 0 0 3 2 343
LYS 346 2EG6.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 0 0 2 0 1 2 343
LYS 8 2EG7.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 0 1 4 0 3 2 343
LYS 26 2EG7.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 7 2 343
LYS 131 2EG7.subA 4 0 2 2 0 1 0 1 1 1 3 0 0 0 0 0 1 0 2 2 343
LYS 172 2EG7.subA 8 0 0 0 0 1 0 1 1 7 1 0 0 0 3 1 3 0 7 2 343
LYS 181 2EG7.subA 6 2 0 2 0 0 1 0 1 3 3 1 0 1 1 0 1 0 3 2 343
LYS 226 2EG7.subA 11 0 2 2 0 1 1 1 2 7 4 2 0 2 1 1 3 0 6 2 343
LYS 259 2EG7.subA 12 0 1 1 0 2 2 2 4 7 5 0 0 0 6 1 0 0 3 2 343
LYS 346 2EG7.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 0 0 3 0 1 2 343
LYS 8 2EG8.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 0 1 3 0 2 2 343
LYS 26 2EG8.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 6 2 343
LYS 131 2EG8.subA 4 0 2 2 0 1 0 1 1 1 3 0 0 0 0 0 1 0 1 2 343
LYS 172 2EG8.subA 8 0 0 0 0 1 0 1 1 7 1 0 0 0 3 2 2 0 5 2 343
LYS 181 2EG8.subA 6 2 0 2 0 0 1 0 1 3 3 1 0 1 1 0 1 0 2 2 343
LYS 226 2EG8.subA 11 0 2 2 0 1 1 1 2 7 4 2 0 2 1 1 3 0 6 2 343
LYS 259 2EG8.subA 12 0 1 1 0 2 2 2 4 7 5 0 0 0 6 1 0 0 3 2 343
LYS 346 2EG8.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 0 0 3 0 0 2 343
LYS 13 2FTW.subA 10 1 1 2 0 0 0 0 0 8 2 0 2 2 1 0 5 0 2 2 484
LYS 25 2FTW.subA 5 0 0 0 0 0 0 0 0 5 0 2 1 3 1 1 0 0 5 2 484
LYS 36 2FTW.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 0 0 4 0 2 2 484
LYS 40 2FTW.subA 3 0 0 0 0 0 0 0 0 3 0 2 1 3 0 0 0 0 3 2 484
LYS 45 2FTW.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 1 0 1 0 0 2 484
LYS 49 2FTW.subA 6 0 0 0 0 0 0 0 0 6 0 1 0 1 0 0 5 0 4 2 484
LYS 56 2FTW.subA 9 3 0 3 0 0 0 0 0 6 3 1 0 1 1 1 3 0 2 2 484
LYS 117 2FTW.subA 7 2 0 2 0 0 0 0 0 5 2 0 0 0 0 4 1 0 3 2 484
LYS 118 2FTW.subA 3 1 0 1 0 0 0 0 0 2 1 0 1 1 0 1 0 0 2 2 484
LYS 123 2FTW.subA 9 2 1 3 1 0 0 1 1 5 4 0 1 1 1 2 1 0 3 2 484
LYS 150 2FTW.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 0 2 0 0 2 484
LYS 164 2FTW.subA 6 1 0 1 0 0 0 0 0 5 1 0 1 1 1 2 1 0 3 2 484
LYS 179 2FTW.subA 5 0 0 0 0 1 1 1 2 3 2 0 0 0 0 2 1 0 3 2 484
LYS 182 2FTW.subA 9 1 1 2 0 0 0 0 0 7 2 1 0 1 2 2 2 0 3 2 484
LYS 203 2FTW.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 0 1 0 0 4 2 484
LYS 204 2FTW.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 2 2 484
LYS 259 2FTW.subA 6 1 1 2 0 0 1 0 1 3 3 0 1 1 1 0 1 0 3 2 484
LYS 262 2FTW.subA 4 0 1 1 0 1 1 1 2 1 3 0 0 0 1 0 0 0 4 2 484
LYS 307 2FTW.subA 7 0 0 0 0 0 0 0 0 7 0 1 1 2 3 0 2 0 8 2 484
LYS 334 2FTW.subA 14 0 0 0 0 1 0 1 1 13 1 0 2 2 7 2 2 0 5 2 484
LYS 338 2FTW.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 2 1 0 0 0 2 484
LYS 343 2FTW.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 1 0 2 0 0 2 484
LYS 367 2FTW.subA 9 1 0 1 0 0 0 0 0 8 1 1 1 2 3 0 3 0 4 2 484
LYS 391 2FTW.subA 11 3 0 3 0 2 0 2 2 6 5 0 0 0 2 1 3 0 6 2 484
LYS 411 2FTW.subA 7 1 1 2 0 0 0 0 0 5 2 2 0 2 0 0 3 0 6 2 484
LYS 415 2FTW.subA 7 2 1 3 0 0 1 0 1 3 4 1 1 2 0 1 0 0 2 2 484
LYS 431 2FTW.subA 6 0 0 0 0 0 0 0 0 6 0 3 0 3 1 0 2 0 3 2 484
LYS 451 2FTW.subA 4 1 0 1 0 0 0 0 0 3 1 1 1 2 0 0 1 0 4 2 484
LYS 456 2FTW.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 0 5 2 484
LYS 479 2FTW.subA 4 1 1 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 1 2 484
LYS 487 2FTW.subA 5 1 1 2 0 1 0 1 1 2 3 0 0 0 0 0 2 0 1 2 484
LYS 9 2FTY.subA 11 2 1 3 0 0 0 0 0 8 3 1 1 2 1 0 5 0 1 2 532
LYS 30 2FTY.subA 10 1 0 1 0 0 0 0 0 9 1 1 3 4 3 0 2 0 1 2 532
LYS 70 2FTY.subA 7 2 1 3 0 0 0 0 0 4 3 0 0 0 2 0 2 0 1 2 532
LYS 104 2FTY.subA 5 1 1 2 0 0 0 0 0 3 2 2 0 2 1 0 0 0 3 2 532
LYS 105 2FTY.subA 3 0 0 0 0 0 0 0 0 3 0 1 0 1 1 0 1 0 4 2 532
LYS 115 2FTY.subA 6 1 1 2 0 0 0 0 0 4 2 0 0 0 1 1 2 0 0 2 532
LYS 141 2FTY.subA 5 1 1 2 0 0 0 0 0 3 2 0 0 0 1 1 1 0 6 2 532
LYS 191 2FTY.subA 7 1 0 1 0 1 0 1 1 5 2 0 2 2 1 2 0 0 2 2 532
LYS 207 2FTY.subA 5 1 1 2 0 0 0 0 0 3 2 0 0 0 1 1 1 0 3 2 532
LYS 266 2FTY.subA 9 0 1 1 0 0 0 0 0 8 1 2 3 5 2 0 1 0 2 2 532
LYS 271 2FTY.subA 5 0 0 0 0 0 0 0 0 5 0 1 1 2 2 0 1 0 7 2 532
LYS 274 2FTY.subA 9 0 0 0 1 0 0 1 1 8 1 1 0 1 3 1 3 0 2 2 532
LYS 327 2FTY.subA 9 0 1 1 0 0 0 0 0 8 1 2 0 2 1 4 1 0 1 2 532
LYS 341 2FTY.subA 5 0 0 0 0 0 0 0 0 5 0 2 1 3 1 0 1 0 2 2 532
LYS 345 2FTY.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 1 1 0 0 2 2 532
LYS 367 2FTY.subA 5 0 1 1 0 0 1 0 1 3 2 2 0 2 0 1 0 0 0 2 532
LYS 373 2FTY.subA 13 0 0 0 0 0 1 0 1 12 1 2 2 4 2 4 2 0 1 2 532
LYS 383 2FTY.subA 10 1 1 2 0 1 0 1 1 7 3 0 3 3 2 2 0 0 0 2 532
LYS 419 2FTY.subA 10 0 1 1 0 0 0 0 0 9 1 2 1 3 1 0 5 0 4 2 532
LYS 430 2FTY.subA 7 0 0 0 0 0 0 0 0 7 0 1 1 2 1 2 2 0 3 2 532
LYS 438 2FTY.subA 12 3 0 3 1 0 0 1 1 8 4 2 1 3 2 0 3 0 4 2 532
LYS 459 2FTY.subA 2 1 0 1 0 0 0 0 0 1 1 1 0 1 0 0 0 0 1 2 532
LYS 460 2FTY.subA 5 0 2 2 0 1 0 1 1 2 3 1 0 1 0 1 0 0 1 2 532
LYS 465 2FTY.subA 7 0 0 0 0 0 0 0 0 7 0 1 1 2 2 3 0 0 0 2 532
LYS 467 2FTY.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 0 2 532
LYS 472 2FTY.subA 3 0 0 0 0 0 0 0 0 3 0 1 1 2 0 0 1 0 2 2 532
LYS 489 2FTY.subA 5 1 0 1 0 1 0 1 1 3 2 0 1 1 0 0 2 0 0 2 532
LYS 498 2FTY.subA 10 4 0 4 2 0 0 2 2 4 6 0 0 0 3 0 1 0 4 2 532
LYS 500 2FTY.subA 8 1 1 2 0 0 0 0 0 6 2 0 1 1 2 0 3 0 3 2 532
LYS 504 2FTY.subA 6 0 2 2 0 1 0 1 1 3 3 0 0 0 0 2 1 0 2 2 532
LYS 510 2FTY.subA 7 0 2 2 0 0 0 0 0 5 2 1 0 1 1 1 2 0 2 2 532
LYS 516 2FTY.subA 5 2 0 2 1 0 0 1 1 2 3 0 0 0 1 1 0 0 1 2 532
LYS 519 2FTY.subA 7 1 0 1 0 1 0 1 1 5 2 0 1 1 1 2 1 0 0 2 532
LYS 522 2FTY.subA 7 0 0 0 0 0 0 0 0 7 0 2 1 3 2 2 0 0 3 2 532
LYS 529 2FTY.subA 7 1 1 2 1 0 0 1 1 4 3 0 1 1 2 0 1 0 2 2 532
LYS 539 2FTY.subA 4 0 0 0 0 1 0 1 1 3 1 0 0 0 1 2 0 0 1 2 532
LYS 9 2FVK.subA 11 2 1 3 0 0 0 0 0 8 3 1 1 2 1 0 5 0 4 2 532
LYS 30 2FVK.subA 10 1 0 1 0 0 0 0 0 9 1 1 3 4 3 0 2 0 2 2 532
LYS 70 2FVK.subA 7 2 1 3 0 0 0 0 0 4 3 0 0 0 2 0 2 0 3 2 532
LYS 104 2FVK.subA 4 1 1 2 0 0 0 0 0 2 2 1 0 1 1 0 0 0 3 2 532
LYS 105 2FVK.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 1 0 1 0 1 2 532
LYS 115 2FVK.subA 6 1 1 2 0 0 0 0 0 4 2 0 0 0 1 1 2 0 1 2 532
LYS 141 2FVK.subA 6 1 1 2 0 0 0 0 0 4 2 1 0 1 1 1 1 0 6 2 532
LYS 191 2FVK.subA 7 1 0 1 0 1 0 1 1 5 2 0 2 2 1 2 0 0 2 2 532
LYS 207 2FVK.subA 5 1 1 2 0 0 0 0 0 3 2 0 0 0 1 1 1 0 2 2 532
LYS 266 2FVK.subA 9 0 1 1 0 0 0 0 0 8 1 2 3 5 2 0 1 0 4 2 532
LYS 271 2FVK.subA 5 0 0 0 0 0 0 0 0 5 0 1 1 2 2 0 1 0 10 2 532
LYS 274 2FVK.subA 11 0 0 0 1 0 0 1 1 10 1 1 0 1 4 2 3 0 4 2 532
LYS 327 2FVK.subA 9 0 1 1 0 0 0 0 0 8 1 2 0 2 1 4 1 0 2 2 532
LYS 341 2FVK.subA 5 0 0 0 1 0 0 1 1 4 1 1 1 2 1 0 1 0 2 2 532
LYS 345 2FVK.subA 5 0 0 0 1 0 0 1 1 4 1 1 1 2 1 1 0 0 4 2 532
LYS 367 2FVK.subA 6 0 1 1 0 0 1 0 1 4 2 2 0 2 0 2 0 0 2 2 532
LYS 373 2FVK.subA 13 0 0 0 0 0 1 0 1 12 1 2 2 4 2 4 2 0 1 2 532
LYS 383 2FVK.subA 8 1 0 1 0 1 0 1 1 6 2 0 3 3 1 2 0 0 1 2 532
LYS 419 2FVK.subA 10 0 1 1 0 0 0 0 0 9 1 2 1 3 1 0 5 0 5 2 532
LYS 430 2FVK.subA 7 0 0 0 0 0 0 0 0 7 0 1 1 2 1 2 2 0 2 2 532
LYS 438 2FVK.subA 12 3 0 3 1 0 0 1 1 8 4 2 1 3 2 0 3 0 5 2 532
LYS 459 2FVK.subA 2 1 0 1 0 0 0 0 0 1 1 1 0 1 0 0 0 0 0 2 532
LYS 460 2FVK.subA 5 0 2 2 0 1 0 1 1 2 3 1 0 1 0 1 0 0 2 2 532
LYS 465 2FVK.subA 7 0 0 0 0 0 0 0 0 7 0 1 1 2 2 3 0 0 4 2 532
LYS 467 2FVK.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 2 2 532
LYS 472 2FVK.subA 3 0 0 0 0 0 0 0 0 3 0 1 1 2 0 0 1 0 3 2 532
LYS 489 2FVK.subA 5 1 0 1 0 1 0 1 1 3 2 0 1 1 0 0 2 0 3 2 532
LYS 498 2FVK.subA 10 4 0 4 2 0 0 2 2 4 6 0 0 0 3 0 1 0 4 2 532
LYS 500 2FVK.subA 6 1 0 1 0 0 0 0 0 5 1 0 1 1 1 0 3 0 3 2 532
LYS 504 2FVK.subA 7 0 2 2 0 1 0 1 1 4 3 0 0 0 0 2 2 0 2 2 532
LYS 510 2FVK.subA 7 0 2 2 0 0 0 0 0 5 2 1 0 1 1 1 2 0 2 2 532
LYS 516 2FVK.subA 5 2 0 2 1 0 0 1 1 2 3 0 0 0 1 1 0 0 2 2 532
LYS 519 2FVK.subA 7 1 0 1 0 1 0 1 1 5 2 0 1 1 1 2 1 0 3 2 532
LYS 522 2FVK.subA 7 0 0 0 0 0 0 0 0 7 0 2 1 3 2 2 0 0 10 2 532
LYS 529 2FVK.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 2 0 0 0 5 2 532
LYS 539 2FVK.subA 4 0 0 0 0 1 0 1 1 3 1 0 0 0 1 2 0 0 1 2 532
LYS 9 2FVM.subA 11 2 1 3 0 0 0 0 0 8 3 1 1 2 1 0 5 0 0 2 532
LYS 30 2FVM.subA 10 1 0 1 0 0 0 0 0 9 1 1 3 4 3 0 2 0 3 2 532
LYS 70 2FVM.subA 7 2 1 3 0 0 0 0 0 4 3 0 0 0 2 0 2 0 3 2 532
LYS 104 2FVM.subA 4 1 1 2 0 0 0 0 0 2 2 1 0 1 1 0 0 0 1 2 532
LYS 105 2FVM.subA 3 0 0 0 0 0 0 0 0 3 0 1 0 1 1 0 1 0 1 2 532
LYS 115 2FVM.subA 6 1 1 2 0 0 0 0 0 4 2 0 0 0 1 1 2 0 0 2 532
LYS 141 2FVM.subA 5 1 1 2 0 0 0 0 0 3 2 0 0 0 1 1 1 0 4 2 532
LYS 191 2FVM.subA 7 1 0 1 0 1 0 1 1 5 2 0 2 2 1 2 0 0 2 2 532
LYS 207 2FVM.subA 5 1 1 2 0 0 0 0 0 3 2 0 0 0 1 1 1 0 1 2 532
LYS 266 2FVM.subA 9 0 1 1 0 0 0 0 0 8 1 2 3 5 2 0 1 0 4 2 532
LYS 271 2FVM.subA 5 0 0 0 0 0 0 0 0 5 0 1 1 2 2 0 1 0 6 2 532
LYS 274 2FVM.subA 10 0 0 0 1 0 0 1 1 9 1 1 0 1 3 2 3 0 3 2 532
LYS 327 2FVM.subA 9 0 1 1 0 0 0 0 0 8 1 2 0 2 1 4 1 0 1 2 532
LYS 341 2FVM.subA 6 0 0 0 1 0 0 1 1 5 1 1 2 3 0 1 1 0 4 2 532
LYS 345 2FVM.subA 5 0 0 0 1 0 0 1 1 4 1 1 1 2 1 1 0 0 2 2 532
LYS 367 2FVM.subA 5 0 1 1 0 0 1 0 1 3 2 2 0 2 0 1 0 0 1 2 532
LYS 373 2FVM.subA 13 0 0 0 0 0 1 0 1 12 1 2 2 4 2 4 2 0 1 2 532
LYS 383 2FVM.subA 8 1 0 1 0 1 0 1 1 6 2 0 3 3 1 2 0 0 0 2 532
LYS 419 2FVM.subA 10 0 1 1 0 0 0 0 0 9 1 2 1 3 1 0 5 0 4 2 532
LYS 430 2FVM.subA 7 0 0 0 0 0 0 0 0 7 0 1 1 2 1 2 2 0 3 2 532
LYS 438 2FVM.subA 12 3 0 3 1 0 0 1 1 8 4 2 1 3 2 0 3 0 5 2 532
LYS 459 2FVM.subA 1 0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 0 0 0 2 532
LYS 460 2FVM.subA 5 0 2 2 0 1 0 1 1 2 3 1 0 1 0 1 0 0 1 2 532
LYS 465 2FVM.subA 6 0 0 0 0 0 0 0 0 6 0 1 1 2 2 2 0 0 2 2 532
LYS 467 2FVM.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 0 2 532
LYS 472 2FVM.subA 3 0 0 0 0 0 0 0 0 3 0 1 1 2 0 0 1 0 2 2 532
LYS 489 2FVM.subA 5 1 0 1 0 1 0 1 1 3 2 0 1 1 0 0 2 0 1 2 532
LYS 498 2FVM.subA 11 4 0 4 2 0 0 2 2 5 6 0 0 0 4 0 1 0 5 2 532
LYS 500 2FVM.subA 5 0 0 0 0 0 0 0 0 5 0 0 1 1 1 0 3 0 4 2 532
LYS 504 2FVM.subA 7 0 2 2 0 1 0 1 1 4 3 0 0 0 0 2 2 0 2 2 532
LYS 510 2FVM.subA 6 0 2 2 0 0 0 0 0 4 2 1 0 1 0 1 2 0 1 2 532
LYS 516 2FVM.subA 5 2 0 2 1 0 0 1 1 2 3 0 0 0 1 1 0 0 4 2 532
LYS 519 2FVM.subA 7 1 0 1 0 1 0 1 1 5 2 0 1 1 1 2 1 0 0 2 532
LYS 522 2FVM.subA 7 0 0 0 0 0 0 0 0 7 0 2 1 3 2 2 0 0 6 2 532
LYS 529 2FVM.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 2 0 0 0 3 2 532
LYS 539 2FVM.subA 4 0 0 0 0 1 0 1 1 3 1 0 0 0 1 2 0 0 1 2 532
LYS 18 2GC5.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 0 1 0 1 2 408
LYS 50 2GC5.subA 13 0 1 1 0 0 0 0 0 12 1 4 1 5 3 0 3 0 5 2 408
LYS 172 2GC5.subA 2 0 0 0 0 0 0 0 0 2 0 0 1 1 0 0 1 0 0 2 408
LYS 183 2GC5.subA 5 0 0 0 0 0 1 0 1 4 1 1 0 1 2 0 1 0 7 2 408
LYS 190 2GC5.subA 10 0 0 0 0 1 1 1 2 8 2 1 1 2 2 0 4 0 4 2 408
LYS 211 2GC5.subA 9 0 0 0 0 0 0 0 0 9 0 2 0 2 3 0 4 0 5 2 408
LYS 230 2GC5.subA 5 1 0 1 0 1 0 1 1 3 2 1 0 1 2 0 0 0 6 2 408
LYS 232 2GC5.subA 3 0 0 0 0 1 0 1 1 2 1 0 0 0 1 0 1 0 2 2 408
LYS 273 2GC5.subA 10 2 0 2 0 0 0 0 0 8 2 0 1 1 1 1 5 0 6 2 408
LYS 277 2GC5.subA 3 0 0 0 0 0 0 0 0 3 0 0 1 1 1 0 1 0 3 2 408
LYS 303 2GC5.subA 8 1 1 2 0 0 0 0 0 6 2 1 0 1 0 1 4 0 0 2 408
LYS 329 2GC5.subA 5 0 0 0 0 0 0 0 0 5 0 1 1 2 0 1 2 0 0 2 408
LYS 388 2GC5.subA 10 1 1 2 0 1 0 1 1 7 3 1 0 1 2 1 3 0 1 2 408
LYS 50 2GC6.subA 11 0 1 1 0 0 1 0 1 9 2 3 0 3 4 0 2 0 3 2 407
LYS 172 2GC6.subA 3 0 0 0 0 0 1 0 1 2 1 0 1 1 0 0 1 0 1 2 407
LYS 183 2GC6.subA 5 0 0 0 0 0 1 0 1 4 1 1 0 1 2 0 1 0 7 2 407
LYS 190 2GC6.subA 9 0 0 0 0 1 0 1 1 8 1 1 1 2 2 0 4 0 7 2 407
LYS 211 2GC6.subA 8 0 0 0 0 0 0 0 0 8 0 1 0 1 3 0 4 0 5 2 407
LYS 230 2GC6.subA 4 0 0 0 0 1 0 1 1 3 1 1 0 1 2 0 0 0 7 2 407
LYS 232 2GC6.subA 3 0 0 0 0 1 0 1 1 2 1 0 0 0 1 0 1 0 4 2 407
LYS 273 2GC6.subA 9 2 0 2 0 0 0 0 0 7 2 0 1 1 1 1 4 0 6 2 407
LYS 277 2GC6.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 1 0 1 0 5 2 407
LYS 303 2GC6.subA 8 1 1 2 0 0 0 0 0 6 2 1 0 1 0 1 4 0 1 2 407
LYS 329 2GC6.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 0 0 2 0 0 2 407
LYS 388 2GC6.subA 10 1 1 2 0 1 0 1 1 7 3 1 0 1 2 1 3 0 1 2 407
LYS 2 2GWN.subA 9 2 0 2 0 0 0 0 0 5 2 1 0 1 1 0 3 0 3 2 451
LYS 15 2GWN.subA 5 0 1 1 0 0 0 0 0 4 1 1 1 2 1 1 0 0 4 2 451
LYS 74 2GWN.subA 8 0 1 1 0 1 2 1 3 4 4 0 0 0 2 0 2 0 1 2 451
LYS 111 2GWN.subA 12 1 2 3 0 2 0 2 2 6 5 2 0 2 2 1 1 0 5 2 451
LYS 138 2GWN.subA 5 1 1 2 0 1 0 1 1 2 3 0 0 0 0 0 2 0 1 2 451
LYS 142 2GWN.subA 6 2 0 2 0 0 1 0 1 3 3 0 1 1 0 0 2 0 7 2 451
LYS 164 2GWN.subA 4 1 2 3 0 0 0 0 0 1 3 0 1 1 0 0 0 0 3 2 451
LYS 169 2GWN.subA 6 0 3 3 0 0 0 0 0 3 3 1 0 1 0 0 2 0 1 2 451
LYS 184 2GWN.subA 8 0 4 4 0 1 0 1 1 3 5 1 0 1 1 0 1 0 4 2 451
LYS 192 2GWN.subA 7 0 1 1 0 2 0 2 2 4 3 0 1 1 0 0 3 0 4 2 451
LYS 196 2GWN.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 2 1 2 0 2 2 451
LYS 198 2GWN.subA 4 0 0 0 0 0 1 0 1 3 1 0 0 0 1 2 0 0 3 2 451
LYS 244 2GWN.subA 6 0 4 4 0 0 0 0 0 2 4 1 0 1 0 1 0 0 3 2 451
LYS 258 2GWN.subA 12 1 0 1 0 3 0 3 3 8 4 2 2 4 2 0 2 0 6 2 451
LYS 284 2GWN.subA 14 0 0 0 1 0 1 1 2 12 2 2 2 4 4 1 3 0 7 2 451
LYS 290 2GWN.subA 12 1 1 2 1 1 0 2 2 8 4 3 0 3 4 0 1 0 7 2 451
LYS 291 2GWN.subA 6 2 2 4 1 0 0 1 1 1 5 1 0 1 0 0 0 0 10 2 451
LYS 322 2GWN.subA 13 0 2 2 1 1 1 2 3 8 5 1 1 2 3 1 2 0 3 2 451
LYS 360 2GWN.subA 11 1 1 2 0 0 0 0 0 9 2 2 0 2 1 0 6 0 4 2 451
LYS 373 2GWN.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 1 3 0 0 6 2 451
LYS 432 2GWN.subA 4 0 0 0 0 0 0 0 0 4 0 0 0 0 3 1 0 0 7 2 451
LYS 10 2ICS.subA 9 0 1 1 0 0 1 0 1 7 2 2 1 3 0 0 4 0 0 2 368
LYS 25 2ICS.subA 6 2 1 3 0 0 0 0 0 3 3 0 0 0 1 1 1 0 0 2 368
LYS 27 2ICS.subA 4 0 1 1 1 0 0 1 1 2 2 0 0 0 2 0 0 0 1 2 368
LYS 28 2ICS.subA 5 1 1 2 1 0 0 1 1 2 3 0 0 0 1 0 1 0 2 2 368
LYS 41 2ICS.subA 4 1 1 2 0 0 0 0 0 2 2 0 0 0 1 1 0 0 0 2 368
LYS 66 2ICS.subA 5 0 2 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 1 2 368
LYS 80 2ICS.subA 5 1 1 2 1 0 0 1 1 2 3 0 0 0 1 0 1 0 1 2 368
LYS 81 2ICS.subA 9 1 1 2 1 1 0 2 2 5 4 0 0 0 2 1 2 0 3 2 368
LYS 108 2ICS.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 1 1 2 0 2 2 368
LYS 119 2ICS.subA 8 0 1 1 1 0 0 1 1 6 2 2 1 3 1 1 1 0 3 2 368
LYS 133 2ICS.subA 9 2 1 3 1 0 0 1 1 5 4 2 0 2 1 0 2 0 2 2 368
LYS 140 2ICS.subA 6 0 1 1 0 0 0 0 0 5 1 1 2 3 1 0 1 0 0 2 368
LYS 141 2ICS.subA 4 0 1 1 0 0 0 0 0 3 1 1 0 1 1 0 1 0 1 2 368
LYS 174 2ICS.subA 10 1 2 3 0 0 0 0 0 6 3 0 2 2 1 0 3 0 0 2 368
LYS 205 2ICS.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 2 0 1 0 3 2 368
LYS 216 2ICS.subA 10 1 1 2 0 0 0 0 0 8 2 2 2 4 3 0 1 0 4 2 368
LYS 227 2ICS.subA 4 1 1 2 0 0 0 0 0 2 2 1 0 1 0 0 1 0 3 2 368
LYS 229 2ICS.subA 10 2 1 3 0 0 0 0 0 7 3 0 1 1 2 1 3 0 1 2 368
LYS 238 2ICS.subA 9 0 1 1 0 0 0 0 0 7 1 0 2 2 3 0 2 0 3 2 368
LYS 265 2ICS.subA 8 0 1 1 1 0 0 1 1 5 2 0 0 0 2 2 1 0 2 2 368
LYS 292 2ICS.subA 13 1 1 2 0 0 1 0 1 10 3 3 0 3 1 2 4 0 2 2 368
LYS 306 2ICS.subA 11 0 2 2 1 0 0 1 1 8 3 1 0 1 3 1 3 0 2 2 368
LYS 309 2ICS.subA 6 0 2 2 0 0 0 0 0 4 2 2 0 2 2 0 0 0 1 2 368
LYS 320 2ICS.subA 9 2 0 2 0 0 0 0 0 7 2 1 2 3 1 1 2 0 7 2 368
LYS 327 2ICS.subA 9 1 1 2 0 0 0 0 0 7 2 2 1 3 2 0 2 0 1 2 368
LYS 341 2ICS.subA 8 0 3 3 0 0 0 0 0 5 3 2 1 3 1 0 1 0 1 2 368
LYS 354 2ICS.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 1 0 1 0 1 2 368
LYS 361 2ICS.subA 8 1 0 1 0 0 0 0 0 7 1 1 0 1 0 1 5 0 0 2 368
LYS 34 2J6V.subA 10 0 2 2 0 2 0 2 2 6 4 1 1 2 2 0 2 0 1 2 281
LYS 151 2J6V.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 1 1 0 4 2 281
LYS 153 2J6V.subA 5 0 1 1 0 2 0 2 2 2 3 0 0 0 2 0 0 0 3 2 281
LYS 188 2J6V.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 0 2 0 3 2 281
LYS 238 2J6V.subA 1 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 2 2 281
LYS 239 2J6V.subA 6 1 0 1 1 1 0 2 2 3 3 0 1 1 2 0 0 0 4 2 281
LYS 271 2J6V.subA 7 1 0 1 1 0 1 1 2 4 3 0 0 0 3 1 0 0 2 2 281
LYS 273 2J6V.subA 10 0 2 2 0 0 0 0 0 8 2 0 2 2 4 1 1 0 3 2 281
LYS 6 2JFF.subA 7 2 0 2 0 0 0 0 0 5 2 0 2 2 1 1 1 0 5 2 433
LYS 45 2JFF.subA 5 2 0 2 0 0 0 0 0 3 2 0 0 0 1 0 2 0 1 2 433
LYS 115 2JFF.subA 14 0 1 1 0 0 0 0 0 13 1 5 2 7 2 0 3 0 9 2 433
LYS 127 2JFF.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 3 0 2 0 8 2 433
LYS 206 2JFF.subA 5 0 0 0 0 0 0 0 0 5 0 0 1 1 2 0 2 0 5 2 433
LYS 251 2JFF.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 2 2 433
LYS 254 2JFF.subA 7 0 1 1 0 1 0 1 1 5 2 0 2 2 0 1 2 0 2 2 433
LYS 259 2JFF.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 6 2 433
LYS 262 2JFF.subA 6 0 0 0 0 0 0 0 0 6 0 3 0 3 0 1 2 0 5 2 433
LYS 291 2JFF.subA 7 0 0 0 0 0 0 0 0 7 0 3 0 3 3 0 1 0 5 2 433
LYS 319 2JFF.subA 11 1 1 2 0 1 1 1 2 7 4 2 0 2 2 1 2 0 7 2 433
LYS 348 2JFF.subA 12 3 0 3 0 0 1 0 1 8 4 2 1 3 4 0 1 0 9 2 433
LYS 420 2JFF.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 1 0 0 2 2 433
LYS 434 2JFF.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 3 2 433
LYS 6 2JFG.subA 7 2 0 2 0 0 0 0 0 5 2 0 2 2 1 1 1 0 3 2 440
LYS 45 2JFG.subA 5 1 0 1 0 0 0 0 0 4 1 0 0 0 2 0 2 0 5 2 440
LYS 115 2JFG.subA 14 0 1 1 0 0 0 0 0 13 1 5 2 7 2 0 3 0 4 2 440
LYS 127 2JFG.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 3 0 2 0 8 2 440
LYS 206 2JFG.subA 4 0 0 0 0 0 0 0 0 4 0 0 0 0 2 0 2 0 4 2 440
LYS 251 2JFG.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 3 2 440
LYS 254 2JFG.subA 7 0 1 1 0 1 0 1 1 5 2 0 2 2 0 1 2 0 3 2 440
LYS 259 2JFG.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 7 2 440
LYS 262 2JFG.subA 7 0 1 1 0 0 0 0 0 6 1 3 0 3 0 1 2 0 3 2 440
LYS 291 2JFG.subA 5 0 0 0 0 0 0 0 0 5 0 2 0 2 2 0 1 0 0 2 440
LYS 319 2JFG.subA 9 1 1 2 0 1 1 1 2 5 4 2 1 3 1 1 0 0 5 2 440
LYS 348 2JFG.subA 12 3 0 3 0 0 1 0 1 8 4 2 1 3 4 0 1 0 10 2 440
LYS 420 2JFG.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 1 0 0 8 2 440
LYS 434 2JFG.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 7 2 440
LYS 6 2JFH.subA 7 2 0 2 0 0 0 0 0 5 2 0 2 2 1 1 1 0 2 2 431
LYS 45 2JFH.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 2 2 431
LYS 115 2JFH.subA 14 0 1 1 0 0 0 0 0 13 1 5 2 7 2 0 3 0 4 2 431
LYS 127 2JFH.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 3 0 2 0 6 2 431
LYS 206 2JFH.subA 4 0 0 0 0 0 0 0 0 4 0 0 0 0 2 0 2 0 5 2 431
LYS 251 2JFH.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 3 2 431
LYS 254 2JFH.subA 6 0 1 1 0 1 0 1 1 4 2 0 1 1 0 1 2 0 2 2 431
LYS 259 2JFH.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 2 2 431
LYS 262 2JFH.subA 6 0 0 0 0 0 0 0 0 6 0 3 0 3 0 1 2 0 1 2 431
LYS 291 2JFH.subA 7 0 0 0 0 0 0 0 0 7 0 3 0 3 3 0 1 0 1 2 431
LYS 319 2JFH.subA 9 1 0 1 0 1 1 1 2 6 3 2 0 2 2 1 1 0 2 2 431
LYS 348 2JFH.subA 12 3 0 3 0 0 1 0 1 8 4 2 1 3 4 0 1 0 6 2 431
LYS 420 2JFH.subA 4 0 0 0 0 0 0 0 0 4 0 0 2 2 1 1 0 0 0 2 431
LYS 434 2JFH.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 1 2 431
LYS 20 2KAU.subC 7 2 1 3 0 0 0 0 0 4 3 0 0 0 1 1 2 0 1 2 566
LYS 44 2KAU.subC 8 0 1 1 1 0 0 1 1 6 2 0 1 1 0 2 3 0 4 2 566
LYS 49 2KAU.subC 7 0 0 0 1 0 0 1 1 6 1 0 0 0 3 1 2 0 2 2 566
LYS 83 2KAU.subC 14 2 0 2 0 0 0 0 0 12 2 1 1 2 4 2 4 0 4 2 566
LYS 89 2KAU.subC 10 3 0 3 0 1 0 1 1 6 4 0 0 0 3 1 2 0 1 2 566
LYS 98 2KAU.subC 15 1 2 3 0 0 0 0 0 12 3 3 2 5 4 0 3 0 4 2 566
LYS 124 2KAU.subC 10 0 1 1 0 0 0 0 0 9 1 1 0 1 5 1 2 0 1 2 566
LYS 196 2KAU.subC 7 0 0 0 0 0 0 0 0 7 0 0 2 2 4 1 0 0 0 2 566
LYS 382 2KAU.subC 10 1 1 2 0 1 1 1 2 6 4 0 0 0 2 1 3 0 0 2 566
LYS 401 2KAU.subC 11 1 1 2 0 1 2 1 3 6 5 0 1 1 1 2 2 0 2 2 566
LYS 406 2KAU.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 2 3 2 0 0 2 566
LYS 428 2KAU.subC 12 0 2 2 0 0 1 0 1 9 3 1 0 1 3 0 5 0 1 2 566
LYS 443 2KAU.subC 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 4 1 0 0 2 566
LYS 449 2KAU.subC 14 1 0 1 0 1 0 1 1 12 2 1 1 2 4 0 6 0 0 2 566
LYS 515 2KAU.subC 7 1 0 1 0 0 1 0 1 5 2 0 1 1 2 0 2 0 3 2 566
LYS 522 2KAU.subC 17 0 0 0 0 0 1 0 1 16 1 2 3 5 7 1 3 0 3 2 566
LYS 77 2O4M.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 10 2 330
LYS 82 2O4M.subA 9 1 1 2 0 1 0 1 1 6 3 0 0 0 2 1 3 0 12 2 330
LYS 175 2O4M.subA 4 0 0 0 0 0 0 0 0 4 0 0 2 2 2 0 0 0 8 2 330
LYS 185 2O4M.subA 9 0 2 2 0 1 0 1 1 6 3 1 0 1 1 0 4 0 6 2 330
LYS 285 2O4M.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 4 2 330
LYS 294 2O4M.subA 4 0 0 0 0 1 0 1 1 3 1 0 2 2 0 0 1 0 9 2 330
LYS 339 2O4M.subA 6 0 1 1 1 0 0 1 1 4 2 0 0 0 1 1 2 0 3 2 330
LYS 77 2O4Q.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 5 2 331
LYS 82 2O4Q.subA 9 1 1 2 0 1 0 1 1 6 3 0 0 0 3 1 2 0 11 2 331
LYS 175 2O4Q.subA 3 0 0 0 0 0 0 0 0 3 0 0 1 1 2 0 0 0 5 2 331
LYS 185 2O4Q.subA 8 0 2 2 0 1 0 1 1 5 3 0 0 0 1 0 4 0 0 2 331
LYS 285 2O4Q.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 4 2 331
LYS 294 2O4Q.subA 4 0 0 0 0 1 0 1 1 3 1 0 2 2 0 0 1 0 6 2 331
LYS 339 2O4Q.subA 7 0 1 1 1 0 0 1 1 5 2 0 0 0 1 1 3 0 1 2 331
LYS 77 2OB3.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 8 2 330
LYS 82 2OB3.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 2 1 2 0 12 2 330
LYS 175 2OB3.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 5 2 330
LYS 185 2OB3.subA 8 0 2 2 0 1 0 1 1 5 3 0 0 0 1 0 4 0 3 2 330
LYS 285 2OB3.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 4 2 330
LYS 294 2OB3.subA 4 0 0 0 0 1 0 1 1 3 1 0 2 2 0 0 1 0 7 2 330
LYS 339 2OB3.subA 6 0 1 1 1 0 0 1 1 4 2 0 0 0 1 1 2 0 3 2 330
LYS 27 2ODO.subA 6 2 1 3 0 0 0 0 0 3 3 0 0 0 2 1 0 0 1 2 356
LYS 115 2ODO.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 0 1 0 0 2 356
LYS 138 2ODO.subA 6 1 0 1 0 0 0 0 0 5 1 0 1 1 2 1 1 0 0 2 356
LYS 152 2ODO.subA 5 0 1 1 0 0 1 0 1 3 2 0 0 0 1 0 2 0 2 2 356
LYS 238 2ODO.subA 6 0 0 0 0 0 0 0 0 6 0 3 0 3 2 0 1 0 0 2 356
LYS 261 2ODO.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 0 2 356
LYS 325 2ODO.subA 6 1 0 1 0 0 0 0 0 5 1 1 1 2 3 0 0 0 1 2 356
LYS 348 2ODO.subA 5 1 0 1 0 1 0 1 1 3 2 0 1 1 1 0 1 0 0 2 356
LYS 18 2OEK.subA 3 1 0 1 0 1 0 1 1 1 2 1 0 1 0 0 0 0 1 2 412
LYS 19 2OEK.subA 10 3 1 4 0 0 0 0 0 6 4 0 0 0 2 0 4 0 6 2 412
LYS 45 2OEK.subA 5 0 0 0 0 1 1 1 2 3 2 0 1 1 0 0 2 0 5 2 412
LYS 47 2OEK.subA 9 0 1 1 0 1 1 1 2 6 3 1 0 1 1 1 3 0 7 2 412
LYS 68 2OEK.subA 6 0 2 2 0 1 0 1 1 3 3 0 0 0 1 1 1 0 2 2 412
LYS 71 2OEK.subA 8 2 2 4 0 1 0 1 1 3 5 1 0 1 0 0 2 0 5 2 412
LYS 76 2OEK.subA 8 0 0 0 0 0 0 0 0 8 0 1 0 1 2 0 5 0 3 2 412
LYS 98 2OEK.subA 4 0 0 0 0 0 0 0 0 4 0 0 0 0 1 0 3 0 6 2 412
LYS 117 2OEK.subA 9 1 1 2 0 1 0 1 1 6 3 0 1 1 2 3 0 0 4 2 412
LYS 147 2OEK.subA 5 0 0 0 0 0 1 0 1 4 1 0 0 0 2 1 1 0 3 2 412
LYS 162 2OEK.subA 7 0 1 1 0 0 0 0 0 6 1 2 0 2 1 0 3 0 3 2 412
LYS 163 2OEK.subA 5 0 1 1 0 0 0 0 0 4 1 1 1 2 1 0 1 0 2 2 412
LYS 188 2OEK.subA 6 0 2 2 0 1 0 1 1 3 3 1 0 1 1 0 1 0 3 2 412
LYS 194 2OEK.subA 6 0 0 0 0 0 0 0 0 6 0 1 0 1 3 0 2 0 5 2 412
LYS 206 2OEK.subA 12 1 0 1 0 3 0 3 3 8 4 2 0 2 3 0 3 0 3 2 412
LYS 217 2OEK.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 2 1 0 0 4 2 412
LYS 222 2OEK.subA 6 2 0 2 0 0 0 0 0 4 2 0 0 0 1 0 3 0 6 2 412
LYS 224 2OEK.subA 9 2 0 2 0 0 0 0 0 7 2 1 0 1 3 1 2 0 10 2 412
LYS 226 2OEK.subA 6 1 2 3 0 1 0 1 1 2 4 0 0 0 1 0 1 0 3 2 412
LYS 287 2OEK.subA 10 1 2 3 0 0 0 0 0 7 3 1 0 1 2 1 3 0 10 2 412
LYS 401 2OEK.subA 7 1 1 2 0 0 1 0 1 4 3 0 1 1 3 0 0 0 0 2 412
LYS 18 2OEL.subA 3 1 1 2 0 1 0 1 1 0 3 0 0 0 0 0 0 0 1 2 412
LYS 19 2OEL.subA 10 3 1 4 0 0 0 0 0 6 4 0 0 0 2 0 4 0 6 2 412
LYS 45 2OEL.subA 5 0 0 0 0 1 1 1 2 3 2 0 1 1 0 0 2 0 4 2 412
LYS 47 2OEL.subA 9 0 1 1 0 1 1 1 2 6 3 1 0 1 1 1 3 0 5 2 412
LYS 68 2OEL.subA 5 0 1 1 0 1 0 1 1 3 2 0 0 0 1 1 1 0 0 2 412
LYS 71 2OEL.subA 8 2 2 4 0 1 0 1 1 3 5 1 0 1 0 0 2 0 3 2 412
LYS 76 2OEL.subA 8 0 0 0 0 0 0 0 0 8 0 1 0 1 2 0 5 0 4 2 412
LYS 98 2OEL.subA 4 0 0 0 0 0 0 0 0 4 0 0 0 0 1 0 3 0 7 2 412
LYS 117 2OEL.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 2 3 0 0 4 2 412
LYS 147 2OEL.subA 5 0 0 0 0 0 1 0 1 4 1 0 0 0 2 1 1 0 2 2 412
LYS 162 2OEL.subA 7 0 1 1 0 0 0 0 0 6 1 2 0 2 1 0 3 0 3 2 412
LYS 163 2OEL.subA 6 0 1 1 0 0 0 0 0 5 1 1 1 2 2 0 1 0 2 2 412
LYS 188 2OEL.subA 6 0 2 2 0 1 0 1 1 3 3 1 0 1 1 0 1 0 2 2 412
LYS 194 2OEL.subA 7 0 0 0 0 0 0 0 0 7 0 1 1 2 3 0 2 0 4 2 412
LYS 206 2OEL.subA 12 1 0 1 0 3 0 3 3 8 4 2 0 2 3 0 3 0 1 2 412
LYS 217 2OEL.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 2 1 0 0 5 2 412
LYS 222 2OEL.subA 6 2 0 2 0 0 0 0 0 4 2 0 0 0 1 0 3 0 7 2 412
LYS 224 2OEL.subA 9 2 0 2 0 0 0 0 0 7 2 1 0 1 3 1 2 0 11 2 412
LYS 226 2OEL.subA 6 1 2 3 0 1 0 1 1 2 4 0 0 0 1 0 1 0 2 2 412
LYS 287 2OEL.subA 10 1 2 3 0 0 0 0 0 7 3 1 0 1 2 1 3 0 10 2 412
LYS 401 2OEL.subA 7 1 1 2 0 0 1 0 1 4 3 0 1 1 3 0 0 0 1 2 412
LYS 18 2OEM.subA 2 1 0 1 0 1 0 1 1 0 2 0 0 0 0 0 0 0 2 2 408
LYS 19 2OEM.subA 10 3 1 4 0 0 0 0 0 6 4 0 0 0 2 0 4 0 6 2 408
LYS 45 2OEM.subA 5 0 0 0 0 1 1 1 2 3 2 0 1 1 0 0 2 0 2 2 408
LYS 47 2OEM.subA 9 0 1 1 0 1 1 1 2 6 3 1 0 1 1 1 3 0 5 2 408
LYS 68 2OEM.subA 5 0 1 1 0 1 0 1 1 3 2 0 0 0 1 1 1 0 0 2 408
LYS 71 2OEM.subA 8 2 2 4 0 1 0 1 1 3 5 1 0 1 0 0 2 0 1 2 408
LYS 76 2OEM.subA 8 0 0 0 0 0 0 0 0 8 0 1 0 1 2 0 5 0 3 2 408
LYS 98 2OEM.subA 4 0 0 0 0 0 0 0 0 4 0 0 0 0 1 0 3 0 4 2 408
LYS 117 2OEM.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 2 3 0 0 4 2 408
LYS 147 2OEM.subA 8 1 0 1 0 0 1 0 1 6 2 0 0 0 3 1 1 1 4 2 408
LYS 162 2OEM.subA 8 0 1 1 0 0 0 0 0 7 1 2 1 3 1 0 3 0 2 2 408
LYS 163 2OEM.subA 5 0 1 1 0 0 0 0 0 4 1 1 1 2 1 0 1 0 3 2 408
LYS 188 2OEM.subA 6 0 2 2 0 1 0 1 1 3 3 1 0 1 1 0 1 0 5 2 408
LYS 194 2OEM.subA 8 0 0 0 0 1 0 1 1 7 1 1 1 2 3 0 2 0 4 2 408
LYS 206 2OEM.subA 11 1 0 1 0 3 0 3 3 7 4 2 0 2 2 0 3 0 6 2 408
LYS 217 2OEM.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 2 1 0 0 5 2 408
LYS 222 2OEM.subA 6 2 0 2 0 0 0 0 0 4 2 0 0 0 1 0 3 0 4 2 408
LYS 224 2OEM.subA 9 2 0 2 0 0 0 0 0 7 2 1 0 1 3 1 2 0 10 2 408
LYS 226 2OEM.subA 6 1 2 3 0 1 0 1 1 2 4 0 0 0 1 0 1 0 3 2 408
LYS 287 2OEM.subA 10 1 2 3 0 0 0 0 0 7 3 1 0 1 2 1 3 0 10 2 408
LYS 401 2OEM.subA 7 1 1 2 0 0 1 0 1 4 3 0 1 1 3 0 0 0 2 2 408
LYS 23 2OGJ.subA 7 0 0 0 0 0 0 0 0 7 0 2 1 3 2 1 1 0 0 2 379
LYS 29 2OGJ.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 2 0 0 0 0 2 379
LYS 45 2OGJ.subA 6 1 0 1 0 1 0 1 1 4 2 0 1 1 2 0 1 0 1 2 379
LYS 130 2OGJ.subA 10 0 1 1 0 1 1 1 2 7 3 1 0 1 1 0 5 0 2 2 379
LYS 153 2OGJ.subA 3 2 0 2 0 0 0 0 0 1 2 0 0 0 0 0 1 0 0 2 379
LYS 192 2OGJ.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 0 0 4 0 0 2 379
LYS 195 2OGJ.subA 8 1 0 1 1 0 0 1 1 5 2 0 0 0 2 0 3 0 0 2 379
LYS 196 2OGJ.subA 4 0 0 0 1 0 0 1 1 3 1 0 0 0 0 0 3 0 0 2 379
LYS 199 2OGJ.subA 4 0 0 0 2 0 0 2 2 2 2 0 0 0 1 0 1 0 0 2 379
LYS 202 2OGJ.subA 4 0 0 0 0 0 0 0 0 4 0 0 0 0 1 0 3 0 0 2 379
LYS 236 2OGJ.subA 7 0 0 0 0 0 0 0 0 7 0 3 2 5 2 0 0 0 1 2 379
LYS 273 2OGJ.subA 4 1 0 1 0 0 0 0 0 3 1 1 0 1 0 1 1 0 0 2 379
LYS 310 2OGJ.subA 14 1 0 1 0 0 1 0 1 12 2 5 0 5 2 2 3 0 0 2 379
LYS 372 2OGJ.subA 5 1 0 1 0 2 0 2 2 2 3 0 0 0 1 0 1 0 0 2 379
LYS 77 2OQL.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 7 2 330
LYS 82 2OQL.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 2 1 2 0 11 2 330
LYS 175 2OQL.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 3 2 330
LYS 185 2OQL.subA 8 0 2 2 0 1 0 1 1 5 3 0 0 0 1 0 4 0 4 2 330
LYS 285 2OQL.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 4 2 330
LYS 294 2OQL.subA 4 0 0 0 0 1 0 1 1 3 1 0 2 2 0 0 1 0 6 2 330
LYS 339 2OQL.subA 6 0 1 1 1 0 0 1 1 4 2 0 0 0 1 1 2 0 1 2 330
LYS 24 2P9V.subA 6 0 1 1 1 0 0 1 1 4 2 0 2 2 0 0 2 0 2 2 355
LYS 37 2P9V.subA 5 0 0 0 0 0 0 0 0 5 0 0 1 1 2 2 0 0 3 2 355
LYS 50 2P9V.subA 3 1 0 1 0 0 0 0 0 2 1 0 1 1 1 0 0 0 7 2 355
LYS 51 2P9V.subA 5 1 0 1 0 0 0 0 0 4 1 0 1 1 1 1 1 0 6 2 355
LYS 67 2P9V.subA 14 0 1 1 0 0 0 0 0 13 1 3 1 4 3 3 2 0 1 2 355
LYS 84 2P9V.subA 5 1 0 1 1 0 0 1 1 3 2 1 0 1 0 0 2 0 7 2 355
LYS 91 2P9V.subA 6 1 0 1 1 0 0 1 1 4 2 2 0 2 1 1 0 0 6 2 355
LYS 99 2P9V.subA 4 0 0 0 0 0 0 0 0 4 0 1 2 3 1 0 0 0 2 2 355
LYS 126 2P9V.subA 5 2 1 3 0 0 0 0 0 2 3 1 0 1 0 0 1 0 3 2 355
LYS 164 2P9V.subA 6 0 0 0 0 0 0 0 0 6 0 0 0 0 2 2 2 0 7 2 355
LYS 183 2P9V.subA 6 0 0 0 0 1 0 1 1 5 1 0 2 2 1 0 2 0 4 2 355
LYS 197 2P9V.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 1 0 0 0 11 2 355
LYS 207 2P9V.subA 5 0 1 1 0 1 0 1 1 3 2 0 0 0 2 0 1 0 3 2 355
LYS 224 2P9V.subA 12 0 1 1 0 0 0 0 0 11 1 1 1 2 1 3 5 0 8 2 355
LYS 239 2P9V.subA 7 1 1 2 0 1 0 1 1 4 3 0 1 1 1 0 2 0 5 2 355
LYS 246 2P9V.subA 4 0 1 1 0 0 0 0 0 3 1 1 2 3 0 0 0 0 4 2 355
LYS 290 2P9V.subA 2 0 0 0 0 0 0 0 0 2 0 0 1 1 0 0 1 0 0 2 355
LYS 299 2P9V.subA 3 0 0 0 0 0 0 0 0 3 0 0 1 1 1 0 1 0 4 2 355
LYS 332 2P9V.subA 6 0 2 2 0 0 0 0 0 4 2 0 1 1 0 1 2 0 2 2 355
LYS 342 2P9V.subA 9 0 0 0 1 0 0 1 1 8 1 0 3 3 2 2 1 0 7 2 355
LYS 5 2QF7.subA 10 1 1 2 1 0 0 1 1 7 3 3 0 3 2 0 2 0 4 2 1076
LYS 28 2QF7.subA 5 1 0 1 1 0 0 1 1 3 2 2 0 2 0 0 1 0 5 2 1076
LYS 38 2QF7.subA 6 1 1 2 0 1 0 1 1 3 3 0 0 0 1 1 1 0 1 2 1076
LYS 45 2QF7.subA 8 0 0 0 0 1 1 1 2 6 2 0 0 0 2 2 2 0 1 2 1076
LYS 79 2QF7.subA 6 0 0 0 0 1 0 1 1 5 1 0 0 0 2 0 3 0 4 2 1076
LYS 105 2QF7.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 4 2 1076
LYS 114 2QF7.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 5 2 1076
LYS 124 2QF7.subA 3 0 0 0 0 0 0 0 0 3 0 0 1 1 1 0 1 0 0 2 1076
LYS 209 2QF7.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 1 2 0 2 2 1076
LYS 245 2QF7.subA 11 0 4 4 0 0 1 0 1 6 5 0 3 3 1 0 2 0 6 2 1076
LYS 269 2QF7.subA 7 0 0 0 0 0 0 0 0 7 0 0 0 0 2 1 4 0 4 2 1076
LYS 292 2QF7.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 1 2 0 0 4 2 1076
LYS 319 2QF7.subA 7 1 1 2 0 0 1 0 1 4 3 0 0 0 2 0 2 0 3 2 1076
LYS 404 2QF7.subA 12 1 1 2 0 1 1 1 2 8 4 1 1 2 1 1 4 0 2 2 1076
LYS 446 2QF7.subA 7 1 0 1 0 0 1 0 1 5 2 1 0 1 1 3 0 0 4 2 1076
LYS 468 2QF7.subA 2 0 0 0 0 1 0 1 1 1 1 0 0 0 0 0 1 0 0 2 1076
LYS 475 2QF7.subA 10 2 0 2 0 1 0 1 1 7 3 3 1 4 1 1 1 0 1 2 1076
LYS 496 2QF7.subA 6 1 0 1 0 1 0 1 1 4 2 0 0 0 3 1 0 0 2 2 1076
LYS 515 2QF7.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 2 0 1 0 4 2 1076
LYS 519 2QF7.subA 11 2 0 2 1 1 0 2 2 7 4 1 2 3 2 0 2 0 6 2 1076
LYS 528 2QF7.subA 3 0 1 1 0 0 1 0 1 1 2 0 0 0 1 0 0 0 5 2 1076
LYS 529 2QF7.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 1 1 1 0 1 2 1076
LYS 538 2QF7.subA 4 0 2 2 0 1 0 1 1 1 3 0 1 1 0 0 0 0 2 2 1076
LYS 637 2QF7.subA 6 1 0 1 0 0 0 0 0 5 1 0 1 1 1 1 2 0 4 2 1076
LYS 645 2QF7.subA 5 0 0 0 0 1 0 1 1 4 1 0 1 1 2 1 0 0 3 2 1076
LYS 675 2QF7.subA 11 2 1 3 1 0 0 1 1 7 4 0 1 1 1 1 4 0 5 2 1076
LYS 694 2QF7.subA 7 1 0 1 0 1 0 1 1 5 2 1 0 1 3 1 0 0 3 2 1076
LYS 698 2QF7.subA 5 1 0 1 0 0 0 0 0 4 1 1 1 2 0 1 1 0 5 2 1076
LYS 709 2QF7.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 1 0 1 0 2 2 1076
LYS 725 2QF7.subA 7 1 0 1 0 0 0 0 0 6 1 0 0 0 3 1 2 0 6 2 1076
LYS 730 2QF7.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 4 0 1 0 6 2 1076
LYS 734 2QF7.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 2 1 1 0 7 2 1076
LYS 829 2QF7.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 2 2 1076
LYS 849 2QF7.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 0 2 1 0 2 2 1076
LYS 880 2QF7.subA 13 1 1 2 0 0 1 0 1 10 3 3 1 4 1 1 4 0 4 2 1076
LYS 886 2QF7.subA 6 1 0 1 0 0 0 0 0 5 1 2 1 3 1 0 1 0 4 2 1076
LYS 923 2QF7.subA 10 1 1 2 0 0 0 0 0 8 2 1 0 1 2 1 4 0 3 2 1076
LYS 939 2QF7.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 1 0 0 0 8 2 1076
LYS 940 2QF7.subA 11 1 1 2 0 1 0 1 1 8 3 1 0 1 2 1 4 0 3 2 1076
LYS 943 2QF7.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 0 1 0 3 2 1076
LYS 946 2QF7.subA 2 0 1 1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 0 2 1076
LYS 957 2QF7.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 0 0 2 0 1 2 1076
LYS 966 2QF7.subA 4 1 0 1 0 1 0 1 1 2 2 0 0 0 1 0 1 0 0 2 1076
LYS 970 2QF7.subA 2 0 1 1 0 0 0 0 0 1 1 0 0 0 0 0 1 0 2 2 1076
LYS 971 2QF7.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 0 0 3 0 4 2 1076
LYS 989 2QF7.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 2 2 1076
LYS 1028 2QF7.subA 2 0 1 1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 4 2 1076
LYS 1030 2QF7.subA 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 0 1 0 1 2 1076
LYS 1062 2QF7.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 0 0 4 0 2 2 1076
LYS 1078 2QF7.subA 3 0 0 0 0 2 0 2 2 1 2 0 0 0 1 0 0 0 0 2 1076
LYS 1128 2QF7.subA 6 1 1 2 0 0 0 0 0 4 2 0 0 0 3 1 0 0 0 2 1076
LYS 1138 2QF7.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 0 2 1076
LYS 1146 2QF7.subA 7 2 0 2 0 0 0 0 0 5 2 0 0 0 4 1 0 0 3 2 1076
LYS 25 2QPX.subA 6 2 0 2 0 0 0 0 0 3 2 0 1 1 1 0 1 0 8 2 376
LYS 42 2QPX.subA 2 2 0 2 0 0 0 0 0 0 2 0 0 0 0 0 0 0 3 2 376
LYS 50 2QPX.subA 8 2 0 2 0 0 0 0 0 6 2 1 1 2 2 1 1 0 7 2 376
LYS 63 2QPX.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 3 0 2 0 14 2 376
LYS 95 2QPX.subA 8 0 2 2 0 0 1 0 1 5 3 0 0 0 1 1 3 0 10 2 376
LYS 124 2QPX.subA 9 0 2 2 0 0 0 0 0 7 2 0 1 1 2 0 4 0 6 2 376
LYS 155 2QPX.subA 6 0 0 0 0 0 0 0 0 6 0 1 3 4 0 1 1 0 8 2 376
LYS 158 2QPX.subA 5 0 0 0 0 0 0 0 0 4 0 0 1 1 2 1 0 0 10 2 376
LYS 193 2QPX.subA 9 2 1 3 0 0 2 0 2 3 5 0 0 0 0 2 1 0 8 2 376
LYS 198 2QPX.subA 6 1 2 3 0 1 0 1 1 1 4 0 0 0 0 0 1 0 10 2 376
LYS 203 2QPX.subA 8 1 1 2 0 0 1 0 1 5 3 2 0 2 1 0 2 0 15 2 376
LYS 250 2QPX.subA 5 1 0 1 0 1 0 1 1 3 2 0 0 0 1 1 1 0 4 2 376
LYS 254 2QPX.subA 4 0 0 0 1 0 0 1 1 3 1 1 0 1 1 1 0 0 2 2 376
LYS 255 2QPX.subA 6 0 0 0 1 0 0 1 1 5 1 0 0 0 2 1 2 0 6 2 376
LYS 258 2QPX.subA 11 1 0 1 0 0 1 0 1 8 2 0 1 1 2 2 3 0 5 2 376
LYS 333 2QPX.subA 13 1 0 1 0 1 1 1 2 10 3 0 2 2 2 2 4 0 4 2 376
LYS 351 2QPX.subA 13 2 0 2 0 0 0 0 0 10 2 0 3 3 1 1 5 0 2 2 376
LYS 352 2QPX.subA 6 0 0 0 0 0 0 0 0 6 0 0 1 1 2 1 2 0 9 2 376
LYS 365 2QPX.subA 5 0 1 1 0 0 0 0 0 4 1 1 1 2 1 0 1 0 9 2 376
LYS 77 2R1K.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 8 2 329
LYS 82 2R1K.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 2 1 2 0 11 2 329
LYS 175 2R1K.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 9 2 329
LYS 185 2R1K.subA 8 0 2 2 0 1 0 1 1 5 3 0 0 0 1 0 4 0 3 2 329
LYS 285 2R1K.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 6 2 329
LYS 293 2R1K.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 0 1 2 0 9 2 329
LYS 339 2R1K.subA 7 0 1 1 1 0 0 1 1 5 2 0 0 0 1 1 3 0 3 2 329
LYS 77 2R1L.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 6 2 329
LYS 82 2R1L.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 2 1 2 0 11 2 329
LYS 175 2R1L.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 7 2 329
LYS 185 2R1L.subA 9 0 2 2 0 1 0 1 1 6 3 1 0 1 1 0 4 0 3 2 329
LYS 285 2R1L.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 3 2 329
LYS 293 2R1L.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 1 1 2 0 8 2 329
LYS 339 2R1L.subA 6 0 1 1 1 0 0 1 1 4 2 0 0 0 1 1 2 0 2 2 329
LYS 77 2R1M.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 2 2 328
LYS 82 2R1M.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 2 1 2 0 5 2 328
LYS 175 2R1M.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 0 2 328
LYS 185 2R1M.subA 9 0 2 2 0 1 0 1 1 6 3 1 0 1 1 0 4 0 3 2 328
LYS 285 2R1M.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 0 2 328
LYS 293 2R1M.subA 7 1 0 1 0 0 0 0 0 6 1 1 0 1 1 1 3 0 8 2 328
LYS 339 2R1M.subA 6 0 1 1 1 0 0 1 1 4 2 0 0 0 1 1 2 0 0 2 328
LYS 77 2R1N.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 6 2 328
LYS 82 2R1N.subA 9 1 1 2 0 1 0 1 1 6 3 0 0 0 3 1 2 0 9 2 328
LYS 175 2R1N.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 7 2 328
LYS 185 2R1N.subA 9 0 2 2 0 1 0 1 1 6 3 1 0 1 1 0 4 0 3 2 328
LYS 285 2R1N.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 2 2 328
LYS 293 2R1N.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 1 1 2 0 9 2 328
LYS 339 2R1N.subA 6 0 1 1 1 0 0 1 1 4 2 0 0 0 1 1 2 0 2 2 328
LYS 77 2R1P.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 7 2 328
LYS 82 2R1P.subA 9 1 1 2 0 1 0 1 1 6 3 0 0 0 3 1 2 0 12 2 328
LYS 175 2R1P.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 2 2 328
LYS 285 2R1P.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 2 2 328
LYS 293 2R1P.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 0 1 3 0 9 2 328
LYS 339 2R1P.subA 6 0 1 1 1 0 0 1 1 4 2 0 0 0 1 1 2 0 2 2 328
LYS 28 2RJG.subA 6 2 0 2 0 0 0 0 0 4 2 1 1 2 0 1 1 0 1 2 359
LYS 34 2RJG.subA 9 0 1 1 0 0 0 0 0 8 1 0 0 0 3 1 4 0 1 2 359
LYS 73 2RJG.subA 10 2 0 2 0 0 0 0 0 8 2 1 0 1 4 1 2 0 3 2 359
LYS 148 2RJG.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 1 0 0 0 2 2 359
LYS 167 2RJG.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 1 2 359
LYS 173 2RJG.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 1 2 359
LYS 185 2RJG.subA 8 0 2 2 0 0 0 0 0 6 2 0 0 0 3 2 1 0 2 2 359
LYS 248 2RJG.subA 4 1 1 2 0 0 1 0 1 1 3 0 0 0 1 0 0 0 1 2 359
LYS 317 2RJG.subA 5 2 0 2 0 0 0 0 0 3 2 0 1 1 1 0 1 0 3 2 359
LYS 339 2RJG.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 1 0 2 0 0 2 359
LYS 356 2RJG.subA 6 0 0 0 0 0 0 0 0 6 0 1 1 2 0 1 3 0 0 2 359
LYS 28 2RJH.subA 5 2 0 2 0 0 0 0 0 3 2 1 0 1 0 1 1 0 5 2 359
LYS 34 2RJH.subA 9 0 1 1 0 0 0 0 0 8 1 0 0 0 3 1 4 0 2 2 359
LYS 73 2RJH.subA 10 2 0 2 0 0 0 0 0 8 2 1 0 1 4 1 2 0 4 2 359
LYS 148 2RJH.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 1 0 0 0 0 2 359
LYS 167 2RJH.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 2 0 0 0 1 2 359
LYS 173 2RJH.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 2 2 359
LYS 185 2RJH.subA 8 0 2 2 0 0 0 0 0 6 2 0 0 0 3 2 1 0 1 2 359
LYS 248 2RJH.subA 4 1 1 2 0 0 1 0 1 1 3 0 0 0 1 0 0 0 2 2 359
LYS 317 2RJH.subA 4 2 0 2 0 0 0 0 0 2 2 0 1 1 1 0 0 0 3 2 359
LYS 339 2RJH.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 1 0 2 0 0 2 359
LYS 356 2RJH.subA 6 0 0 0 0 0 0 0 0 6 0 1 0 1 0 1 4 0 0 2 359
LYS 30 2RL3.subA 6 0 1 1 0 0 0 0 0 5 1 3 1 4 0 1 0 0 1 2 246
LYS 45 2RL3.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 1 1 1 0 3 2 246
LYS 49 2RL3.subA 3 0 0 0 0 0 0 0 0 3 0 3 0 3 0 0 0 0 0 2 246
LYS 61 2RL3.subA 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 6 2 246
LYS 84 2RL3.subA 5 0 0 0 0 0 0 0 0 5 0 0 2 2 1 0 2 0 5 2 246
LYS 91 2RL3.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 2 2 246
LYS 95 2RL3.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 1 2 246
LYS 100 2RL3.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 3 2 246
LYS 134 2RL3.subA 5 0 0 0 1 1 0 2 2 3 2 0 1 1 0 1 1 0 6 2 246
LYS 137 2RL3.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 10 2 246
LYS 138 2RL3.subA 4 0 0 0 1 0 0 1 1 3 1 1 0 1 0 2 0 0 5 2 246
LYS 152 2RL3.subA 7 1 1 2 0 0 0 0 0 5 2 0 1 1 2 1 1 0 4 2 246
LYS 177 2RL3.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 0 0 2 0 2 2 246
LYS 182 2RL3.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 0 1 0 3 2 246
LYS 189 2RL3.subA 10 0 1 1 0 0 1 0 1 8 2 0 2 2 0 2 4 0 2 2 246
LYS 205 2RL3.subA 11 0 0 0 0 0 0 0 0 11 0 5 0 5 1 0 4 0 4 2 246
LYS 228 2RL3.subA 10 0 4 4 0 0 0 0 0 6 4 0 0 0 1 2 3 0 1 2 246
LYS 246 2RL3.subA 7 0 1 1 0 0 0 0 0 6 1 1 1 2 0 0 4 0 6 2 246
LYS 251 2RL3.subA 6 0 0 0 0 1 0 1 1 5 1 3 0 3 1 1 0 0 6 2 246
LYS 256 2RL3.subA 6 0 2 2 0 0 0 0 0 4 2 2 0 2 0 0 2 0 5 2 246
LYS 39 2SFP.subA 9 0 1 1 0 0 0 0 0 8 1 0 0 0 3 1 4 0 1 2 379
LYS 76 2SFP.subA 7 1 1 2 0 1 0 1 1 4 3 0 0 0 2 0 2 0 0 2 379
LYS 140 2SFP.subA 4 2 0 2 0 0 0 0 0 2 2 0 0 0 1 0 1 0 2 2 379
LYS 146 2SFP.subA 5 0 2 2 0 1 0 1 1 2 3 1 0 1 0 1 0 0 0 2 379
LYS 234 2SFP.subA 7 1 0 1 0 0 0 0 0 6 1 1 0 1 3 0 2 0 1 2 379
LYS 242 2SFP.subA 7 1 1 2 0 0 0 0 0 5 2 0 2 2 2 0 1 0 1 2 379
LYS 255 2SFP.subA 5 0 0 0 0 0 1 0 1 4 1 1 0 1 0 0 3 0 1 2 379
LYS 256 2SFP.subA 5 0 2 2 0 0 0 0 0 3 2 0 1 1 0 1 1 0 1 2 379
LYS 262 2SFP.subA 6 0 1 1 0 0 0 0 0 5 1 1 0 1 2 1 1 0 0 2 379
LYS 303 2SFP.subA 8 0 1 1 0 0 1 0 1 6 2 0 1 1 1 0 4 0 0 2 379
LYS 328 2SFP.subA 9 0 0 0 0 1 1 1 2 7 2 3 1 4 1 1 1 0 3 2 379
LYS 372 2SFP.subA 3 1 0 1 0 1 1 1 2 0 3 0 0 0 0 0 0 0 1 2 379
LYS 6 2UAG.subA 9 2 1 3 0 0 0 0 0 6 3 0 2 2 1 1 2 0 3 2 429
LYS 45 2UAG.subA 5 1 0 1 0 0 0 0 0 4 1 0 0 0 2 0 2 0 3 2 429
LYS 115 2UAG.subA 14 0 1 1 0 0 0 0 0 13 1 5 2 7 2 0 3 2 5 2 429
LYS 127 2UAG.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 3 0 2 0 7 2 429
LYS 206 2UAG.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 0 2 0 7 2 429
LYS 251 2UAG.subA 5 1 0 1 0 1 0 1 1 3 2 0 0 0 1 1 1 0 2 2 429
LYS 254 2UAG.subA 5 0 1 1 0 1 0 1 1 3 2 0 0 0 0 1 2 0 1 2 429
LYS 259 2UAG.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 4 2 429
LYS 262 2UAG.subA 7 0 1 1 0 0 0 0 0 6 1 3 0 3 0 1 2 0 3 2 429
LYS 291 2UAG.subA 6 0 0 0 0 0 0 0 0 6 0 3 0 3 2 0 1 0 3 2 429
LYS 319 2UAG.subA 7 0 1 1 0 1 1 1 2 4 3 1 0 1 2 1 0 0 8 2 429
LYS 348 2UAG.subA 11 3 0 3 0 0 1 0 1 7 4 2 1 3 3 0 1 1 10 2 429
LYS 420 2UAG.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 1 0 0 0 2 429
LYS 434 2UAG.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 4 2 429
LYS 2 2UBP.subC 2 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 2 0 0 2 570
LYS 33 2UBP.subC 4 1 1 2 0 0 0 0 0 2 2 0 0 0 0 2 0 0 12 2 570
LYS 48 2UBP.subC 7 0 0 0 0 0 0 0 0 7 0 0 1 1 3 1 2 0 7 2 570
LYS 84 2UBP.subC 12 1 0 1 0 0 0 0 0 11 1 1 1 2 4 3 2 0 4 2 570
LYS 90 2UBP.subC 11 2 1 3 0 1 0 1 1 7 4 1 1 2 0 1 4 0 9 2 570
LYS 99 2UBP.subC 14 3 1 4 0 0 0 0 0 10 4 3 0 3 3 1 3 0 5 2 570
LYS 127 2UBP.subC 12 1 2 3 0 1 0 1 1 8 4 0 0 0 4 1 3 0 5 2 570
LYS 169 2UBP.subC 7 0 1 1 0 0 0 0 0 6 1 1 0 1 4 0 1 0 8 2 570
LYS 182 2UBP.subC 5 0 1 1 1 0 0 1 1 3 2 0 1 1 0 1 1 0 2 2 570
LYS 185 2UBP.subC 6 0 1 1 1 0 0 1 1 4 2 1 0 1 0 0 3 0 4 2 570
LYS 199 2UBP.subC 8 0 0 0 0 0 1 0 1 7 1 0 1 1 3 1 1 0 0 2 570
LYS 326 2UBP.subC 6 1 0 1 0 0 1 0 1 4 2 0 2 2 0 0 2 0 2 2 570
LYS 383 2UBP.subC 12 2 1 3 0 0 0 0 0 9 3 1 2 3 3 1 2 0 2 2 570
LYS 385 2UBP.subC 10 3 0 3 0 1 0 1 1 6 4 0 0 0 3 1 2 0 7 2 570
LYS 386 2UBP.subC 5 2 1 3 0 0 0 0 0 2 3 0 1 1 0 0 1 0 8 2 570
LYS 395 2UBP.subC 7 0 3 3 1 0 0 1 1 3 4 1 2 3 0 0 0 0 5 2 570
LYS 404 2UBP.subC 12 1 1 2 1 1 0 2 2 8 4 0 1 1 1 3 3 0 5 2 570
LYS 409 2UBP.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 1 2 4 0 0 2 570
LYS 431 2UBP.subC 12 0 3 3 0 0 1 0 1 8 4 1 0 1 3 1 3 0 7 2 570
LYS 441 2UBP.subC 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 1 0 0 4 2 570
LYS 446 2UBP.subC 7 0 0 0 0 0 0 0 0 7 0 0 1 1 2 3 1 0 2 2 570
LYS 452 2UBP.subC 14 1 0 1 0 0 0 0 0 13 1 2 2 4 3 1 5 0 2 2 570
LYS 497 2UBP.subC 8 1 0 1 1 0 0 1 1 6 2 3 1 4 1 0 1 0 11 2 570
LYS 507 2UBP.subC 5 0 1 1 0 0 0 0 0 4 1 0 0 0 2 0 2 0 12 2 570
LYS 511 2UBP.subC 5 0 1 1 0 1 1 1 2 2 3 0 0 0 1 0 1 0 10 2 570
LYS 518 2UBP.subC 7 1 0 1 1 0 1 1 2 4 3 1 1 2 1 0 1 0 11 2 570
LYS 525 2UBP.subC 16 1 0 1 0 0 1 0 1 14 2 2 1 3 6 0 5 0 3 2 570
LYS 526 2UBP.subC 2 1 0 1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 7 2 570
LYS 529 2UBP.subC 8 1 0 1 0 0 0 0 0 7 1 1 1 2 0 2 3 0 10 2 570
LYS 547 2UBP.subC 6 1 1 2 0 0 0 0 0 4 2 0 0 0 1 0 3 0 3 2 570
LYS 559 2UBP.subC 3 0 1 1 0 0 0 0 0 2 1 0 0 0 0 0 2 0 5 2 570
LYS 6 2UUO.subA 7 2 0 2 0 0 0 0 0 5 2 0 2 2 1 1 1 0 2 2 430
LYS 45 2UUO.subA 5 1 0 1 0 0 0 0 0 4 1 0 0 0 1 0 3 0 2 2 430
LYS 115 2UUO.subA 14 0 1 1 0 0 0 0 0 13 1 5 2 7 2 0 3 0 5 2 430
LYS 127 2UUO.subA 7 0 1 1 0 0 0 0 0 6 1 0 0 0 3 0 3 0 6 2 430
LYS 206 2UUO.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 0 2 0 2 2 430
LYS 251 2UUO.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 2 2 430
LYS 254 2UUO.subA 7 0 1 1 0 1 0 1 1 5 2 0 1 1 0 1 3 0 2 2 430
LYS 259 2UUO.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 3 2 430
LYS 262 2UUO.subA 7 0 1 1 0 0 0 0 0 6 1 3 0 3 0 1 2 0 2 2 430
LYS 291 2UUO.subA 7 0 0 0 0 0 0 0 0 7 0 3 0 3 3 0 1 0 2 2 430
LYS 319 2UUO.subA 10 1 0 1 0 1 1 1 2 7 3 2 0 2 2 1 2 0 2 2 430
LYS 348 2UUO.subA 12 3 0 3 0 0 1 0 1 8 4 2 1 3 4 0 1 0 7 2 430
LYS 420 2UUO.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 1 0 0 2 2 430
LYS 434 2UUO.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 1 2 430
LYS 6 2UUP.subA 7 2 0 2 0 0 0 0 0 5 2 0 2 2 1 1 1 0 3 2 440
LYS 45 2UUP.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 2 2 440
LYS 115 2UUP.subA 14 0 1 1 0 0 0 0 0 13 1 5 2 7 2 0 3 0 5 2 440
LYS 127 2UUP.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 3 0 2 0 6 2 440
LYS 206 2UUP.subA 5 0 0 0 0 0 0 0 0 5 0 0 1 1 2 0 2 0 5 2 440
LYS 251 2UUP.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 3 2 440
LYS 254 2UUP.subA 7 0 1 1 0 1 0 1 1 5 2 0 2 2 0 1 2 0 1 2 440
LYS 259 2UUP.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 5 2 440
LYS 262 2UUP.subA 7 0 1 1 0 0 0 0 0 6 1 3 0 3 0 1 2 0 5 2 440
LYS 291 2UUP.subA 5 0 0 0 0 0 0 0 0 5 0 2 0 2 2 0 1 0 2 2 440
LYS 319 2UUP.subA 10 1 1 2 0 1 1 1 2 6 4 2 0 2 2 1 1 0 7 2 440
LYS 348 2UUP.subA 12 3 0 3 0 0 1 0 1 8 4 2 1 3 4 0 1 0 8 2 440
LYS 420 2UUP.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 1 0 0 1 2 440
LYS 434 2UUP.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 4 2 440
LYS 2 2UYN.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 2 0 1 1 127
LYS 3 2UYN.subA 2 1 0 1 0 0 0 0 0 1 1 0 0 0 0 0 1 0 3 1 127
LYS 73 2UYN.subA 7 0 0 0 0 1 0 1 1 6 1 3 0 3 0 0 3 0 2 1 127
LYS 92 2UYN.subA 6 1 1 2 0 0 0 0 0 4 2 0 2 2 0 2 0 0 4 1 127
LYS 115 2UYN.subA 8 1 0 1 0 0 0 0 0 7 1 0 0 0 2 0 5 0 4 1 127
LYS 118 2UYN.subA 8 2 0 2 0 0 0 0 0 6 2 1 1 2 1 0 3 0 7 1 127
LYS 14 2V63.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 3 2 466
LYS 18 2V63.subA 5 1 0 1 0 0 0 0 0 4 1 2 0 2 0 0 2 0 6 2 466
LYS 81 2V63.subA 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 3 2 466
LYS 128 2V63.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 0 2 466
LYS 146 2V63.subA 5 0 0 0 0 0 0 0 0 4 0 1 0 1 2 0 1 0 2 2 466
LYS 161 2V63.subA 7 1 1 2 0 0 1 0 1 4 3 0 0 0 1 0 3 0 3 2 466
LYS 164 2V63.subA 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 4 2 466
LYS 175 2V63.subA 8 1 0 1 1 0 0 1 1 6 2 1 0 1 3 0 2 1 1 2 466
LYS 177 2V63.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 3 2 466
LYS 183 2V63.subA 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 3 2 466
LYS 227 2V63.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 2 0 0 1 2 466
LYS 236 2V63.subA 12 1 1 2 0 1 0 1 1 9 3 1 1 2 3 0 4 0 3 2 466
LYS 252 2V63.subA 5 0 2 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 7 2 466
LYS 258 2V63.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 2 1 2 0 5 2 466
LYS 316 2V63.subA 13 1 0 1 0 2 0 2 2 10 3 0 1 1 2 2 5 0 1 2 466
LYS 334 2V63.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 3 2 466
LYS 356 2V63.subA 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 4 2 466
LYS 450 2V63.subA 4 0 0 0 0 1 0 1 1 3 1 1 0 1 1 1 0 0 2 2 466
LYS 463 2V63.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 2 0 0 1 2 466
LYS 466 2V63.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 5 2 466
LYS 474 2V63.subA 5 2 0 2 0 1 0 1 1 2 3 1 0 1 0 0 1 0 2 2 466
LYS 14 2V67.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 4 2 466
LYS 18 2V67.subA 5 1 0 1 0 0 0 0 0 4 1 2 0 2 0 0 2 0 5 2 466
LYS 81 2V67.subA 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 3 2 466
LYS 128 2V67.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 1 0 0 0 2 466
LYS 146 2V67.subA 5 0 0 0 0 0 0 0 0 4 0 1 0 1 2 0 1 0 4 2 466
LYS 161 2V67.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 1 0 2 0 3 2 466
LYS 164 2V67.subA 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 4 2 466
LYS 175 2V67.subA 8 1 0 1 1 0 0 1 1 6 2 1 0 1 3 0 2 1 0 2 466
LYS 177 2V67.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 2 2 466
LYS 183 2V67.subA 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 4 2 466
LYS 227 2V67.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 2 0 0 3 2 466
LYS 236 2V67.subA 13 1 1 2 0 2 0 2 2 9 4 1 1 2 3 0 4 0 2 2 466
LYS 252 2V67.subA 5 0 2 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 6 2 466
LYS 258 2V67.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 2 1 2 0 6 2 466
LYS 316 2V67.subA 12 1 0 1 0 2 0 2 2 9 3 0 1 1 2 2 4 0 1 2 466
LYS 334 2V67.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 3 2 466
LYS 356 2V67.subA 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 3 2 466
LYS 450 2V67.subA 4 0 0 0 0 1 0 1 1 3 1 1 0 1 1 1 0 0 2 2 466
LYS 463 2V67.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 2 0 0 1 2 466
LYS 466 2V67.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 4 2 466
LYS 474 2V67.subA 6 2 0 2 0 1 0 1 1 3 3 1 0 1 0 0 2 0 2 2 466
LYS 14 2V68.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 3 2 466
LYS 18 2V68.subA 6 1 0 1 0 0 0 0 0 5 1 2 0 2 1 0 2 0 6 2 466
LYS 81 2V68.subA 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 2 2 466
LYS 128 2V68.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 1 0 0 0 2 466
LYS 146 2V68.subA 5 0 0 0 0 0 0 0 0 4 0 1 0 1 2 0 1 0 2 2 466
LYS 161 2V68.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 1 0 2 0 3 2 466
LYS 164 2V68.subA 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 2 2 466
LYS 175 2V68.subA 8 1 0 1 1 0 0 1 1 6 2 1 0 1 3 0 2 1 1 2 466
LYS 177 2V68.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 2 2 466
LYS 183 2V68.subA 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 4 2 466
LYS 227 2V68.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 2 0 0 5 2 466
LYS 236 2V68.subA 13 1 1 2 0 2 0 2 2 9 4 1 1 2 3 0 4 0 2 2 466
LYS 252 2V68.subA 5 0 2 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 5 2 466
LYS 258 2V68.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 2 1 2 0 5 2 466
LYS 316 2V68.subA 13 1 0 1 0 2 0 2 2 10 3 0 1 1 2 2 5 0 1 2 466
LYS 334 2V68.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 4 2 466
LYS 356 2V68.subA 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 4 2 466
LYS 450 2V68.subA 4 0 0 0 0 1 0 1 1 3 1 1 0 1 1 1 0 0 0 2 466
LYS 463 2V68.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 2 0 0 2 2 466
LYS 466 2V68.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 2 2 466
LYS 474 2V68.subA 6 2 0 2 0 1 0 1 1 3 3 1 0 1 0 0 2 0 1 2 466
LYS 14 2V69.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 0 2 459
LYS 18 2V69.subA 6 1 0 1 0 0 0 0 0 5 1 2 0 2 1 0 2 0 0 2 459
LYS 81 2V69.subA 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 0 2 459
LYS 128 2V69.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 1 0 0 0 2 459
LYS 146 2V69.subA 5 0 0 0 0 0 0 0 0 4 0 1 0 1 2 0 1 0 0 2 459
LYS 161 2V69.subA 7 1 1 2 0 0 1 0 1 4 3 0 0 0 1 0 3 0 0 2 459
LYS 164 2V69.subA 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 0 2 459
LYS 175 2V69.subA 8 1 0 1 1 0 0 1 1 6 2 1 0 1 3 0 2 1 0 2 459
LYS 177 2V69.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 1 2 459
LYS 183 2V69.subA 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 1 2 459
LYS 227 2V69.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 2 0 0 1 2 459
LYS 236 2V69.subA 12 1 1 2 0 1 0 1 1 9 3 1 1 2 3 0 4 0 1 2 459
LYS 252 2V69.subA 5 0 2 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 0 2 459
LYS 258 2V69.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 2 1 2 0 1 2 459
LYS 316 2V69.subA 12 1 0 1 0 2 0 2 2 9 3 0 1 1 2 2 4 0 0 2 459
LYS 334 2V69.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 0 2 459
LYS 356 2V69.subA 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 1 2 459
LYS 450 2V69.subA 4 0 0 0 0 1 0 1 1 3 1 1 0 1 1 1 0 0 0 2 459
LYS 463 2V69.subA 6 0 2 2 0 1 0 1 1 3 3 0 0 0 1 2 0 0 0 2 459
LYS 466 2V69.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 0 2 459
LYS 14 2V6A.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 5 2 467
LYS 18 2V6A.subA 6 1 0 1 0 0 0 0 0 5 1 2 0 2 1 0 2 0 9 2 467
LYS 81 2V6A.subA 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 7 2 467
LYS 128 2V6A.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 1 0 0 0 2 467
LYS 146 2V6A.subA 5 0 0 0 0 0 0 0 0 4 0 1 0 1 2 0 1 0 2 2 467
LYS 161 2V6A.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 1 0 2 0 3 2 467
LYS 164 2V6A.subA 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 4 2 467
LYS 175 2V6A.subA 8 1 0 1 1 0 0 1 1 6 2 1 0 1 3 0 2 1 0 2 467
LYS 177 2V6A.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 2 2 467
LYS 183 2V6A.subA 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 4 2 467
LYS 227 2V6A.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 2 0 0 2 2 467
LYS 236 2V6A.subA 12 1 1 2 0 1 0 1 1 9 3 1 1 2 3 0 4 0 2 2 467
LYS 252 2V6A.subA 5 0 2 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 8 2 467
LYS 258 2V6A.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 2 1 2 0 6 2 467
LYS 316 2V6A.subA 13 1 0 1 0 2 0 2 2 10 3 0 1 1 2 2 5 0 1 2 467
LYS 334 2V6A.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 2 2 467
LYS 356 2V6A.subA 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 4 2 467
LYS 450 2V6A.subA 4 0 0 0 0 1 0 1 1 3 1 1 0 1 1 1 0 0 1 2 467
LYS 463 2V6A.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 2 0 0 1 2 467
LYS 466 2V6A.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 6 2 467
LYS 474 2V6A.subA 6 2 0 2 0 1 0 1 1 3 3 1 0 1 0 0 2 0 2 2 467
LYS 8 2VC5.subA 8 1 0 1 0 0 0 0 0 7 1 0 1 1 2 0 4 0 1 2 314
LYS 14 2VC5.subA 5 1 1 2 1 1 0 2 2 1 4 1 0 1 0 0 0 0 0 2 314
LYS 54 2VC5.subA 6 0 0 0 0 1 0 1 1 5 1 1 1 2 1 0 2 0 0 2 314
LYS 62 2VC5.subA 8 1 0 1 0 0 0 0 0 7 1 1 1 2 2 1 2 0 1 2 314
LYS 81 2VC5.subA 7 1 2 3 0 1 0 1 1 3 4 0 0 0 0 2 1 0 0 2 314
LYS 84 2VC5.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 1 2 314
LYS 123 2VC5.subA 6 0 2 2 0 0 1 0 1 3 3 1 0 1 0 0 2 0 2 2 314
LYS 132 2VC5.subA 8 1 1 2 0 0 0 0 0 6 2 0 1 1 3 0 2 0 0 2 314
LYS 147 2VC5.subA 3 1 0 1 1 0 0 1 1 1 2 1 0 1 0 0 0 0 1 2 314
LYS 151 2VC5.subA 10 1 1 2 1 1 0 2 2 6 4 0 1 1 0 1 4 0 1 2 314
LYS 161 2VC5.subA 7 0 2 2 0 0 0 0 0 5 2 0 1 1 2 0 2 0 2 2 314
LYS 164 2VC5.subA 5 0 1 1 0 1 0 1 1 3 2 1 1 2 0 0 1 0 0 2 314
LYS 194 2VC5.subA 11 1 0 1 0 1 0 1 1 9 2 1 1 2 2 2 3 0 0 2 314
LYS 210 2VC5.subA 8 3 0 3 0 1 0 1 1 4 4 0 0 0 0 1 3 0 0 2 314
LYS 211 2VC5.subA 6 1 0 1 1 0 0 1 1 4 2 0 0 0 0 1 3 0 0 2 314
LYS 215 2VC5.subA 6 1 0 1 1 0 0 1 1 4 2 0 1 1 1 0 2 0 0 2 314
LYS 234 2VC5.subA 5 2 0 2 0 1 0 1 1 2 3 0 0 0 1 0 1 0 0 2 314
LYS 244 2VC5.subA 4 1 0 1 0 1 0 1 1 2 2 0 0 0 0 0 2 0 2 2 314
LYS 250 2VC5.subA 8 1 0 1 0 0 0 0 0 7 1 1 0 1 2 2 2 0 1 2 314
LYS 267 2VC5.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 3 2 0 0 1 2 314
LYS 271 2VC5.subA 8 0 0 0 0 1 0 1 1 7 1 0 0 0 4 3 0 0 1 2 314
LYS 273 2VC5.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 1 1 1 0 3 2 314
LYS 292 2VC5.subA 7 0 1 1 0 1 0 1 1 5 2 0 0 0 1 0 4 0 5 2 314
LYS 306 2VC5.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 2 2 0 0 2 2 314
LYS 310 2VC5.subA 6 0 1 1 1 0 0 1 1 4 2 2 0 2 1 0 1 0 0 2 314
LYS 311 2VC5.subA 6 1 1 2 0 0 0 0 0 4 2 1 1 2 0 2 0 0 0 2 314
LYS 8 2VC7.subA 8 1 0 1 0 0 0 0 0 7 1 0 1 1 2 0 4 0 4 2 314
LYS 14 2VC7.subA 4 1 1 2 1 0 0 1 1 1 3 1 0 1 0 0 0 0 2 2 314
LYS 54 2VC7.subA 6 0 0 0 0 1 0 1 1 5 1 1 1 2 1 0 2 0 3 2 314
LYS 62 2VC7.subA 8 1 0 1 0 0 0 0 0 7 1 1 1 2 2 1 2 0 6 2 314
LYS 81 2VC7.subA 6 0 2 2 0 1 0 1 1 3 3 0 0 0 0 2 1 0 3 2 314
LYS 84 2VC7.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 4 2 314
LYS 123 2VC7.subA 6 0 2 2 0 0 1 0 1 3 3 1 0 1 0 0 2 0 3 2 314
LYS 132 2VC7.subA 8 1 1 2 0 0 0 0 0 6 2 0 1 1 3 0 2 0 3 2 314
LYS 147 2VC7.subA 5 1 1 2 1 0 0 1 1 2 3 1 0 1 0 0 1 0 4 2 314
LYS 151 2VC7.subA 7 1 1 2 1 0 0 1 1 4 3 0 0 0 0 0 4 0 4 2 314
LYS 161 2VC7.subA 7 0 2 2 0 0 0 0 0 5 2 0 1 1 2 0 2 0 2 2 314
LYS 164 2VC7.subA 5 0 1 1 0 1 0 1 1 3 2 1 1 2 0 0 1 0 2 2 314
LYS 194 2VC7.subA 11 1 0 1 0 1 0 1 1 9 2 1 1 2 2 2 3 0 3 2 314
LYS 210 2VC7.subA 8 3 0 3 0 1 0 1 1 4 4 0 0 0 0 1 3 0 1 2 314
LYS 211 2VC7.subA 6 1 0 1 1 0 0 1 1 4 2 0 0 0 0 1 3 0 2 2 314
LYS 215 2VC7.subA 6 1 0 1 1 0 0 1 1 4 2 0 1 1 1 0 2 0 3 2 314
LYS 234 2VC7.subA 6 2 1 3 0 1 0 1 1 2 4 0 0 0 1 0 1 0 1 2 314
LYS 244 2VC7.subA 5 1 0 1 0 1 0 1 1 3 2 0 1 1 0 0 2 0 4 2 314
LYS 250 2VC7.subA 8 1 0 1 0 0 0 0 0 7 1 1 0 1 2 2 2 0 2 2 314
LYS 267 2VC7.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 3 2 0 0 1 2 314
LYS 271 2VC7.subA 10 0 1 1 0 1 0 1 1 8 2 1 0 1 4 3 0 0 2 2 314
LYS 273 2VC7.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 1 1 1 0 0 2 314
LYS 292 2VC7.subA 8 0 1 1 0 1 0 1 1 6 2 0 1 1 1 0 4 0 4 2 314
LYS 306 2VC7.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 2 2 0 0 4 2 314
LYS 310 2VC7.subA 7 0 1 1 1 0 0 1 1 5 2 2 1 3 1 0 1 0 3 2 314
LYS 311 2VC7.subA 6 1 1 2 0 0 0 0 0 4 2 1 1 2 0 2 0 0 2 2 314
LYS 14 2VDH.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 4 2 465
LYS 18 2VDH.subA 5 1 0 1 0 0 0 0 0 4 1 2 0 2 0 0 2 0 4 2 465
LYS 81 2VDH.subA 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 1 2 465
LYS 128 2VDH.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 0 2 465
LYS 146 2VDH.subA 5 0 0 0 0 0 0 0 0 4 0 1 0 1 2 0 1 0 1 2 465
LYS 161 2VDH.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 1 0 2 0 3 2 465
LYS 164 2VDH.subA 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 1 2 465
LYS 175 2VDH.subA 8 1 0 1 1 0 0 1 1 6 2 1 0 1 3 0 2 1 0 2 465
LYS 177 2VDH.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 1 2 465
LYS 183 2VDH.subA 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 3 2 465
LYS 227 2VDH.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 2 0 0 2 2 465
LYS 236 2VDH.subA 11 1 1 2 0 0 0 0 0 9 2 1 1 2 3 0 4 0 3 2 465
LYS 252 2VDH.subA 5 0 2 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 3 2 465
LYS 258 2VDH.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 2 1 2 0 4 2 465
LYS 316 2VDH.subA 12 1 0 1 0 2 0 2 2 9 3 0 1 1 2 2 4 0 1 2 465
LYS 334 2VDH.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 2 2 465
LYS 356 2VDH.subA 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 3 2 465
LYS 450 2VDH.subA 4 0 0 0 0 1 0 1 1 3 1 1 0 1 1 1 0 0 0 2 465
LYS 463 2VDH.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 2 0 0 2 2 465
LYS 466 2VDH.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 2 2 465
LYS 474 2VDH.subA 6 2 0 2 0 1 0 1 1 3 3 1 0 1 0 0 2 0 1 2 465
LYS 14 2VDI.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 2 2 465
LYS 18 2VDI.subA 6 1 0 1 0 0 0 0 0 5 1 2 0 2 1 0 2 0 3 2 465
LYS 81 2VDI.subA 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 1 2 465
LYS 128 2VDI.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 1 0 0 0 2 465
LYS 146 2VDI.subA 5 0 0 0 0 0 0 0 0 4 0 1 0 1 2 0 1 0 1 2 465
LYS 161 2VDI.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 1 0 2 0 2 2 465
LYS 164 2VDI.subA 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 0 2 465
LYS 175 2VDI.subA 8 1 0 1 1 0 0 1 1 6 2 1 0 1 3 0 2 1 0 2 465
LYS 177 2VDI.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 1 2 465
LYS 183 2VDI.subA 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 2 2 465
LYS 227 2VDI.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 2 0 0 3 2 465
LYS 236 2VDI.subA 12 1 1 2 0 1 0 1 1 9 3 1 1 2 3 0 4 0 2 2 465
LYS 252 2VDI.subA 5 0 2 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 4 2 465
LYS 258 2VDI.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 2 1 2 0 5 2 465
LYS 316 2VDI.subA 12 1 0 1 0 2 0 2 2 9 3 0 1 1 2 1 5 0 0 2 465
LYS 334 2VDI.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 2 2 465
LYS 356 2VDI.subA 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 1 2 465
LYS 450 2VDI.subA 4 0 0 0 0 1 0 1 1 3 1 1 0 1 1 1 0 0 0 2 465
LYS 463 2VDI.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 2 0 0 4 2 465
LYS 466 2VDI.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 2 2 465
LYS 474 2VDI.subA 6 2 0 2 0 1 0 1 1 3 3 1 0 1 0 0 2 0 0 2 465
LYS 6 2VTD.subA 7 2 0 2 0 0 0 0 0 5 2 0 2 2 1 1 1 0 5 2 435
LYS 45 2VTD.subA 5 1 0 1 0 0 0 0 0 4 1 0 0 0 1 0 3 0 3 2 435
LYS 115 2VTD.subA 14 0 1 1 0 0 0 0 0 13 1 5 2 7 2 0 3 0 5 2 435
LYS 127 2VTD.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 3 0 2 0 7 2 435
LYS 206 2VTD.subA 4 0 0 0 0 0 0 0 0 4 0 0 0 0 2 0 2 0 2 2 435
LYS 251 2VTD.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 6 2 435
LYS 254 2VTD.subA 7 0 1 1 0 1 0 1 1 5 2 0 1 1 0 1 3 0 3 2 435
LYS 259 2VTD.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 5 2 435
LYS 262 2VTD.subA 7 0 1 1 0 0 0 0 0 6 1 3 0 3 0 1 2 0 2 2 435
LYS 291 2VTD.subA 6 0 0 0 0 0 0 0 0 6 0 2 0 2 3 0 1 0 1 2 435
LYS 319 2VTD.subA 10 1 0 1 0 1 1 1 2 7 3 2 0 2 2 1 2 0 3 2 435
LYS 348 2VTD.subA 12 3 0 3 0 0 1 0 1 8 4 2 1 3 4 0 1 0 8 2 435
LYS 420 2VTD.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 1 0 0 2 2 435
LYS 434 2VTD.subA 5 0 1 1 0 1 0 1 1 3 2 0 0 0 2 0 1 0 3 2 435
LYS 6 2VTE.subA 7 2 0 2 0 0 0 0 0 5 2 0 2 2 1 1 1 0 0 2 434
LYS 45 2VTE.subA 5 1 0 1 0 0 0 0 0 4 1 0 0 0 2 0 2 0 3 2 434
LYS 115 2VTE.subA 13 0 1 1 0 0 0 0 0 12 1 4 2 6 2 0 3 0 1 2 434
LYS 127 2VTE.subA 7 0 1 1 0 0 0 0 0 6 1 0 0 0 3 0 3 0 3 2 434
LYS 206 2VTE.subA 4 0 0 0 0 0 0 0 0 4 0 0 0 0 2 0 2 0 2 2 434
LYS 251 2VTE.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 3 2 434
LYS 254 2VTE.subA 7 0 1 1 0 1 0 1 1 5 2 0 2 2 0 1 2 0 1 2 434
LYS 259 2VTE.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 3 2 434
LYS 262 2VTE.subA 6 0 0 0 0 0 0 0 0 6 0 3 0 3 0 1 2 0 2 2 434
LYS 291 2VTE.subA 7 0 0 0 0 0 0 0 0 7 0 3 0 3 3 0 1 0 0 2 434
LYS 319 2VTE.subA 9 1 0 1 0 1 1 1 2 6 3 2 1 3 1 2 0 0 4 2 434
LYS 348 2VTE.subA 12 3 0 3 0 0 1 0 1 8 4 2 1 3 4 0 1 0 4 2 434
LYS 420 2VTE.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 1 0 0 2 2 434
LYS 434 2VTE.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 1 2 434
LYS 30 2WGW.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 1 0 0 0 2 246
LYS 45 2WGW.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 1 1 1 0 6 2 246
LYS 49 2WGW.subA 3 0 0 0 0 0 0 0 0 3 0 3 0 3 0 0 0 0 0 2 246
LYS 61 2WGW.subA 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 9 2 246
LYS 70 2WGW.subA 12 0 0 0 1 0 0 1 1 11 1 4 0 4 2 2 2 0 2 2 246
LYS 84 2WGW.subA 5 0 0 0 0 0 0 0 0 5 0 0 2 2 1 0 2 0 6 2 246
LYS 91 2WGW.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 5 2 246
LYS 95 2WGW.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 1 2 246
LYS 100 2WGW.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 3 2 246
LYS 134 2WGW.subA 5 0 0 0 1 1 0 2 2 3 2 0 1 1 0 1 1 0 8 2 246
LYS 137 2WGW.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 5 2 246
LYS 138 2WGW.subA 4 0 0 0 1 0 0 1 1 3 1 1 0 1 0 2 0 0 5 2 246
LYS 152 2WGW.subA 7 1 1 2 0 0 0 0 0 5 2 0 1 1 2 1 1 0 4 2 246
LYS 177 2WGW.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 0 0 2 0 7 2 246
LYS 182 2WGW.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 0 1 0 5 2 246
LYS 189 2WGW.subA 10 0 1 1 0 0 1 0 1 8 2 0 2 2 0 2 4 0 2 2 246
LYS 205 2WGW.subA 11 0 0 0 0 0 0 0 0 11 0 5 0 5 1 0 4 0 4 2 246
LYS 228 2WGW.subA 10 0 4 4 0 0 0 0 0 6 4 0 0 0 1 2 3 0 3 2 246
LYS 246 2WGW.subA 6 0 0 0 0 0 0 0 0 6 0 1 1 2 0 0 4 0 8 2 246
LYS 251 2WGW.subA 6 0 0 0 0 1 0 1 1 5 1 3 0 3 1 1 0 0 8 2 246
LYS 256 2WGW.subA 6 0 2 2 0 0 0 0 0 4 2 2 0 2 0 0 2 0 6 2 246
LYS 6 2WJP.subA 7 2 0 2 0 0 0 0 0 5 2 0 2 2 1 1 1 0 4 2 435
LYS 45 2WJP.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 3 2 435
LYS 115 2WJP.subA 12 0 1 1 0 0 0 0 0 11 1 4 1 5 2 0 3 0 6 2 435
LYS 127 2WJP.subA 7 0 1 1 0 0 0 0 0 6 1 0 0 0 3 0 3 0 10 2 435
LYS 206 2WJP.subA 5 0 0 0 0 0 0 0 0 5 0 0 1 1 2 0 2 0 7 2 435
LYS 251 2WJP.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 4 2 435
LYS 254 2WJP.subA 7 0 1 1 0 1 0 1 1 5 2 0 2 2 0 1 2 0 2 2 435
LYS 259 2WJP.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 4 2 435
LYS 262 2WJP.subA 6 0 0 0 0 0 0 0 0 6 0 3 0 3 0 1 2 0 6 2 435
LYS 291 2WJP.subA 6 0 0 0 0 0 0 0 0 6 0 3 0 3 2 0 1 0 3 2 435
LYS 319 2WJP.subA 7 1 0 1 0 1 1 1 2 4 3 2 0 2 1 1 0 0 6 2 435
LYS 348 2WJP.subA 12 3 0 3 0 0 1 0 1 8 4 2 1 3 4 0 1 0 8 2 435
LYS 420 2WJP.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 1 0 0 4 2 435
LYS 434 2WJP.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 6 2 435
LYS 157 2WTZ.subA 13 0 1 1 0 0 0 0 0 12 1 7 0 7 2 0 2 1 0 2 504
LYS 396 2WTZ.subA 6 1 1 2 0 0 2 0 2 2 4 0 0 0 2 0 0 0 0 2 504
LYS 428 2WTZ.subA 9 1 0 1 0 1 1 1 2 6 3 0 0 0 5 0 1 0 0 2 504
LYS 503 2WTZ.subA 12 3 1 4 0 1 0 1 1 7 5 1 2 3 3 1 0 0 0 2 504
LYS 30 2X01.subA 5 0 1 1 0 0 0 0 0 4 1 1 1 2 1 1 0 0 1 2 244
LYS 45 2X01.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 1 1 1 0 1 2 244
LYS 49 2X01.subA 3 0 0 0 0 0 0 0 0 3 0 3 0 3 0 0 0 0 2 2 244
LYS 61 2X01.subA 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 3 2 244
LYS 70 2X01.subA 12 0 0 0 1 0 0 1 1 11 1 2 0 2 3 2 3 0 3 2 244
LYS 84 2X01.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 4 2 244
LYS 91 2X01.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 3 2 244
LYS 95 2X01.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 1 2 244
LYS 100 2X01.subA 4 0 1 1 0 1 0 1 1 2 2 0 1 1 0 0 1 0 2 2 244
LYS 134 2X01.subA 5 0 0 0 1 1 0 2 2 3 2 0 1 1 0 1 1 0 5 2 244
LYS 137 2X01.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 5 2 244
LYS 138 2X01.subA 4 0 0 0 1 0 0 1 1 3 1 1 0 1 0 2 0 0 4 2 244
LYS 152 2X01.subA 7 1 1 2 0 0 0 0 0 5 2 0 1 1 2 1 1 0 1 2 244
LYS 177 2X01.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 0 0 2 0 1 2 244
LYS 182 2X01.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 0 1 0 3 2 244
LYS 189 2X01.subA 10 0 1 1 0 0 1 0 1 8 2 0 2 2 0 2 4 0 2 2 244
LYS 205 2X01.subA 11 0 0 0 0 0 0 0 0 11 0 4 0 4 2 0 4 0 1 2 244
LYS 228 2X01.subA 10 0 4 4 0 0 0 0 0 6 4 0 0 0 1 2 3 0 0 2 244
LYS 246 2X01.subA 8 0 1 1 0 0 0 0 0 7 1 1 1 2 1 0 4 0 1 2 244
LYS 251 2X01.subA 6 0 0 0 0 1 0 1 1 5 1 3 0 3 1 1 0 0 3 2 244
LYS 256 2X01.subA 6 0 2 2 0 0 0 0 0 4 2 2 0 2 0 0 2 0 2 2 244
LYS 30 2X02.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 1 0 0 0 2 245
LYS 45 2X02.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 1 1 1 0 7 2 245
LYS 49 2X02.subA 3 0 0 0 0 0 0 0 0 3 0 3 0 3 0 0 0 0 3 2 245
LYS 61 2X02.subA 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 10 2 245
LYS 84 2X02.subA 5 0 0 0 0 0 0 0 0 5 0 0 2 2 1 0 2 0 6 2 245
LYS 91 2X02.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 9 2 245
LYS 95 2X02.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 1 2 245
LYS 100 2X02.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 6 2 245
LYS 134 2X02.subA 5 0 0 0 1 1 0 2 2 3 2 0 1 1 0 1 1 0 6 2 245
LYS 137 2X02.subA 4 0 0 0 0 0 0 0 0 4 0 0 2 2 0 1 1 0 7 2 245
LYS 138 2X02.subA 4 0 0 0 1 0 0 1 1 3 1 1 0 1 0 2 0 0 6 2 245
LYS 152 2X02.subA 7 1 1 2 0 0 0 0 0 5 2 0 1 1 2 1 1 0 8 2 245
LYS 177 2X02.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 0 0 2 0 5 2 245
LYS 182 2X02.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 0 1 0 6 2 245
LYS 189 2X02.subA 10 0 1 1 0 0 1 0 1 8 2 0 2 2 0 2 4 0 2 2 245
LYS 205 2X02.subA 11 0 0 0 0 0 0 0 0 11 0 5 0 5 1 0 4 0 5 2 245
LYS 228 2X02.subA 9 0 4 4 0 0 0 0 0 5 4 0 0 0 1 2 2 0 3 2 245
LYS 246 2X02.subA 6 0 0 0 0 0 0 0 0 6 0 1 1 2 0 0 4 0 8 2 245
LYS 251 2X02.subA 6 0 0 0 0 1 0 1 1 5 1 3 0 3 1 1 0 0 10 2 245
LYS 256 2X02.subA 6 0 2 2 0 0 0 0 0 4 2 2 0 2 0 0 2 0 7 2 245
LYS 6 2X5O.subA 7 2 0 2 0 0 0 0 0 5 2 0 2 2 1 1 1 0 6 2 439
LYS 45 2X5O.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 0 2 439
LYS 115 2X5O.subA 13 0 1 1 0 0 0 0 0 12 1 4 2 6 2 0 3 0 7 2 439
LYS 127 2X5O.subA 7 0 1 1 0 0 0 0 0 6 1 0 0 0 3 0 3 0 10 2 439
LYS 206 2X5O.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 0 2 0 5 2 439
LYS 251 2X5O.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 1 1 0 4 2 439
LYS 254 2X5O.subA 7 0 1 1 0 1 0 1 1 5 2 0 2 2 0 1 2 0 2 2 439
LYS 259 2X5O.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 3 2 439
LYS 262 2X5O.subA 6 0 0 0 0 0 0 0 0 6 0 3 0 3 0 1 2 0 5 2 439
LYS 291 2X5O.subA 6 0 0 0 0 0 0 0 0 6 0 3 0 3 2 0 1 0 1 2 439
LYS 319 2X5O.subA 9 1 1 2 0 1 1 1 2 5 4 3 0 3 1 1 0 0 8 2 439
LYS 348 2X5O.subA 12 3 0 3 0 0 1 0 1 8 4 2 1 3 4 0 1 0 8 2 439
LYS 420 2X5O.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 1 0 0 3 2 439
LYS 434 2X5O.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 8 2 439
LYS 157 2XJA.subA 12 0 1 1 0 0 0 0 0 11 1 6 0 6 2 0 2 2 1 2 502
LYS 396 2XJA.subA 5 1 1 2 0 0 1 0 1 2 3 0 0 0 2 0 0 0 0 2 502
LYS 428 2XJA.subA 9 1 0 1 0 2 1 2 3 5 4 0 0 0 4 0 1 0 0 2 502
LYS 503 2XJA.subA 11 3 1 4 0 1 0 1 1 6 5 0 2 2 3 1 0 0 0 2 502
LYS 6 2XPC.subA 7 2 0 2 0 0 0 0 0 5 2 0 2 2 1 1 1 0 6 2 435
LYS 45 2XPC.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 0 2 435
LYS 115 2XPC.subA 14 0 1 1 0 0 0 0 0 13 1 5 2 7 2 0 3 0 7 2 435
LYS 127 2XPC.subA 7 0 1 1 0 0 0 0 0 6 1 0 0 0 3 0 3 0 8 2 435
LYS 206 2XPC.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 0 2 0 6 2 435
LYS 251 2XPC.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 2 2 435
LYS 254 2XPC.subA 8 0 1 1 0 1 0 1 1 6 2 0 2 2 0 1 3 0 2 2 435
LYS 259 2XPC.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 7 2 435
LYS 262 2XPC.subA 7 0 1 1 0 0 0 0 0 6 1 3 0 3 0 1 2 0 5 2 435
LYS 291 2XPC.subA 6 0 0 0 0 0 0 0 0 6 0 3 0 3 2 0 1 0 1 2 435
LYS 319 2XPC.subA 10 1 0 1 0 1 1 1 2 7 3 2 0 2 2 1 2 0 5 2 435
LYS 348 2XPC.subA 12 3 0 3 0 0 1 0 1 8 4 2 1 3 4 0 1 0 8 2 435
LYS 420 2XPC.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 1 0 0 4 2 435
LYS 434 2XPC.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 7 2 435
LYS 6 2Y66.subA 7 2 0 2 0 0 0 0 0 5 2 0 2 2 1 1 1 0 3 2 434
LYS 45 2Y66.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 2 2 434
LYS 115 2Y66.subA 13 0 1 1 0 0 0 0 0 12 1 5 1 6 2 0 3 0 4 2 434
LYS 127 2Y66.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 3 0 2 0 7 2 434
LYS 206 2Y66.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 0 2 0 5 2 434
LYS 251 2Y66.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 3 2 434
LYS 254 2Y66.subA 7 0 1 1 0 1 0 1 1 5 2 0 2 2 0 1 2 0 2 2 434
LYS 259 2Y66.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 4 2 434
LYS 262 2Y66.subA 6 0 0 0 0 0 0 0 0 6 0 3 0 3 0 1 2 0 4 2 434
LYS 291 2Y66.subA 6 0 0 0 0 0 0 0 0 6 0 3 0 3 2 0 1 0 3 2 434
LYS 319 2Y66.subA 9 1 1 2 0 1 1 1 2 5 4 3 0 3 1 1 0 0 4 2 434
LYS 348 2Y66.subA 12 3 0 3 0 0 1 0 1 8 4 2 1 3 4 0 1 0 6 2 434
LYS 420 2Y66.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 1 0 0 1 2 434
LYS 434 2Y66.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 6 2 434
LYS 6 2Y67.subA 7 2 0 2 0 0 0 0 0 5 2 0 2 2 1 1 1 0 1 2 435
LYS 45 2Y67.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 1 2 435
LYS 115 2Y67.subA 14 0 1 1 0 0 0 0 0 13 1 5 2 7 2 0 3 0 4 2 435
LYS 127 2Y67.subA 7 0 1 1 0 0 0 0 0 6 1 0 0 0 3 0 3 0 4 2 435
LYS 206 2Y67.subA 4 0 0 0 0 0 0 0 0 4 0 0 0 0 2 0 2 0 3 2 435
LYS 251 2Y67.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 1 2 435
LYS 254 2Y67.subA 8 0 1 1 0 1 0 1 1 6 2 0 2 2 0 1 3 0 2 2 435
LYS 259 2Y67.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 1 2 435
LYS 262 2Y67.subA 6 0 0 0 0 0 0 0 0 6 0 3 0 3 0 1 2 0 2 2 435
LYS 291 2Y67.subA 6 0 0 0 0 0 0 0 0 6 0 3 0 3 2 0 1 0 0 2 435
LYS 319 2Y67.subA 8 1 0 1 0 1 1 1 2 5 3 3 0 3 1 1 0 0 2 2 435
LYS 348 2Y67.subA 11 3 0 3 0 0 1 0 1 7 4 2 1 3 3 0 1 0 6 2 435
LYS 420 2Y67.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 1 0 0 0 2 435
LYS 434 2Y67.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 1 2 435
LYS 6 2Y68.subA 7 2 0 2 0 0 0 0 0 5 2 0 2 2 1 1 1 0 6 2 435
LYS 45 2Y68.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 1 2 435
LYS 115 2Y68.subA 13 0 1 1 0 0 0 0 0 12 1 4 2 6 2 0 3 0 6 2 435
LYS 127 2Y68.subA 7 0 1 1 0 0 0 0 0 6 1 0 0 0 3 0 3 0 9 2 435
LYS 206 2Y68.subA 6 0 0 0 0 0 0 0 0 6 0 0 2 2 2 0 2 0 7 2 435
LYS 251 2Y68.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 3 2 435
LYS 254 2Y68.subA 7 0 1 1 0 1 0 1 1 5 2 0 2 2 0 1 2 0 0 2 435
LYS 259 2Y68.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 3 2 435
LYS 262 2Y68.subA 6 0 0 0 0 0 0 0 0 6 0 3 0 3 0 1 2 0 5 2 435
LYS 291 2Y68.subA 6 0 0 0 0 0 0 0 0 6 0 3 0 3 2 0 1 0 2 2 435
LYS 319 2Y68.subA 7 1 0 1 0 1 1 1 2 4 3 2 0 2 1 1 0 0 8 2 435
LYS 348 2Y68.subA 12 3 0 3 0 0 1 0 1 8 4 2 1 3 4 0 1 0 9 2 435
LYS 420 2Y68.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 1 0 0 1 2 435
LYS 434 2Y68.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 4 2 435
LYS 8 2Z24.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 0 1 4 0 2 2 343
LYS 26 2Z24.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 8 2 343
LYS 131 2Z24.subA 4 0 2 2 0 1 0 1 1 1 3 0 0 0 0 0 1 0 3 2 343
LYS 172 2Z24.subA 8 0 0 0 0 1 0 1 1 7 1 0 0 0 3 1 3 0 8 2 343
LYS 181 2Z24.subA 6 2 0 2 0 0 1 0 1 3 3 1 0 1 1 0 1 0 4 2 343
LYS 226 2Z24.subA 11 0 2 2 0 1 1 1 2 7 4 2 0 2 1 1 3 0 6 2 343
LYS 259 2Z24.subA 12 0 1 1 0 2 2 2 4 7 5 0 0 0 6 1 0 0 3 2 343
LYS 346 2Z24.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 0 0 2 0 1 2 343
LYS 8 2Z25.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 0 1 4 0 5 2 343
LYS 26 2Z25.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 11 2 343
LYS 131 2Z25.subA 4 0 2 2 0 1 0 1 1 1 3 0 0 0 0 0 1 0 3 2 343
LYS 172 2Z25.subA 8 0 0 0 0 1 0 1 1 7 1 0 0 0 3 1 3 0 8 2 343
LYS 181 2Z25.subA 6 2 0 2 0 0 1 0 1 3 3 1 0 1 1 0 1 0 5 2 343
LYS 226 2Z25.subA 11 0 2 2 0 1 1 1 2 7 4 2 0 2 1 1 3 0 6 2 343
LYS 259 2Z25.subA 12 0 1 1 0 2 2 2 4 7 5 0 0 0 6 1 0 0 3 2 343
LYS 346 2Z25.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 0 0 3 0 2 2 343
LYS 8 2Z26.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 0 1 3 0 3 2 344
LYS 26 2Z26.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 12 2 344
LYS 131 2Z26.subA 4 0 2 2 0 1 0 1 1 1 3 0 0 0 0 0 1 0 3 2 344
LYS 172 2Z26.subA 8 0 0 0 0 1 0 1 1 7 1 0 0 0 3 1 3 0 10 2 344
LYS 181 2Z26.subA 6 2 0 2 0 0 1 0 1 3 3 1 0 1 1 0 1 0 5 2 344
LYS 226 2Z26.subA 11 0 2 2 0 1 1 1 2 7 4 2 0 2 1 1 3 0 6 2 344
LYS 259 2Z26.subA 12 0 1 1 0 2 2 2 4 7 5 0 0 0 6 1 0 0 3 2 344
LYS 346 2Z26.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 0 0 2 0 1 2 344
LYS 8 2Z27.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 0 1 4 0 3 2 343
LYS 26 2Z27.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 10 2 343
LYS 131 2Z27.subA 4 0 2 2 0 1 0 1 1 1 3 0 0 0 0 0 1 0 2 2 343
LYS 172 2Z27.subA 8 0 0 0 0 1 0 1 1 7 1 0 0 0 3 1 3 0 10 2 343
LYS 181 2Z27.subA 6 2 0 2 0 0 1 0 1 3 3 1 0 1 1 0 1 0 4 2 343
LYS 226 2Z27.subA 11 0 2 2 0 1 1 1 2 7 4 2 0 2 1 1 3 0 6 2 343
LYS 259 2Z27.subA 12 0 1 1 0 2 2 2 4 7 5 0 0 0 6 1 0 0 3 2 343
LYS 346 2Z27.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 0 0 3 0 1 2 343
LYS 8 2Z28.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 0 1 4 0 2 2 343
LYS 26 2Z28.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 9 2 343
LYS 131 2Z28.subA 4 0 2 2 0 1 0 1 1 1 3 0 0 0 0 0 1 0 2 2 343
LYS 172 2Z28.subA 8 0 0 0 0 1 0 1 1 7 1 0 0 0 3 1 3 0 6 2 343
LYS 181 2Z28.subA 6 2 0 2 0 0 1 0 1 3 3 1 0 1 1 0 1 0 4 2 343
LYS 226 2Z28.subA 11 0 2 2 0 1 1 1 2 7 4 2 0 2 1 1 3 0 6 2 343
LYS 259 2Z28.subA 12 0 1 1 0 2 2 2 4 7 5 0 0 0 6 1 0 0 3 2 343
LYS 346 2Z28.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 0 0 2 0 1 2 343
LYS 8 2Z29.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 0 1 4 0 3 2 343
LYS 26 2Z29.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 9 2 343
LYS 131 2Z29.subA 4 0 2 2 0 1 0 1 1 1 3 0 0 0 0 0 1 0 2 2 343
LYS 172 2Z29.subA 8 0 0 0 0 1 0 1 1 7 1 0 0 0 3 1 3 0 9 2 343
LYS 181 2Z29.subA 6 2 0 2 0 0 1 0 1 3 3 1 0 1 1 0 1 0 3 2 343
LYS 226 2Z29.subA 11 0 2 2 0 1 1 1 2 7 4 2 0 2 1 1 3 0 6 2 343
LYS 259 2Z29.subA 12 0 1 1 0 2 2 2 4 7 5 0 0 0 6 1 0 0 3 2 343
LYS 346 2Z29.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 0 0 3 0 1 2 343
LYS 8 2Z2A.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 0 1 4 0 0 2 343
LYS 26 2Z2A.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 9 2 343
LYS 131 2Z2A.subA 4 0 2 2 0 1 0 1 1 1 3 0 0 0 0 0 1 0 3 2 343
LYS 172 2Z2A.subA 8 0 0 0 0 1 0 1 1 7 1 0 0 0 3 1 3 0 7 2 343
LYS 181 2Z2A.subA 6 2 0 2 0 0 1 0 1 3 3 1 0 1 1 0 1 0 4 2 343
LYS 226 2Z2A.subA 11 0 2 2 0 1 1 1 2 7 4 2 0 2 1 1 3 0 6 2 343
LYS 259 2Z2A.subA 12 0 1 1 0 2 2 2 4 7 5 0 0 0 6 1 0 0 3 2 343
LYS 346 2Z2A.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 0 0 2 0 1 2 343
LYS 8 2Z2B.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 0 1 3 0 0 2 331
LYS 26 2Z2B.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 0 2 331
LYS 121 2Z2B.subA 4 0 2 2 0 1 0 1 1 1 3 0 0 0 0 0 1 0 2 2 331
LYS 162 2Z2B.subA 8 0 0 0 0 1 0 1 1 7 1 0 0 0 3 1 3 0 7 2 331
LYS 171 2Z2B.subA 6 2 0 2 0 0 1 0 1 3 3 1 0 1 1 0 1 0 2 2 331
LYS 216 2Z2B.subA 10 0 1 1 0 1 1 1 2 7 3 2 0 2 1 1 3 0 6 2 331
LYS 249 2Z2B.subA 12 0 1 1 0 2 2 2 4 7 5 0 0 0 6 1 0 0 0 2 331
LYS 336 2Z2B.subA 4 0 1 1 0 0 0 0 0 3 1 1 0 1 0 0 2 0 0 2 331
LYS 208 2ZC1.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 0 3 0 0 11 2 333
LYS 385 2ZC1.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 2 2 333
LYS 15 3A12.subA 6 1 1 2 0 0 0 0 0 4 2 0 0 0 1 2 1 0 4 2 437
LYS 21 3A12.subA 3 0 0 0 0 0 0 0 0 3 0 1 0 1 1 1 0 0 1 2 437
LYS 22 3A12.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 0 1 0 0 3 2 437
LYS 73 3A12.subA 10 1 0 1 0 0 0 0 0 9 1 1 0 1 4 3 1 0 2 2 437
LYS 116 3A12.subA 3 0 0 0 0 1 0 1 1 2 1 0 0 0 1 0 1 0 1 2 437
LYS 119 3A12.subA 6 0 1 1 0 0 0 0 0 5 1 1 0 1 3 0 1 0 2 2 437
LYS 131 3A12.subA 6 0 3 3 0 0 0 0 0 3 3 0 0 0 1 1 1 0 4 2 437
LYS 148 3A12.subA 4 0 1 1 0 1 0 1 1 2 2 0 0 0 1 0 1 0 1 2 437
LYS 153 3A12.subA 2 1 0 1 0 0 0 0 0 1 1 0 0 0 0 0 1 0 1 2 437
LYS 163 3A12.subA 7 1 0 1 1 0 0 1 1 5 2 0 0 0 4 0 1 1 1 2 437
LYS 165 3A12.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 3 2 437
LYS 175 3A12.subA 4 0 3 3 0 0 0 0 0 1 3 0 0 0 0 0 1 0 0 2 437
LYS 211 3A12.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 0 2 0 1 2 437
LYS 215 3A12.subA 5 1 1 2 0 0 0 0 0 3 2 0 0 0 0 1 2 0 1 2 437
LYS 223 3A12.subA 5 0 3 3 0 0 0 0 0 2 3 0 0 0 0 0 2 0 3 2 437
LYS 224 3A12.subA 11 1 2 3 0 1 0 1 1 7 4 2 0 2 2 0 3 0 2 2 437
LYS 250 3A12.subA 8 0 1 1 0 0 1 0 1 6 2 1 0 1 1 0 4 0 1 2 437
LYS 303 3A12.subA 12 1 0 1 0 0 0 0 0 11 1 0 1 1 1 4 5 0 1 2 437
LYS 322 3A12.subA 6 0 0 0 0 0 0 0 0 6 0 0 0 0 4 0 2 0 4 2 437
LYS 327 3A12.subA 6 0 0 0 0 0 0 0 0 6 0 0 0 0 4 1 1 0 5 2 437
LYS 343 3A12.subA 5 0 1 1 0 1 1 1 2 2 3 0 0 0 1 1 0 0 2 2 437
LYS 355 3A12.subA 5 0 0 0 0 0 1 0 1 4 1 1 1 2 0 2 0 0 3 2 437
LYS 360 3A12.subA 7 0 0 0 0 0 0 0 0 7 0 1 0 1 4 1 1 0 4 2 437
LYS 426 3A12.subA 5 1 1 2 0 0 0 0 0 3 2 1 0 1 1 1 0 0 0 2 437
LYS 429 3A12.subA 2 0 1 1 0 0 1 0 1 0 2 0 0 0 0 0 0 0 1 2 437
LYS 437 3A12.subA 7 0 1 1 0 1 1 1 2 4 3 0 0 0 3 1 0 0 1 2 437
LYS 15 3A13.subA 6 1 1 2 0 0 0 0 0 4 2 0 0 0 1 2 1 0 4 2 436
LYS 21 3A13.subA 5 0 1 1 1 0 0 1 1 3 2 1 0 1 1 1 0 0 3 2 436
LYS 22 3A13.subA 5 0 1 1 1 1 0 2 2 2 3 1 0 1 0 1 0 0 3 2 436
LYS 73 3A13.subA 10 1 0 1 0 0 0 0 0 9 1 1 0 1 3 4 1 0 3 2 436
LYS 116 3A13.subA 3 0 0 0 0 1 0 1 1 2 1 0 0 0 1 0 1 0 0 2 436
LYS 119 3A13.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 3 0 1 0 3 2 436
LYS 131 3A13.subA 6 0 3 3 0 0 0 0 0 3 3 0 0 0 1 1 1 0 6 2 436
LYS 148 3A13.subA 5 0 1 1 0 1 0 1 1 3 2 0 0 0 1 0 2 0 3 2 436
LYS 153 3A13.subA 2 1 0 1 0 0 0 0 0 1 1 0 0 0 0 0 1 0 3 2 436
LYS 163 3A13.subA 9 1 0 1 2 0 0 2 2 6 3 0 0 0 4 0 2 1 2 2 436
LYS 165 3A13.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 3 2 436
LYS 175 3A13.subA 5 1 3 4 0 0 0 0 0 1 4 0 0 0 0 0 1 0 1 2 436
LYS 211 3A13.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 0 2 0 2 2 436
LYS 215 3A13.subA 5 1 1 2 0 0 0 0 0 3 2 0 0 0 0 1 2 0 1 2 436
LYS 223 3A13.subA 6 0 3 3 1 0 0 1 1 2 4 0 0 0 0 0 2 0 2 2 436
LYS 224 3A13.subA 11 1 2 3 0 1 0 1 1 7 4 2 0 2 2 0 3 0 2 2 436
LYS 250 3A13.subA 9 0 1 1 1 0 1 1 2 6 3 1 0 1 1 0 4 0 2 2 436
LYS 303 3A13.subA 12 1 0 1 0 0 0 0 0 11 1 0 1 1 1 4 5 0 1 2 436
LYS 322 3A13.subA 7 0 0 0 1 0 0 1 1 6 1 0 0 0 4 0 2 0 2 2 436
LYS 343 3A13.subA 4 0 1 1 0 0 1 0 1 2 2 0 0 0 1 1 0 0 1 2 436
LYS 355 3A13.subA 5 0 0 0 0 0 1 0 1 4 1 1 1 2 0 2 0 0 2 2 436
LYS 360 3A13.subA 7 0 0 0 0 0 0 0 0 7 0 1 0 1 4 1 1 0 4 2 436
LYS 426 3A13.subA 4 1 1 2 0 0 0 0 0 2 2 1 0 1 1 0 0 0 0 2 436
LYS 429 3A13.subA 2 0 1 1 0 0 1 0 1 0 2 0 0 0 0 0 0 0 0 2 436
LYS 437 3A13.subA 7 0 1 1 0 1 1 1 2 4 3 0 0 0 3 1 0 0 0 2 436
LYS 77 3A3W.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 6 2 329
LYS 82 3A3W.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 2 1 2 0 13 2 329
LYS 175 3A3W.subA 3 0 0 0 0 0 0 0 0 3 0 0 1 1 2 0 0 0 0 2 329
LYS 285 3A3W.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 3 2 329
LYS 293 3A3W.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 1 1 2 0 8 2 329
LYS 339 3A3W.subA 7 0 1 1 1 0 0 1 1 5 2 0 0 0 1 1 3 0 1 2 329
LYS 77 3A3X.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 6 2 329
LYS 82 3A3X.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 2 1 2 0 13 2 329
LYS 175 3A3X.subA 3 0 0 0 0 0 0 0 0 3 0 0 1 1 2 0 0 0 0 2 329
LYS 285 3A3X.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 3 2 329
LYS 293 3A3X.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 0 1 2 0 8 2 329
LYS 339 3A3X.subA 7 0 1 1 1 0 0 1 1 5 2 0 0 0 1 1 3 0 1 2 329
LYS 77 3A4J.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 6 2 328
LYS 82 3A4J.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 2 1 2 0 11 2 328
LYS 175 3A4J.subA 3 0 0 0 0 0 0 0 0 3 0 0 1 1 2 0 0 0 0 2 328
LYS 285 3A4J.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 6 2 328
LYS 293 3A4J.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 0 1 3 0 12 2 328
LYS 339 3A4J.subA 7 0 1 1 1 1 0 2 2 4 3 0 0 0 1 1 2 0 4 2 328
LYS 28 3B8T.subA 6 2 0 2 0 0 0 0 0 4 2 1 1 2 0 1 1 0 2 2 359
LYS 34 3B8T.subA 9 0 1 1 0 0 0 0 0 8 1 0 0 0 3 1 4 0 0 2 359
LYS 73 3B8T.subA 10 2 0 2 0 0 0 0 0 8 2 1 0 1 4 1 2 0 2 2 359
LYS 148 3B8T.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 1 0 0 0 2 2 359
LYS 167 3B8T.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 1 2 359
LYS 173 3B8T.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 1 2 359
LYS 185 3B8T.subA 8 0 2 2 0 0 0 0 0 6 2 0 0 0 3 2 1 0 2 2 359
LYS 248 3B8T.subA 4 1 1 2 0 0 1 0 1 1 3 0 0 0 1 0 0 0 1 2 359
LYS 317 3B8T.subA 5 2 0 2 0 0 0 0 0 3 2 0 1 1 1 0 1 0 1 2 359
LYS 339 3B8T.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 1 0 2 0 0 2 359
LYS 356 3B8T.subA 6 0 0 0 0 0 0 0 0 6 0 1 0 1 0 1 4 0 0 2 359
LYS 28 3B8U.subA 5 2 0 2 0 0 0 0 0 3 2 1 0 1 0 1 1 0 1 2 359
LYS 34 3B8U.subA 9 0 1 1 0 0 0 0 0 8 1 0 0 0 3 1 4 0 1 2 359
LYS 73 3B8U.subA 9 2 0 2 0 0 0 0 0 7 2 1 0 1 4 1 1 0 1 2 359
LYS 148 3B8U.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 1 0 0 0 1 2 359
LYS 167 3B8U.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 2 2 359
LYS 173 3B8U.subA 5 0 1 1 0 0 0 0 0 4 1 1 1 2 2 0 0 0 0 2 359
LYS 185 3B8U.subA 8 0 2 2 0 0 0 0 0 6 2 0 0 0 3 2 1 0 1 2 359
LYS 248 3B8U.subA 4 1 1 2 0 0 1 0 1 1 3 0 0 0 1 0 0 0 1 2 359
LYS 317 3B8U.subA 5 2 0 2 0 0 0 0 0 3 2 0 1 1 1 0 1 0 0 2 359
LYS 339 3B8U.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 1 0 2 0 0 2 359
LYS 356 3B8U.subA 6 0 0 0 0 0 0 0 0 6 0 1 0 1 1 1 3 0 1 2 359
LYS 28 3B8V.subA 6 2 0 2 0 0 0 0 0 4 2 1 1 2 0 1 1 0 2 2 359
LYS 34 3B8V.subA 9 0 1 1 0 0 0 0 0 8 1 0 0 0 3 1 4 0 1 2 359
LYS 73 3B8V.subA 9 2 0 2 0 0 0 0 0 7 2 1 0 1 4 1 1 0 4 2 359
LYS 148 3B8V.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 1 0 0 0 2 2 359
LYS 167 3B8V.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 1 2 359
LYS 173 3B8V.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 1 2 359
LYS 185 3B8V.subA 8 0 2 2 0 0 0 0 0 6 2 0 0 0 3 2 1 0 2 2 359
LYS 221 3B8V.subA 3 1 0 1 0 0 0 0 0 2 1 1 0 1 0 0 1 0 3 2 359
LYS 248 3B8V.subA 4 1 1 2 0 0 1 0 1 1 3 0 0 0 1 0 0 0 1 2 359
LYS 317 3B8V.subA 4 2 0 2 0 0 0 0 0 2 2 0 1 1 1 0 0 0 4 2 359
LYS 339 3B8V.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 1 0 2 0 0 2 359
LYS 356 3B8V.subA 7 0 0 0 0 0 0 0 0 7 0 1 1 2 0 1 4 0 0 2 359
LYS 28 3B8W.subA 6 2 0 2 0 0 0 0 0 4 2 1 1 2 0 1 1 0 3 2 359
LYS 34 3B8W.subA 9 0 1 1 0 0 0 0 0 8 1 0 0 0 3 1 4 0 1 2 359
LYS 73 3B8W.subA 9 2 0 2 0 0 0 0 0 7 2 1 0 1 4 1 1 0 2 2 359
LYS 148 3B8W.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 1 0 0 0 0 2 359
LYS 167 3B8W.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 1 2 359
LYS 173 3B8W.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 1 2 359
LYS 185 3B8W.subA 8 0 2 2 0 0 0 0 0 6 2 0 0 0 3 2 1 0 2 2 359
LYS 248 3B8W.subA 4 1 1 2 0 0 1 0 1 1 3 0 0 0 1 0 0 0 1 2 359
LYS 317 3B8W.subA 4 2 0 2 0 0 0 0 0 2 2 0 1 1 1 0 0 0 2 2 359
LYS 339 3B8W.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 1 0 2 0 0 2 359
LYS 356 3B8W.subA 8 0 0 0 0 0 0 0 0 8 0 1 1 2 1 1 4 0 1 2 359
LYS 499 3BG3.subA 8 0 0 0 0 1 0 1 1 7 1 1 2 3 2 1 1 0 0 2 680
LYS 519 3BG3.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 0 1 0 0 2 680
LYS 588 3BG3.subA 7 0 2 2 0 0 1 0 1 4 3 1 0 1 0 0 3 0 0 2 680
LYS 589 3BG3.subA 5 1 0 1 0 0 1 0 1 3 2 0 0 0 0 1 2 0 0 2 680
LYS 600 3BG3.subA 5 0 1 1 0 0 0 0 0 4 1 2 0 2 0 0 2 0 0 2 680
LYS 661 3BG3.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 0 2 3 0 0 2 680
LYS 667 3BG3.subA 10 1 2 3 0 0 0 0 0 7 3 0 0 0 5 1 1 0 0 2 680
LYS 717 3BG3.subA 7 1 0 1 0 0 1 0 1 5 2 1 0 1 2 1 1 0 0 2 680
LYS 748 3BG3.subA 8 1 0 1 0 0 0 0 0 7 1 1 0 1 3 1 2 0 0 2 680
LYS 855 3BG3.subA 5 0 0 0 0 0 0 0 0 5 0 3 0 3 1 0 1 0 0 2 680
LYS 886 3BG3.subA 9 0 1 1 0 0 0 0 0 8 1 1 1 2 3 1 2 0 0 2 680
LYS 888 3BG3.subA 7 1 1 2 2 0 0 2 2 3 4 1 1 2 0 1 0 0 0 2 680
LYS 891 3BG3.subA 9 1 0 1 2 0 0 2 2 6 3 0 2 2 0 2 2 0 0 2 680
LYS 892 3BG3.subA 6 0 3 3 2 0 0 2 2 1 5 0 0 0 1 0 0 0 0 2 680
LYS 906 3BG3.subA 12 0 1 1 0 0 0 0 0 11 1 3 2 5 1 1 4 0 0 2 680
LYS 912 3BG3.subA 7 1 0 1 0 0 0 0 0 6 1 3 0 3 1 1 1 0 0 2 680
LYS 966 3BG3.subA 6 0 2 2 0 0 0 0 0 4 2 1 0 1 1 1 1 0 0 2 680
LYS 969 3BG3.subA 7 1 2 3 0 0 0 0 0 4 3 1 0 1 1 0 2 0 0 2 680
LYS 992 3BG3.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 1 0 0 0 0 2 680
LYS 1021 3BG3.subA 9 2 1 3 0 0 1 0 1 5 4 1 0 1 1 2 1 0 0 2 680
LYS 1043 3BG3.subA 5 0 1 1 0 0 0 0 0 4 1 0 1 1 2 0 1 0 0 2 680
LYS 1056 3BG3.subA 2 0 0 0 0 0 0 0 0 2 0 1 0 1 1 0 0 0 0 2 680
LYS 1061 3BG3.subA 7 0 3 3 0 0 0 0 0 4 3 0 0 0 2 0 2 0 0 2 680
LYS 1090 3BG3.subA 6 1 0 1 0 0 0 0 0 5 1 0 1 1 2 0 2 0 0 2 680
LYS 1096 3BG3.subA 5 0 1 1 0 0 0 0 0 4 1 1 1 2 1 0 1 0 0 2 680
LYS 1103 3BG3.subA 6 1 1 2 0 0 0 0 0 4 2 1 0 1 1 0 2 0 0 2 680
LYS 1106 3BG3.subA 4 1 0 1 1 0 0 1 1 2 2 0 1 1 0 0 1 0 0 2 680
LYS 1109 3BG3.subA 5 1 0 1 1 0 0 1 1 3 2 0 1 1 1 0 1 0 0 2 680
LYS 1119 3BG3.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 1 0 3 0 0 2 680
LYS 1124 3BG3.subA 8 1 0 1 0 0 0 0 0 7 1 0 0 0 2 0 5 0 0 2 680
LYS 1130 3BG3.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 2 0 1 0 0 2 680
LYS 1133 3BG3.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 3 0 0 0 0 2 680
LYS 1144 3BG3.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 2 0 0 2 680
LYS 1159 3BG3.subA 2 0 0 0 0 1 0 1 1 1 1 0 0 0 0 0 1 0 0 2 680
LYS 1164 3BG3.subA 4 1 0 1 0 0 0 0 0 3 1 1 0 1 0 0 2 0 0 2 680
LYS 34 3C0Q.subA 9 0 2 2 0 2 0 2 2 5 4 1 0 1 2 0 2 0 1 2 277
LYS 151 3C0Q.subA 7 0 2 2 0 0 0 0 0 5 2 0 0 0 3 1 1 0 0 2 277
LYS 153 3C0Q.subA 6 0 1 1 0 2 0 2 2 3 3 0 0 0 2 1 0 0 0 2 277
LYS 188 3C0Q.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 1 0 1 0 0 2 277
LYS 238 3C0Q.subA 3 1 0 1 1 0 0 1 1 1 2 0 0 0 1 0 0 0 0 2 277
LYS 239 3C0Q.subA 4 1 0 1 1 1 0 2 2 1 3 0 0 0 1 0 0 0 0 2 277
LYS 271 3C0Q.subA 6 0 0 0 1 0 1 1 2 4 2 0 0 0 3 1 0 0 1 2 277
LYS 273 3C0Q.subA 11 0 2 2 1 0 0 1 1 8 3 0 2 2 4 1 1 0 1 2 277
LYS 34 3C0S.subA 10 0 2 2 0 2 0 2 2 6 4 1 1 2 2 0 2 0 3 2 281
LYS 151 3C0S.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 1 1 0 3 2 281
LYS 153 3C0S.subA 5 0 1 1 0 2 0 2 2 2 3 0 0 0 2 0 0 0 4 2 281
LYS 188 3C0S.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 0 2 0 7 2 281
LYS 238 3C0S.subA 1 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 2 2 281
LYS 239 3C0S.subA 5 1 0 1 0 1 0 1 1 3 2 0 1 1 2 0 0 0 4 2 281
LYS 271 3C0S.subA 7 1 0 1 0 0 1 0 1 5 2 0 0 0 4 1 0 0 8 2 281
LYS 273 3C0S.subA 10 0 2 2 0 0 0 0 0 8 2 0 2 2 4 1 1 0 5 2 281
LYS 77 3C86.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 8 2 328
LYS 82 3C86.subA 10 1 1 2 0 1 0 1 1 7 3 0 0 0 3 1 3 0 9 2 328
LYS 175 3C86.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 8 2 328
LYS 185 3C86.subA 9 0 2 2 0 1 0 1 1 6 3 1 0 1 1 0 4 0 3 2 328
LYS 285 3C86.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 2 2 328
LYS 293 3C86.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 1 1 2 0 9 2 328
LYS 339 3C86.subA 6 0 1 1 1 0 0 1 1 4 2 0 0 0 1 1 2 0 2 2 328
LYS 77 3CAK.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 12 2 330
LYS 82 3CAK.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 2 1 2 0 12 2 330
LYS 175 3CAK.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 7 2 330
LYS 185 3CAK.subA 8 0 2 2 0 1 0 1 1 5 3 0 0 0 1 0 4 0 4 2 330
LYS 285 3CAK.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 4 2 330
LYS 294 3CAK.subA 4 0 0 0 0 1 0 1 1 3 1 0 2 2 0 0 1 0 8 2 330
LYS 339 3CAK.subA 6 0 1 1 1 0 0 1 1 4 2 0 0 0 1 1 2 0 1 2 330
LYS 77 3CS2.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 5 2 331
LYS 82 3CS2.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 2 1 2 0 11 2 331
LYS 175 3CS2.subA 3 0 0 0 0 0 0 0 0 3 0 0 1 1 2 0 0 0 5 2 331
LYS 185 3CS2.subA 8 0 2 2 0 1 0 1 1 5 3 0 0 0 1 0 4 0 1 2 331
LYS 285 3CS2.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 4 2 331
LYS 294 3CS2.subA 4 0 0 0 0 1 0 1 1 3 1 0 2 2 0 0 1 0 5 2 331
LYS 339 3CS2.subA 7 0 1 1 1 0 0 1 1 5 2 0 0 0 1 1 3 0 1 2 331
LYS 6 3DC8.subA 8 1 0 1 0 0 0 0 0 7 1 1 1 2 1 0 4 0 5 2 483
LYS 18 3DC8.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 3 1 0 0 9 2 483
LYS 22 3DC8.subA 10 0 2 2 0 0 0 0 0 8 2 3 0 3 0 0 5 0 2 2 483
LYS 109 3DC8.subA 10 0 1 1 0 0 0 0 0 9 1 1 1 2 1 3 3 0 7 2 483
LYS 139 3DC8.subA 5 1 1 2 0 0 0 0 0 3 2 1 0 1 0 0 2 0 10 2 483
LYS 141 3DC8.subA 7 1 0 1 0 0 0 0 0 6 1 1 0 1 1 0 4 0 12 2 483
LYS 153 3DC8.subA 6 1 0 1 0 0 0 0 0 5 1 0 0 0 3 1 1 0 12 2 483
LYS 193 3DC8.subA 5 0 1 1 0 0 0 0 0 4 1 0 1 1 1 0 2 0 6 2 483
LYS 252 3DC8.subA 5 0 0 0 0 1 0 1 1 4 1 0 0 0 3 0 1 0 9 2 483
LYS 275 3DC8.subA 3 2 0 2 0 0 0 0 0 1 2 0 0 0 0 1 0 0 4 2 483
LYS 292 3DC8.subA 5 1 0 1 0 0 0 0 0 4 1 1 2 3 0 0 1 0 11 2 483
LYS 322 3DC8.subA 17 0 0 0 0 3 0 3 3 14 3 3 2 5 3 3 3 0 1 2 483
LYS 371 3DC8.subA 5 0 0 0 0 0 0 0 0 5 0 1 1 2 1 1 1 0 10 2 483
LYS 378 3DC8.subA 6 1 0 1 1 1 0 2 2 3 3 0 0 0 1 0 2 0 12 2 483
LYS 379 3DC8.subA 11 3 0 3 1 1 0 2 2 6 5 0 0 0 2 0 4 0 5 2 483
LYS 396 3DC8.subA 3 1 0 1 0 1 0 1 1 1 2 0 0 0 1 0 0 0 10 2 483
LYS 399 3DC8.subA 8 1 1 2 0 0 0 0 0 6 2 3 0 3 0 0 3 0 11 2 483
LYS 404 3DC8.subA 3 0 0 0 0 0 0 0 0 3 0 2 0 2 1 0 0 0 1 2 483
LYS 418 3DC8.subA 9 0 1 1 0 0 0 0 0 8 1 2 0 2 1 2 3 0 11 2 483
LYS 441 3DC8.subA 7 0 1 1 0 0 0 0 0 6 1 3 0 3 0 0 3 0 3 2 483
LYS 466 3DC8.subA 7 1 1 2 0 0 0 0 0 5 2 2 0 2 1 1 1 0 13 2 483
LYS 9 3DUG.subA 8 1 0 1 0 0 0 0 0 7 1 2 0 2 0 0 5 0 0 2 394
LYS 11 3DUG.subA 8 1 0 1 1 0 0 1 1 6 2 1 0 1 2 0 3 0 0 2 394
LYS 23 3DUG.subA 5 0 0 0 1 0 0 1 1 4 1 1 0 1 2 0 1 0 1 2 394
LYS 32 3DUG.subA 8 1 0 1 0 1 0 1 1 6 2 0 1 1 3 1 1 0 0 2 394
LYS 38 3DUG.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 0 1 0 0 2 394
LYS 42 3DUG.subA 3 1 0 1 0 0 0 0 0 2 1 1 1 2 0 0 0 0 0 2 394
LYS 71 3DUG.subA 5 1 0 1 0 0 0 0 0 4 1 1 1 2 2 0 0 0 0 2 394
LYS 93 3DUG.subA 7 0 1 1 0 0 1 0 1 5 2 1 1 2 1 1 1 0 0 2 394
LYS 168 3DUG.subA 4 0 1 1 0 1 0 1 1 2 2 0 0 0 0 1 1 0 0 2 394
LYS 172 3DUG.subA 5 0 0 0 0 1 0 1 1 4 1 0 1 1 0 2 1 0 2 2 394
LYS 175 3DUG.subA 3 0 0 0 0 2 0 2 2 1 2 0 0 0 0 1 0 0 1 2 394
LYS 199 3DUG.subA 10 0 0 0 0 0 0 0 0 10 0 1 1 2 4 1 3 0 0 2 394
LYS 207 3DUG.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 2 0 2 0 0 2 394
LYS 219 3DUG.subA 8 1 0 1 0 0 0 0 0 7 1 1 0 1 1 1 4 0 1 2 394
LYS 231 3DUG.subA 7 0 0 0 1 0 0 1 1 6 1 1 0 1 1 0 4 0 0 2 394
LYS 235 3DUG.subA 6 0 0 0 1 0 0 1 1 5 1 0 0 0 2 0 3 0 0 2 394
LYS 257 3DUG.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 0 0 2 0 0 2 394
LYS 278 3DUG.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 2 0 0 0 0 2 394
LYS 288 3DUG.subA 11 0 2 2 0 1 0 1 1 8 3 3 2 5 1 0 2 0 1 2 394
LYS 294 3DUG.subA 2 0 1 1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 0 2 394
LYS 295 3DUG.subA 7 0 1 1 0 0 0 0 0 6 1 1 2 3 0 1 2 0 0 2 394
LYS 326 3DUG.subA 6 1 1 2 0 0 0 0 0 4 2 0 2 2 1 1 0 0 1 2 394
LYS 348 3DUG.subA 7 0 0 0 0 0 0 0 0 7 0 2 1 3 1 0 3 0 0 2 394
LYS 362 3DUG.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 0 1 1 0 0 2 394
LYS 391 3DUG.subA 10 1 1 2 0 0 0 0 0 8 2 2 0 2 2 1 3 1 0 2 394
LYS 394 3DUG.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 1 1 1 0 0 2 394
LYS 397 3DUG.subA 7 0 0 0 0 1 0 1 1 6 1 0 0 0 1 2 3 0 1 2 394
LYS 77 3E3H.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 5 2 336
LYS 82 3E3H.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 2 1 2 0 9 2 336
LYS 175 3E3H.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 4 2 336
LYS 185 3E3H.subA 8 0 2 2 0 1 0 1 1 5 3 0 0 0 1 0 4 0 2 2 336
LYS 285 3E3H.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 4 2 336
LYS 294 3E3H.subA 4 0 0 0 0 1 0 1 1 3 1 0 2 2 0 0 1 0 3 2 336
LYS 339 3E3H.subA 6 0 1 1 1 0 0 1 1 4 2 0 0 0 1 1 2 0 3 2 336
LYS 8 3E74.subA 7 1 0 1 0 0 0 0 0 5 1 0 1 1 1 0 3 0 1 2 429
LYS 26 3E74.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 2 1 1 0 1 2 429
LYS 29 3E74.subA 6 1 0 1 0 0 0 0 0 5 1 0 1 1 3 0 1 0 2 2 429
LYS 41 3E74.subA 5 1 1 2 0 0 0 0 0 3 2 1 0 1 1 1 0 0 1 2 429
LYS 82 3E74.subA 10 1 1 2 0 1 0 1 1 7 3 0 1 1 3 2 1 0 5 2 429
LYS 108 3E74.subA 12 0 2 2 0 0 0 0 0 9 2 2 2 4 0 3 2 0 2 2 429
LYS 113 3E74.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 2 1 0 0 2 2 429
LYS 115 3E74.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 3 1 1 0 1 2 429
LYS 171 3E74.subA 5 0 0 0 0 0 0 0 0 5 0 0 1 1 1 2 1 0 1 2 429
LYS 175 3E74.subA 10 1 1 2 0 1 1 1 2 6 4 0 1 1 1 0 4 0 0 2 429
LYS 201 3E74.subA 6 0 1 1 0 1 0 1 1 3 2 0 0 0 2 0 1 0 1 2 429
LYS 232 3E74.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 1 1 3 0 2 2 429
LYS 286 3E74.subA 13 0 0 0 1 0 0 1 1 11 1 1 0 1 6 3 1 0 2 2 429
LYS 298 3E74.subA 5 0 1 1 0 0 0 0 0 4 1 0 2 2 1 0 1 0 0 2 429
LYS 303 3E74.subA 11 0 4 4 0 0 0 0 0 6 4 0 1 1 2 0 3 0 1 2 429
LYS 324 3E74.subA 12 0 0 0 1 0 0 1 1 10 1 1 0 1 7 1 1 0 2 2 429
LYS 330 3E74.subA 5 0 0 0 0 0 0 0 0 4 0 0 1 1 3 0 0 0 1 2 429
LYS 352 3E74.subA 4 1 1 2 0 1 0 1 1 1 3 0 1 1 0 0 0 0 4 2 429
LYS 362 3E74.subA 7 0 1 1 0 0 0 0 0 5 1 1 1 2 2 0 1 0 7 2 429
LYS 377 3E74.subA 11 3 0 3 0 1 0 1 1 7 4 1 1 2 2 0 3 0 5 2 429
LYS 384 3E74.subA 8 1 0 1 0 1 0 1 1 6 2 0 2 2 3 0 1 0 3 2 429
LYS 409 3E74.subA 8 1 1 2 0 1 1 1 2 4 4 0 1 1 0 1 2 0 1 2 429
LYS 424 3E74.subA 9 1 0 1 0 0 0 0 0 7 1 2 0 2 0 2 3 0 0 2 429
LYS 444 3E74.subA 8 1 0 1 0 1 0 1 1 6 2 0 0 0 3 0 3 0 6 2 429
LYS 450 3E74.subA 9 0 0 0 0 0 1 0 1 8 1 0 2 2 2 0 4 0 6 2 429
LYS 101 3F4C.subA 2 0 1 1 0 0 0 0 0 1 1 1 0 1 0 0 0 0 1 2 324
LYS 117 3F4C.subA 9 0 0 0 0 0 1 0 1 8 1 1 1 2 2 2 2 0 1 2 324
LYS 155 3F4C.subA 7 0 3 3 0 1 0 1 1 3 4 0 0 0 1 0 2 0 3 2 324
LYS 157 3F4C.subA 9 0 1 1 0 1 0 1 1 7 2 1 1 2 2 0 3 0 2 2 324
LYS 238 3F4C.subA 9 0 3 3 0 0 0 0 0 6 3 1 0 1 3 0 2 0 6 2 324
LYS 248 3F4C.subA 5 0 1 1 0 0 0 0 0 4 1 2 1 3 1 0 0 0 3 2 324
LYS 256 3F4C.subA 7 0 2 2 0 1 0 1 1 4 3 1 0 1 0 1 2 0 2 2 324
LYS 267 3F4C.subA 7 0 1 1 0 1 1 1 2 4 3 0 1 1 3 0 0 0 3 2 324
LYS 298 3F4C.subA 4 1 0 1 1 0 0 1 1 2 2 0 0 0 2 0 0 0 6 2 324
LYS 299 3F4C.subA 12 1 0 1 1 0 0 1 1 10 2 1 1 2 4 1 3 0 1 2 324
LYS 316 3F4C.subA 6 1 0 1 0 1 1 1 2 3 3 1 0 1 1 1 0 0 0 2 324
LYS 357 3F4C.subA 5 0 1 1 0 0 0 0 0 4 1 0 3 3 0 1 0 0 6 2 324
LYS 385 3F4C.subA 5 0 1 1 0 0 0 0 0 4 1 0 1 1 1 0 2 0 2 2 324
LYS 401 3F4C.subA 7 1 0 1 0 0 0 0 0 6 1 0 1 1 1 0 4 0 3 2 324
LYS 3 3F4D.subA 3 0 1 1 0 0 0 0 0 2 1 1 0 1 1 0 0 0 1 2 324
LYS 19 3F4D.subA 9 0 0 0 0 0 1 0 1 8 1 1 1 2 2 2 2 0 4 2 324
LYS 57 3F4D.subA 7 0 3 3 0 1 0 1 1 3 4 0 0 0 1 0 2 0 10 2 324
LYS 59 3F4D.subA 8 0 1 1 0 1 0 1 1 6 2 1 0 1 2 0 3 0 2 2 324
LYS 140 3F4D.subA 9 0 3 3 0 0 0 0 0 6 3 1 0 1 3 0 2 0 3 2 324
LYS 150 3F4D.subA 5 0 1 1 0 0 0 0 0 4 1 2 1 3 1 0 0 0 3 2 324
LYS 158 3F4D.subA 7 0 2 2 0 1 0 1 1 4 3 1 0 1 0 1 2 0 5 2 324
LYS 169 3F4D.subA 7 0 1 1 0 1 1 1 2 4 3 0 1 1 3 0 0 0 3 2 324
LYS 200 3F4D.subA 4 1 0 1 1 0 0 1 1 2 2 0 0 0 2 0 0 0 4 2 324
LYS 201 3F4D.subA 9 1 0 1 1 0 0 1 1 7 2 0 1 1 3 0 3 0 2 2 324
LYS 218 3F4D.subA 6 1 0 1 0 1 1 1 2 3 3 1 0 1 1 1 0 0 4 2 324
LYS 259 3F4D.subA 5 0 1 1 0 0 0 0 0 4 1 0 3 3 0 1 0 0 7 2 324
LYS 287 3F4D.subA 5 0 1 1 0 0 0 0 0 4 1 0 1 1 1 0 2 0 5 2 324
LYS 303 3F4D.subA 7 1 0 1 0 0 0 0 0 6 1 0 1 1 1 0 4 0 3 2 324
LYS 108 3FDK.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 0 3 0 0 7 2 322
LYS 285 3FDK.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 5 2 322
LYS 40 3FV7.subA 4 0 1 1 0 0 0 0 0 3 1 0 2 2 1 0 0 0 4 2 244
LYS 43 3FV7.subA 6 1 1 2 0 0 0 0 0 4 2 1 1 2 1 0 1 0 4 2 244
LYS 58 3FV7.subA 9 1 1 2 0 0 0 0 0 7 2 0 2 2 2 0 3 0 2 2 244
LYS 61 3FV7.subA 2 0 0 0 0 0 0 0 0 2 0 0 1 1 1 0 0 0 1 2 244
LYS 75 3FV7.subA 5 0 1 1 0 1 0 1 1 3 2 0 1 1 1 1 0 0 4 2 244
LYS 96 3FV7.subA 8 0 0 0 0 1 1 1 2 6 2 1 1 2 2 0 2 0 0 2 244
LYS 104 3FV7.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 5 2 244
LYS 108 3FV7.subA 2 1 0 1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 4 2 244
LYS 109 3FV7.subA 2 0 1 1 0 1 0 1 1 0 2 0 0 0 0 0 0 0 2 2 244
LYS 117 3FV7.subA 7 1 2 3 0 0 0 0 0 4 3 0 0 0 0 1 3 0 0 2 244
LYS 147 3FV7.subA 5 0 2 2 0 1 0 1 1 2 3 0 1 1 0 0 1 0 2 2 244
LYS 150 3FV7.subA 6 0 0 0 0 1 0 1 1 5 1 1 3 4 0 0 1 0 4 2 244
LYS 173 3FV7.subA 8 0 1 1 0 0 0 0 0 7 1 0 1 1 2 1 3 0 6 2 244
LYS 194 3FV7.subA 6 0 2 2 0 0 1 0 1 3 3 1 0 1 0 1 1 0 1 2 244
LYS 202 3FV7.subA 10 0 1 1 1 0 1 1 2 7 3 0 2 2 1 2 2 0 0 2 244
LYS 203 3FV7.subA 4 0 3 3 0 0 0 0 0 1 3 0 0 0 0 0 1 0 1 2 244
LYS 208 3FV7.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 0 1 3 0 0 2 244
LYS 214 3FV7.subA 6 0 2 2 0 0 0 0 0 4 2 1 0 1 0 1 2 0 0 2 244
LYS 218 3FV7.subA 11 0 0 0 0 0 0 0 0 11 0 5 0 5 2 0 3 0 2 2 244
LYS 242 3FV7.subA 4 0 0 0 0 0 0 0 0 4 0 0 2 2 1 0 1 0 2 2 244
LYS 243 3FV7.subA 8 0 1 1 1 0 1 1 2 5 3 0 0 0 2 1 2 0 1 2 244
LYS 253 3FV7.subA 7 0 2 2 0 0 0 0 0 5 2 1 2 3 0 0 2 0 0 2 244
LYS 267 3FV7.subA 6 0 2 2 0 0 0 0 0 4 2 1 0 1 0 2 1 0 2 2 244
LYS 40 3FYZ.subA 4 0 1 1 0 0 0 0 0 3 1 0 2 2 1 0 0 0 2 2 244
LYS 43 3FYZ.subA 6 1 1 2 0 0 0 0 0 4 2 1 1 2 1 0 1 0 3 2 244
LYS 58 3FYZ.subA 9 1 1 2 0 0 0 0 0 7 2 0 2 2 2 0 3 0 1 2 244
LYS 61 3FYZ.subA 2 0 0 0 0 0 0 0 0 2 0 0 1 1 1 0 0 0 0 2 244
LYS 75 3FYZ.subA 5 0 1 1 0 1 0 1 1 3 2 0 1 1 1 1 0 0 4 2 244
LYS 96 3FYZ.subA 7 0 0 0 0 1 1 1 2 5 2 1 1 2 2 0 1 0 2 2 244
LYS 104 3FYZ.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 2 2 244
LYS 108 3FYZ.subA 2 1 0 1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 3 2 244
LYS 109 3FYZ.subA 2 0 1 1 0 1 0 1 1 0 2 0 0 0 0 0 0 0 2 2 244
LYS 117 3FYZ.subA 7 1 2 3 0 0 0 0 0 4 3 0 0 0 0 1 3 0 0 2 244
LYS 147 3FYZ.subA 5 0 2 2 0 1 0 1 1 2 3 0 1 1 0 0 1 0 1 2 244
LYS 150 3FYZ.subA 5 0 0 0 0 1 0 1 1 4 1 1 2 3 0 0 1 0 0 2 244
LYS 173 3FYZ.subA 8 0 1 1 0 0 0 0 0 7 1 0 1 1 2 1 3 0 6 2 244
LYS 194 3FYZ.subA 5 0 2 2 0 0 0 0 0 3 2 1 0 1 0 1 1 0 1 2 244
LYS 202 3FYZ.subA 10 0 1 1 1 0 1 1 2 7 3 0 2 2 1 2 2 0 0 2 244
LYS 203 3FYZ.subA 5 0 3 3 0 0 0 0 0 2 3 0 0 0 0 0 2 0 1 2 244
LYS 208 3FYZ.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 0 1 3 0 0 2 244
LYS 214 3FYZ.subA 6 0 2 2 0 0 0 0 0 4 2 1 0 1 0 1 2 0 0 2 244
LYS 218 3FYZ.subA 11 0 0 0 0 0 0 0 0 11 0 5 0 5 2 0 3 0 2 2 244
LYS 242 3FYZ.subA 4 0 0 0 0 0 0 0 0 4 0 0 2 2 1 0 1 0 1 2 244
LYS 243 3FYZ.subA 9 0 1 1 1 0 1 1 2 6 3 0 1 1 2 1 2 0 1 2 244
LYS 253 3FYZ.subA 7 0 2 2 0 0 0 0 0 5 2 1 2 3 0 0 2 0 0 2 244
LYS 267 3FYZ.subA 6 0 2 2 0 0 0 0 0 4 2 1 0 1 0 2 1 0 1 2 244
LYS 40 3FZC.subA 4 0 1 1 0 0 0 0 0 3 1 0 2 2 1 0 0 0 3 2 244
LYS 43 3FZC.subA 6 1 1 2 0 0 0 0 0 4 2 1 1 2 1 0 1 0 4 2 244
LYS 58 3FZC.subA 9 1 1 2 0 0 0 0 0 7 2 0 2 2 2 0 3 0 2 2 244
LYS 61 3FZC.subA 2 0 0 0 0 0 0 0 0 2 0 0 1 1 1 0 0 0 1 2 244
LYS 75 3FZC.subA 5 0 1 1 0 1 0 1 1 3 2 0 1 1 1 1 0 0 3 2 244
LYS 96 3FZC.subA 8 0 0 0 0 1 1 1 2 6 2 1 1 2 2 0 2 0 2 2 244
LYS 104 3FZC.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 5 2 244
LYS 108 3FZC.subA 2 1 0 1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 3 2 244
LYS 109 3FZC.subA 2 0 1 1 0 1 0 1 1 0 2 0 0 0 0 0 0 0 2 2 244
LYS 117 3FZC.subA 7 1 2 3 0 0 0 0 0 4 3 0 0 0 0 1 3 0 0 2 244
LYS 147 3FZC.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 0 2 244
LYS 150 3FZC.subA 6 0 0 0 0 1 0 1 1 5 1 1 3 4 0 0 1 0 0 2 244
LYS 173 3FZC.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 1 1 3 0 4 2 244
LYS 194 3FZC.subA 6 0 2 2 0 0 1 0 1 3 3 1 0 1 0 1 1 0 1 2 244
LYS 202 3FZC.subA 10 0 1 1 1 0 1 1 2 7 3 0 2 2 1 2 2 0 1 2 244
LYS 203 3FZC.subA 3 0 2 2 0 0 0 0 0 1 2 0 0 0 0 0 1 0 2 2 244
LYS 208 3FZC.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 0 1 3 0 1 2 244
LYS 214 3FZC.subA 7 0 2 2 0 0 0 0 0 5 2 1 0 1 1 1 2 0 0 2 244
LYS 218 3FZC.subA 11 0 0 0 0 0 0 0 0 11 0 5 0 5 2 0 3 0 3 2 244
LYS 242 3FZC.subA 5 0 1 1 0 0 0 0 0 4 1 0 2 2 1 0 1 0 2 2 244
LYS 243 3FZC.subA 8 0 1 1 1 0 1 1 2 5 3 0 0 0 2 1 2 0 2 2 244
LYS 253 3FZC.subA 7 0 2 2 0 0 0 0 0 5 2 1 2 3 0 0 2 0 0 2 244
LYS 267 3FZC.subA 6 0 2 2 0 0 0 0 0 4 2 1 0 1 0 2 1 0 1 2 244
LYS 40 3G4P.subA 4 0 1 1 0 0 0 0 0 3 1 0 2 2 1 0 0 0 5 2 244
LYS 43 3G4P.subA 6 1 1 2 0 0 0 0 0 4 2 1 1 2 1 0 1 0 3 2 244
LYS 58 3G4P.subA 10 1 1 2 0 1 0 1 1 7 3 0 2 2 2 0 3 0 3 2 244
LYS 61 3G4P.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 1 0 0 0 6 2 244
LYS 75 3G4P.subA 5 0 1 1 0 1 0 1 1 3 2 0 1 1 1 1 0 0 3 2 244
LYS 96 3G4P.subA 8 0 0 0 0 1 1 1 2 6 2 1 1 2 2 0 2 0 3 2 244
LYS 104 3G4P.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 7 2 244
LYS 108 3G4P.subA 2 1 0 1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 4 2 244
LYS 109 3G4P.subA 2 0 1 1 0 1 0 1 1 0 2 0 0 0 0 0 0 0 4 2 244
LYS 117 3G4P.subA 7 1 2 3 0 0 0 0 0 4 3 0 0 0 0 1 3 0 1 2 244
LYS 147 3G4P.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 4 2 244
LYS 150 3G4P.subA 6 0 0 0 0 1 0 1 1 5 1 1 3 4 0 0 1 0 4 2 244
LYS 173 3G4P.subA 8 0 1 1 0 0 0 0 0 7 1 0 1 1 2 1 3 0 5 2 244
LYS 194 3G4P.subA 6 0 2 2 0 0 1 0 1 3 3 1 0 1 0 1 1 0 2 2 244
LYS 202 3G4P.subA 10 0 1 1 1 0 1 1 2 7 3 0 2 2 1 2 2 0 1 2 244
LYS 203 3G4P.subA 4 0 3 3 0 0 0 0 0 1 3 0 0 0 0 0 1 0 3 2 244
LYS 208 3G4P.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 0 1 3 0 1 2 244
LYS 214 3G4P.subA 7 0 2 2 0 0 0 0 0 5 2 1 0 1 1 1 2 0 0 2 244
LYS 218 3G4P.subA 11 0 0 0 0 0 0 0 0 11 0 5 0 5 2 0 3 0 6 2 244
LYS 242 3G4P.subA 6 0 1 1 0 0 0 0 0 5 1 0 2 2 2 0 1 0 2 2 244
LYS 243 3G4P.subA 8 0 1 1 1 0 1 1 2 5 3 0 0 0 2 1 2 0 1 2 244
LYS 253 3G4P.subA 7 0 2 2 0 0 0 0 0 5 2 1 2 3 0 0 2 0 1 2 244
LYS 267 3G4P.subA 6 0 2 2 0 0 0 0 0 4 2 1 0 1 0 2 1 0 2 2 244
LYS 242 3GTF.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 0 3 0 0 7 2 331
LYS 419 3GTF.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 6 2 331
LYS 208 3GTH.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 0 3 0 0 7 2 331
LYS 385 3GTH.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 8 2 331
LYS 208 3GTI.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 0 3 0 0 5 2 331
LYS 385 3GTI.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 5 2 331
LYS 208 3GTX.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 0 3 0 0 7 2 333
LYS 385 3GTX.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 7 2 333
LYS 208 3GU1.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 0 3 0 0 7 2 331
LYS 385 3GU1.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 9 2 331
LYS 208 3GU2.subA 7 0 1 1 0 0 0 0 0 6 1 2 0 2 1 3 0 0 3 2 330
LYS 385 3GU2.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 5 2 330
LYS 208 3GU9.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 0 3 0 0 8 2 331
LYS 385 3GU9.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 4 2 331
LYS 29 3HBR.subA 5 0 1 1 0 0 0 0 0 4 1 1 3 4 0 0 0 0 3 2 240
LYS 39 3HBR.subA 3 0 0 0 0 0 1 0 1 2 1 2 0 2 0 0 0 0 3 2 240
LYS 51 3HBR.subA 4 0 0 0 0 0 0 0 0 4 0 0 4 4 0 0 0 0 4 2 240
LYS 60 3HBR.subA 4 0 0 0 0 1 0 1 1 3 1 0 2 2 0 0 1 0 2 2 240
LYS 87 3HBR.subA 3 1 0 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 3 2 240
LYS 94 3HBR.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 5 2 240
LYS 116 3HBR.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 0 1 2 0 5 2 240
LYS 137 3HBR.subA 4 0 0 0 0 1 0 1 1 3 1 1 0 1 1 0 1 0 2 2 240
LYS 175 3HBR.subA 8 1 0 1 1 1 0 2 2 5 3 1 0 1 0 2 2 0 4 2 240
LYS 180 3HBR.subA 5 0 0 0 1 1 1 2 3 2 3 0 1 1 0 0 1 0 4 2 240
LYS 192 3HBR.subA 11 0 0 0 0 2 1 2 3 8 3 0 3 3 0 2 3 0 3 2 240
LYS 208 3HBR.subA 10 0 0 0 0 0 0 0 0 10 0 4 0 4 2 0 3 0 6 2 240
LYS 218 3HBR.subA 9 0 1 1 0 0 0 0 0 8 1 2 0 2 2 1 3 0 6 2 240
LYS 255 3HBR.subA 11 0 2 2 1 0 0 1 1 8 3 1 1 2 1 2 3 0 5 2 240
LYS 259 3HBR.subA 6 0 1 1 1 0 0 1 1 4 2 0 1 1 1 0 2 0 5 2 240
LYS 262 3HBR.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 4 2 240
LYS 208 3HTW.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 0 3 0 0 4 2 322
LYS 385 3HTW.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 1 2 322
LYS 36 3ICJ.subA 7 0 1 1 0 0 0 0 0 6 1 0 0 0 1 1 4 0 1 2 468
LYS 51 3ICJ.subA 7 0 1 1 0 0 0 0 0 6 1 2 0 2 2 0 2 0 2 2 468
LYS 52 3ICJ.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 0 1 3 0 3 2 468
LYS 87 3ICJ.subA 4 1 0 1 1 0 0 1 1 2 2 0 0 0 1 0 1 0 2 2 468
LYS 89 3ICJ.subA 7 1 0 1 1 0 0 1 1 5 2 0 0 0 1 1 3 0 2 2 468
LYS 118 3ICJ.subA 5 0 2 2 0 0 0 0 0 3 2 1 0 1 1 0 1 0 0 2 468
LYS 128 3ICJ.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 0 0 4 0 0 2 468
LYS 129 3ICJ.subA 3 0 1 1 0 1 0 1 1 1 2 0 0 0 1 0 0 0 0 2 468
LYS 176 3ICJ.subA 6 1 0 1 0 1 0 1 1 4 2 1 1 2 1 0 1 0 5 2 468
LYS 184 3ICJ.subA 3 0 0 0 0 0 0 0 0 3 0 0 1 1 1 0 1 0 0 2 468
LYS 187 3ICJ.subA 3 1 1 2 0 0 0 0 0 1 2 1 0 1 0 0 0 0 3 2 468
LYS 206 3ICJ.subA 5 0 3 3 0 1 0 1 1 1 4 0 0 0 0 0 1 0 2 2 468
LYS 211 3ICJ.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 0 0 3 0 1 2 468
LYS 216 3ICJ.subA 3 1 0 1 0 0 0 0 0 2 1 1 0 1 0 0 1 0 4 2 468
LYS 219 3ICJ.subA 9 0 3 3 0 1 1 1 2 4 5 0 0 0 1 1 2 0 4 2 468
LYS 245 3ICJ.subA 9 0 2 2 0 0 0 0 0 7 2 1 1 2 1 1 3 0 3 2 468
LYS 248 3ICJ.subA 5 0 1 1 1 0 0 1 1 3 2 0 0 0 1 0 2 0 2 2 468
LYS 260 3ICJ.subA 7 0 1 1 0 1 0 1 1 5 2 0 1 1 1 0 3 0 4 2 468
LYS 274 3ICJ.subA 5 1 1 2 1 0 0 1 1 2 3 0 0 0 0 0 2 0 6 2 468
LYS 282 3ICJ.subA 8 0 0 0 0 0 0 0 0 8 0 0 0 0 2 2 4 0 9 2 468
LYS 326 3ICJ.subA 7 2 1 3 0 0 0 0 0 4 3 0 1 1 1 0 2 0 4 2 468
LYS 337 3ICJ.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 3 1 1 0 4 2 468
LYS 351 3ICJ.subA 7 2 0 2 0 0 0 0 0 5 2 1 0 1 1 0 3 0 2 2 468
LYS 384 3ICJ.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 0 0 3 0 1 2 468
LYS 387 3ICJ.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 0 0 3 0 5 2 468
LYS 413 3ICJ.subA 6 0 2 2 0 1 0 1 1 3 3 0 0 0 1 1 1 0 10 2 468
LYS 419 3ICJ.subA 8 0 2 2 0 1 0 1 1 5 3 2 0 2 1 0 2 0 2 2 468
LYS 426 3ICJ.subA 5 0 0 0 0 1 0 1 1 4 1 1 1 2 0 0 2 0 8 2 468
LYS 483 3ICJ.subA 9 0 1 1 0 0 2 0 2 6 3 3 0 3 1 1 1 0 4 2 468
LYS 501 3ICJ.subA 2 1 0 1 0 0 0 0 0 1 1 0 0 0 0 0 1 0 0 2 468
LYS 32 3IF6.subA 4 0 0 0 0 0 0 0 0 4 0 2 0 2 0 2 0 0 2 2 228
LYS 33 3IF6.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 0 2 0 0 3 2 228
LYS 53 3IF6.subA 2 0 0 0 0 0 1 0 1 1 1 0 0 0 1 0 0 0 0 2 228
LYS 66 3IF6.subA 5 0 1 1 0 2 0 2 2 2 3 0 0 0 1 1 0 0 5 2 228
LYS 75 3IF6.subA 12 0 0 0 1 0 0 1 1 11 1 5 1 6 1 2 2 0 4 2 228
LYS 130 3IF6.subA 5 2 1 3 0 0 0 0 0 2 3 0 0 0 1 0 1 0 1 2 228
LYS 136 3IF6.subA 7 1 0 1 0 1 0 1 1 5 2 0 0 0 3 1 1 0 2 2 228
LYS 142 3IF6.subA 5 0 0 0 0 2 0 2 2 3 2 0 1 1 0 0 2 0 0 2 228
LYS 165 3IF6.subA 7 0 0 0 0 1 0 1 1 6 1 1 2 3 0 1 2 0 3 2 228
LYS 177 3IF6.subA 8 1 0 1 0 1 0 1 1 6 2 0 1 1 0 1 4 0 5 2 228
LYS 186 3IF6.subA 4 0 1 1 0 0 1 0 1 2 2 0 0 0 0 1 1 0 4 2 228
LYS 194 3IF6.subA 10 1 0 1 0 2 0 2 2 7 3 0 2 2 0 2 3 0 2 2 228
LYS 210 3IF6.subA 9 0 0 0 1 0 0 1 1 8 1 3 0 3 2 0 3 0 1 2 228
LYS 56 3ISG.subA 8 0 2 2 0 0 0 0 0 6 2 0 2 2 2 1 1 0 4 2 249
LYS 58 3ISG.subA 12 0 1 1 0 0 0 0 0 11 1 2 2 4 5 0 2 0 3 2 249
LYS 87 3ISG.subA 3 1 0 1 0 0 0 0 0 2 1 1 1 2 0 0 0 0 3 2 249
LYS 91 3ISG.subA 4 0 0 0 0 0 0 0 0 4 0 0 1 1 0 2 1 0 6 2 249
LYS 94 3ISG.subA 3 1 0 1 0 0 0 0 0 2 1 1 0 1 0 1 0 0 4 2 249
LYS 97 3ISG.subA 9 1 0 1 0 0 0 0 0 8 1 1 1 2 2 3 1 0 10 2 249
LYS 109 3ISG.subA 10 0 1 1 0 0 0 0 0 9 1 2 3 5 1 1 2 0 9 2 249
LYS 126 3ISG.subA 5 0 1 1 0 0 0 0 0 4 1 0 1 1 0 0 3 0 1 2 249
LYS 131 3ISG.subA 7 0 0 0 0 0 0 0 0 7 0 0 2 2 2 1 2 0 3 2 249
LYS 133 3ISG.subA 5 0 0 0 1 0 0 1 1 4 1 0 2 2 0 0 2 0 7 2 249
LYS 137 3ISG.subA 6 1 0 1 1 0 0 1 1 4 2 0 2 2 0 1 1 0 5 2 249
LYS 166 3ISG.subA 8 0 0 0 0 0 0 0 0 8 0 2 2 4 1 0 3 0 8 2 249
LYS 178 3ISG.subA 8 1 0 1 0 1 0 1 1 6 2 0 2 2 0 2 2 0 9 2 249
LYS 187 3ISG.subA 8 1 1 2 0 0 0 0 0 6 2 1 1 2 2 1 1 0 9 2 249
LYS 208 3ISG.subA 8 1 0 1 0 0 0 0 0 7 1 2 1 3 0 2 2 0 1 2 249
LYS 212 3ISG.subA 11 0 1 1 0 0 0 0 0 10 1 4 0 4 2 0 3 0 4 2 249
LYS 236 3ISG.subA 5 0 0 0 0 0 0 0 0 5 0 4 1 5 0 0 0 0 2 2 249
LYS 240 3ISG.subA 7 0 0 0 0 0 1 0 1 6 1 0 0 0 2 1 3 0 4 2 249
LYS 260 3ISG.subA 11 0 1 1 0 0 0 0 0 10 1 3 1 4 1 1 4 0 2 2 249
LYS 262 3ISG.subA 8 0 0 0 0 0 0 0 0 8 0 2 1 3 1 1 3 0 11 2 249
LYS 263 3ISG.subA 3 0 0 0 0 0 0 0 0 3 0 1 1 2 0 0 1 0 8 2 249
LYS 9 3JZE.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 0 1 4 0 2 2 344
LYS 27 3JZE.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 5 2 344
LYS 132 3JZE.subA 4 0 2 2 0 1 0 1 1 1 3 0 0 0 0 0 1 0 1 2 344
LYS 173 3JZE.subA 9 0 0 0 0 0 0 0 0 9 0 0 0 0 4 3 2 0 4 2 344
LYS 182 3JZE.subA 6 1 0 1 0 0 1 0 1 4 2 1 1 2 1 0 1 0 4 2 344
LYS 227 3JZE.subA 11 0 2 2 0 1 1 1 2 7 4 2 0 2 1 1 3 0 6 2 344
LYS 260 3JZE.subA 12 0 1 1 0 2 2 2 4 7 5 1 0 1 5 1 0 0 3 2 344
LYS 347 3JZE.subA 5 0 1 1 1 0 0 1 1 3 2 1 0 1 0 0 2 0 4 2 344
LYS 348 3JZE.subA 4 0 1 1 1 0 0 1 1 2 2 0 1 1 0 0 1 0 4 2 344
LYS 15 3KDN.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 1 0 0 0 2 437
LYS 21 3KDN.subA 3 0 0 0 1 0 0 1 1 2 1 1 0 1 0 1 0 0 1 2 437
LYS 22 3KDN.subA 5 0 1 1 1 1 0 2 2 2 3 1 0 1 0 1 0 0 2 2 437
LYS 73 3KDN.subA 10 1 0 1 0 0 0 0 0 9 1 1 0 1 4 3 1 0 3 2 437
LYS 116 3KDN.subA 3 0 0 0 0 1 0 1 1 2 1 0 0 0 1 0 1 0 1 2 437
LYS 119 3KDN.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 3 0 1 0 3 2 437
LYS 131 3KDN.subA 6 0 3 3 0 0 0 0 0 3 3 0 0 0 1 1 1 0 5 2 437
LYS 148 3KDN.subA 4 0 1 1 0 1 0 1 1 2 2 0 0 0 1 0 1 0 1 2 437
LYS 153 3KDN.subA 2 1 0 1 0 0 0 0 0 1 1 0 0 0 0 0 1 0 1 2 437
LYS 163 3KDN.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 4 0 2 1 0 2 437
LYS 165 3KDN.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 0 2 437
LYS 175 3KDN.subA 4 0 3 3 0 0 0 0 0 1 3 0 0 0 0 0 1 0 1 2 437
LYS 211 3KDN.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 2 2 437
LYS 215 3KDN.subA 5 1 1 2 0 0 0 0 0 3 2 0 0 0 0 1 2 0 1 2 437
LYS 223 3KDN.subA 5 0 3 3 0 0 0 0 0 2 3 0 0 0 0 0 2 0 1 2 437
LYS 224 3KDN.subA 10 1 2 3 0 0 0 0 0 7 3 2 0 2 2 0 3 0 1 2 437
LYS 250 3KDN.subA 8 0 1 1 0 0 1 0 1 6 2 1 0 1 1 0 4 0 2 2 437
LYS 303 3KDN.subA 12 1 0 1 0 0 0 0 0 11 1 0 1 1 1 4 5 0 0 2 437
LYS 322 3KDN.subA 6 0 0 0 0 0 0 0 0 6 0 0 0 0 4 0 2 0 1 2 437
LYS 343 3KDN.subA 4 0 1 1 0 0 1 0 1 2 2 0 0 0 1 1 0 0 1 2 437
LYS 355 3KDN.subA 5 0 0 0 0 0 1 0 1 4 1 1 1 2 0 2 0 0 2 2 437
LYS 360 3KDN.subA 7 0 0 0 0 0 0 0 0 7 0 1 0 1 4 1 1 0 4 2 437
LYS 426 3KDN.subA 4 1 1 2 0 0 0 0 0 2 2 1 0 1 1 0 0 0 0 2 437
LYS 429 3KDN.subA 3 0 1 1 0 1 1 1 2 0 3 0 0 0 0 0 0 0 0 2 437
LYS 437 3KDN.subA 7 0 1 1 0 1 1 1 2 4 3 0 0 0 3 1 0 0 0 2 437
LYS 15 3KDO.subA 7 1 1 2 0 0 0 0 0 5 2 0 0 0 2 2 1 0 0 2 436
LYS 21 3KDO.subA 5 0 1 1 1 0 0 1 1 3 2 1 0 1 1 1 0 0 0 2 436
LYS 22 3KDO.subA 5 0 1 1 1 1 0 2 2 2 3 1 0 1 0 1 0 0 2 2 436
LYS 73 3KDO.subA 10 1 0 1 0 0 0 0 0 9 1 1 0 1 4 3 1 0 1 2 436
LYS 116 3KDO.subA 3 0 0 0 0 1 0 1 1 2 1 0 0 0 1 0 1 0 0 2 436
LYS 119 3KDO.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 3 0 1 0 0 2 436
LYS 131 3KDO.subA 6 0 3 3 0 0 0 0 0 3 3 0 0 0 1 1 1 0 2 2 436
LYS 148 3KDO.subA 4 0 1 1 0 1 0 1 1 2 2 0 0 0 1 0 1 0 1 2 436
LYS 153 3KDO.subA 2 1 0 1 0 0 0 0 0 1 1 0 0 0 0 0 1 0 3 2 436
LYS 163 3KDO.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 4 0 2 1 0 2 436
LYS 165 3KDO.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 1 2 436
LYS 175 3KDO.subA 5 1 3 4 0 0 0 0 0 1 4 0 0 0 0 0 1 0 2 2 436
LYS 211 3KDO.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 1 0 2 0 1 2 436
LYS 215 3KDO.subA 6 1 1 2 0 0 0 0 0 4 2 0 0 0 0 1 3 0 1 2 436
LYS 223 3KDO.subA 5 0 3 3 0 0 0 0 0 2 3 0 0 0 0 0 2 0 1 2 436
LYS 224 3KDO.subA 12 1 2 3 0 1 0 1 1 8 4 2 0 2 2 1 3 0 0 2 436
LYS 250 3KDO.subA 7 0 1 1 0 0 1 0 1 5 2 1 0 1 1 0 3 0 1 2 436
LYS 303 3KDO.subA 12 1 0 1 0 0 0 0 0 11 1 0 1 1 1 4 5 0 0 2 436
LYS 322 3KDO.subA 6 0 0 0 0 0 0 0 0 6 0 0 0 0 4 0 2 0 1 2 436
LYS 343 3KDO.subA 4 0 1 1 0 0 1 0 1 2 2 0 0 0 1 1 0 0 0 2 436
LYS 355 3KDO.subA 5 0 0 0 0 0 1 0 1 4 1 1 1 2 0 2 0 0 0 2 436
LYS 360 3KDO.subA 8 0 1 1 0 0 0 0 0 7 1 1 0 1 4 1 1 0 1 2 436
LYS 426 3KDO.subA 5 1 1 2 0 0 0 0 0 3 2 1 0 1 1 1 0 0 0 2 436
LYS 429 3KDO.subA 2 0 1 1 0 0 1 0 1 0 2 0 0 0 0 0 0 0 1 2 436
LYS 437 3KDO.subA 7 0 1 1 0 1 1 1 2 4 3 0 0 0 3 1 0 0 0 2 436
LYS 4 3KZC.subA 5 0 1 1 0 0 1 0 1 3 2 0 0 0 0 1 2 0 0 2 328
LYS 28 3KZC.subA 8 0 2 2 0 1 1 1 2 4 4 0 1 1 2 1 0 0 4 2 328
LYS 31 3KZC.subA 5 0 0 0 1 0 0 1 1 4 1 0 2 2 0 1 1 0 5 2 328
LYS 37 3KZC.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 2 0 2 0 4 2 328
LYS 39 3KZC.subA 5 1 1 2 0 0 0 0 0 3 2 1 0 1 1 0 1 0 1 2 328
LYS 116 3KZC.subA 5 1 0 1 1 0 0 1 1 3 2 0 0 0 1 1 1 0 0 2 328
LYS 122 3KZC.subA 6 2 1 3 1 0 0 1 1 2 4 1 1 2 0 0 0 0 0 2 328
LYS 130 3KZC.subA 7 1 1 2 0 2 0 2 2 3 4 1 1 2 0 0 1 0 3 2 328
LYS 134 3KZC.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 1 1 1 0 1 2 328
LYS 171 3KZC.subA 9 2 0 2 0 1 0 1 1 6 3 1 0 1 1 1 3 0 0 2 328
LYS 172 3KZC.subA 6 2 0 2 0 0 0 0 0 4 2 0 0 0 3 1 0 0 0 2 328
LYS 182 3KZC.subA 3 0 0 0 0 0 1 0 1 2 1 0 0 0 2 0 0 0 0 2 328
LYS 252 3KZC.subA 12 0 0 0 0 0 1 0 1 11 1 2 1 3 3 1 3 0 1 2 328
LYS 267 3KZC.subA 5 1 2 3 0 0 0 0 0 2 3 0 0 0 1 1 0 0 0 2 328
LYS 282 3KZC.subA 9 1 0 1 0 1 1 1 2 6 3 0 0 0 0 2 4 0 0 2 328
LYS 327 3KZC.subA 11 0 0 0 1 0 1 1 2 9 2 0 2 2 2 1 4 0 2 2 328
LYS 4 3KZK.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 0 1 2 0 3 2 334
LYS 28 3KZK.subA 8 0 2 2 0 1 1 1 2 4 4 0 1 1 2 1 0 0 6 2 334
LYS 31 3KZK.subA 5 0 0 0 1 0 0 1 1 4 1 0 2 2 0 1 1 0 3 2 334
LYS 37 3KZK.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 2 0 2 0 5 2 334
LYS 39 3KZK.subA 6 1 1 2 0 0 0 0 0 4 2 1 1 2 1 0 1 0 3 2 334
LYS 74 3KZK.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 0 2 334
LYS 116 3KZK.subA 6 1 0 1 1 0 0 1 1 4 2 0 1 1 1 1 1 0 3 2 334
LYS 122 3KZK.subA 4 2 0 2 1 0 0 1 1 1 3 1 0 1 0 0 0 0 0 2 334
LYS 130 3KZK.subA 5 1 1 2 0 0 0 0 0 3 2 1 1 2 0 0 1 0 2 2 334
LYS 134 3KZK.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 1 1 1 0 1 2 334
LYS 171 3KZK.subA 7 2 0 2 0 1 0 1 1 4 3 0 0 0 1 1 2 0 0 2 334
LYS 172 3KZK.subA 11 2 0 2 0 0 0 0 0 9 2 2 0 2 4 2 1 0 1 2 334
LYS 182 3KZK.subA 4 0 0 0 0 0 1 0 1 3 1 0 0 0 2 0 1 0 0 2 334
LYS 252 3KZK.subA 12 0 0 0 0 0 1 0 1 11 1 2 1 3 3 1 3 0 3 2 334
LYS 267 3KZK.subA 5 1 2 3 0 0 0 0 0 2 3 0 0 0 1 1 0 0 0 2 334
LYS 282 3KZK.subA 10 1 0 1 0 1 1 1 2 7 3 0 0 0 0 2 5 0 1 2 334
LYS 327 3KZK.subA 11 0 0 0 1 0 1 1 2 9 2 0 2 2 2 1 4 0 2 2 334
LYS 4 3KZM.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 0 1 2 0 2 2 332
LYS 28 3KZM.subA 8 0 2 2 0 1 1 1 2 4 4 0 1 1 2 1 0 0 5 2 332
LYS 31 3KZM.subA 5 0 0 0 1 0 0 1 1 4 1 0 2 2 0 1 1 0 4 2 332
LYS 37 3KZM.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 2 0 2 0 5 2 332
LYS 39 3KZM.subA 6 1 1 2 0 0 0 0 0 4 2 1 0 1 1 0 2 0 3 2 332
LYS 74 3KZM.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 1 2 332
LYS 116 3KZM.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 1 1 0 3 2 332
LYS 122 3KZM.subA 3 2 0 2 0 0 0 0 0 1 2 1 0 1 0 0 0 0 0 2 332
LYS 130 3KZM.subA 6 1 1 2 0 1 0 1 1 3 3 1 1 2 0 0 1 0 4 2 332
LYS 134 3KZM.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 1 1 1 0 1 2 332
LYS 171 3KZM.subA 7 2 0 2 0 1 0 1 1 4 3 0 0 0 1 1 2 0 0 2 332
LYS 172 3KZM.subA 6 2 0 2 0 0 0 0 0 4 2 0 0 0 3 1 0 0 1 2 332
LYS 182 3KZM.subA 3 0 0 0 0 0 1 0 1 2 1 0 0 0 2 0 0 0 0 2 332
LYS 252 3KZM.subA 12 0 0 0 0 0 1 0 1 11 1 2 1 3 3 1 3 0 5 2 332
LYS 267 3KZM.subA 5 1 2 3 0 0 0 0 0 2 3 0 0 0 1 1 0 0 1 2 332
LYS 282 3KZM.subA 9 1 0 1 0 1 1 1 2 6 3 0 0 0 0 2 4 0 1 2 332
LYS 327 3KZM.subA 11 0 0 0 1 0 1 1 2 9 2 0 2 2 2 1 4 0 2 2 332
LYS 4 3KZN.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 0 1 2 0 3 2 332
LYS 28 3KZN.subA 8 0 2 2 0 1 1 1 2 4 4 0 1 1 2 1 0 0 8 2 332
LYS 31 3KZN.subA 5 0 0 0 1 0 0 1 1 4 1 0 2 2 0 1 1 0 5 2 332
LYS 37 3KZN.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 2 0 2 0 5 2 332
LYS 39 3KZN.subA 5 1 1 2 0 0 0 0 0 3 2 1 0 1 1 0 1 0 1 2 332
LYS 74 3KZN.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 1 2 332
LYS 116 3KZN.subA 5 1 0 1 0 0 0 0 0 4 1 0 0 0 1 2 1 0 5 2 332
LYS 122 3KZN.subA 3 2 0 2 0 0 0 0 0 1 2 1 0 1 0 0 0 0 6 2 332
LYS 130 3KZN.subA 5 1 1 2 0 0 0 0 0 3 2 1 1 2 0 0 1 0 5 2 332
LYS 134 3KZN.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 1 1 1 0 2 2 332
LYS 171 3KZN.subA 8 2 0 2 0 1 0 1 1 5 3 1 0 1 1 1 2 0 0 2 332
LYS 172 3KZN.subA 11 2 0 2 0 0 0 0 0 9 2 2 0 2 4 2 1 0 5 2 332
LYS 182 3KZN.subA 4 0 0 0 0 0 1 0 1 3 1 0 0 0 2 0 1 0 0 2 332
LYS 252 3KZN.subA 12 0 0 0 0 0 1 0 1 11 1 2 1 3 3 1 3 0 3 2 332
LYS 267 3KZN.subA 5 1 2 3 0 0 0 0 0 2 3 0 0 0 1 1 0 0 3 2 332
LYS 282 3KZN.subA 9 1 0 1 0 1 1 1 2 6 3 0 0 0 0 2 4 0 4 2 332
LYS 327 3KZN.subA 11 0 0 0 1 0 1 1 2 9 2 0 2 2 2 1 4 0 2 2 332
LYS 4 3KZO.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 0 1 2 0 3 2 332
LYS 28 3KZO.subA 9 1 2 3 0 1 1 1 2 4 5 0 1 1 2 1 0 0 6 2 332
LYS 31 3KZO.subA 5 0 0 0 1 0 0 1 1 4 1 0 2 2 0 1 1 0 5 2 332
LYS 37 3KZO.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 2 0 2 0 6 2 332
LYS 39 3KZO.subA 5 1 1 2 0 0 0 0 0 3 2 1 0 1 1 0 1 0 2 2 332
LYS 74 3KZO.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 0 2 332
LYS 116 3KZO.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 1 1 0 4 2 332
LYS 122 3KZO.subA 3 2 0 2 0 0 0 0 0 1 2 1 0 1 0 0 0 0 2 2 332
LYS 130 3KZO.subA 5 1 1 2 0 0 0 0 0 3 2 1 1 2 0 0 1 0 5 2 332
LYS 134 3KZO.subA 3 0 0 0 0 0 0 0 0 3 0 1 0 1 1 1 0 0 1 2 332
LYS 171 3KZO.subA 7 2 0 2 0 1 0 1 1 4 3 0 0 0 1 1 2 0 0 2 332
LYS 172 3KZO.subA 6 2 0 2 0 0 0 0 0 4 2 0 0 0 3 1 0 0 2 2 332
LYS 182 3KZO.subA 4 0 0 0 0 0 1 0 1 3 1 0 0 0 2 0 1 0 0 2 332
LYS 252 3KZO.subA 12 0 0 0 0 0 1 0 1 11 1 2 1 3 3 1 3 0 3 2 332
LYS 267 3KZO.subA 5 1 2 3 0 0 0 0 0 2 3 0 0 0 1 1 0 0 2 2 332
LYS 282 3KZO.subA 10 2 0 2 0 1 1 1 2 6 4 0 0 0 0 2 4 0 1 2 332
LYS 327 3KZO.subA 11 0 0 0 1 0 1 1 2 9 2 0 2 2 2 1 4 0 2 2 332
LYS 4 3L02.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 0 1 2 0 1 2 332
LYS 28 3L02.subA 9 0 2 2 0 1 1 1 2 5 4 0 1 1 2 1 1 0 3 2 332
LYS 31 3L02.subA 5 0 0 0 1 0 0 1 1 4 1 0 2 2 0 1 1 0 3 2 332
LYS 37 3L02.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 2 0 2 0 5 2 332
LYS 39 3L02.subA 5 1 1 2 0 0 0 0 0 3 2 1 0 1 1 0 1 0 0 2 332
LYS 74 3L02.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 1 2 332
LYS 116 3L02.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 1 1 0 0 2 332
LYS 122 3L02.subA 3 2 0 2 0 0 0 0 0 1 2 1 0 1 0 0 0 0 1 2 332
LYS 130 3L02.subA 5 1 1 2 0 0 0 0 0 3 2 1 1 2 0 0 1 0 7 2 332
LYS 134 3L02.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 1 1 1 0 2 2 332
LYS 171 3L02.subA 7 2 0 2 0 1 0 1 1 4 3 0 0 0 1 1 2 0 0 2 332
LYS 172 3L02.subA 6 2 0 2 0 0 0 0 0 4 2 0 0 0 3 1 0 0 1 2 332
LYS 182 3L02.subA 3 0 0 0 0 0 1 0 1 2 1 0 0 0 2 0 0 0 1 2 332
LYS 252 3L02.subA 12 0 0 0 0 0 1 0 1 11 1 2 1 3 3 1 3 0 3 2 332
LYS 267 3L02.subA 5 1 2 3 0 0 0 0 0 2 3 0 0 0 1 1 0 0 1 2 332
LYS 282 3L02.subA 9 1 0 1 0 1 1 1 2 6 3 0 0 0 0 2 4 0 1 2 332
LYS 327 3L02.subA 11 0 0 0 1 0 1 1 2 9 2 0 2 2 2 1 4 0 2 2 332
LYS 4 3L04.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 0 1 2 0 3 2 332
LYS 28 3L04.subA 8 0 2 2 0 1 1 1 2 4 4 0 1 1 2 1 0 0 2 2 332
LYS 31 3L04.subA 5 0 0 0 1 0 0 1 1 4 1 0 2 2 0 1 1 0 1 2 332
LYS 37 3L04.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 2 0 2 0 3 2 332
LYS 39 3L04.subA 7 1 1 2 0 0 0 0 0 5 2 1 0 1 2 0 2 0 1 2 332
LYS 74 3L04.subA 3 1 0 1 0 0 0 0 0 2 1 0 1 1 1 0 0 0 0 2 332
LYS 116 3L04.subA 5 1 0 1 0 0 0 0 0 4 1 0 0 0 1 2 1 0 6 2 332
LYS 122 3L04.subA 3 2 0 2 0 0 0 0 0 1 2 1 0 1 0 0 0 0 2 2 332
LYS 130 3L04.subA 5 1 1 2 0 0 0 0 0 3 2 1 1 2 0 0 1 0 3 2 332
LYS 134 3L04.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 1 1 1 0 1 2 332
LYS 171 3L04.subA 6 1 0 1 0 1 0 1 1 4 2 0 0 0 1 1 2 0 0 2 332
LYS 172 3L04.subA 6 2 0 2 0 0 0 0 0 4 2 0 0 0 3 1 0 0 2 2 332
LYS 182 3L04.subA 3 0 0 0 0 0 1 0 1 2 1 0 0 0 2 0 0 0 1 2 332
LYS 252 3L04.subA 12 0 0 0 0 0 1 0 1 11 1 2 1 3 3 1 3 0 3 2 332
LYS 267 3L04.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 2 1 0 0 3 2 332
LYS 282 3L04.subA 9 1 0 1 0 1 1 1 2 6 3 0 0 0 0 2 4 0 2 2 332
LYS 327 3L04.subA 11 0 0 0 1 0 1 1 2 9 2 0 2 2 2 1 4 0 2 2 332
LYS 4 3L05.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 0 1 2 0 1 2 332
LYS 28 3L05.subA 9 0 2 2 0 1 1 1 2 5 4 0 1 1 2 1 1 0 0 2 332
LYS 31 3L05.subA 5 0 0 0 1 0 0 1 1 4 1 0 2 2 0 1 1 0 1 2 332
LYS 37 3L05.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 2 0 2 0 2 2 332
LYS 39 3L05.subA 6 1 1 2 0 0 0 0 0 4 2 1 0 1 2 0 1 0 0 2 332
LYS 74 3L05.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 0 2 332
LYS 116 3L05.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 1 1 0 1 2 332
LYS 122 3L05.subA 3 2 0 2 0 0 0 0 0 1 2 1 0 1 0 0 0 0 1 2 332
LYS 130 3L05.subA 6 1 1 2 0 1 0 1 1 3 3 1 1 2 0 0 1 0 0 2 332
LYS 134 3L05.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 1 1 1 0 0 2 332
LYS 171 3L05.subA 6 1 0 1 0 1 0 1 1 4 2 0 0 0 1 1 2 0 0 2 332
LYS 172 3L05.subA 6 2 0 2 0 0 0 0 0 4 2 0 0 0 3 1 0 0 1 2 332
LYS 182 3L05.subA 3 0 0 0 0 0 1 0 1 2 1 0 0 0 2 0 0 0 1 2 332
LYS 252 3L05.subA 12 0 0 0 0 0 1 0 1 11 1 2 1 3 3 1 3 0 3 2 332
LYS 267 3L05.subA 5 1 2 3 0 0 0 0 0 2 3 0 0 0 1 1 0 0 2 2 332
LYS 282 3L05.subA 9 1 0 1 0 1 1 1 2 6 3 0 0 0 0 2 4 0 1 2 332
LYS 327 3L05.subA 12 0 0 0 1 0 1 1 2 10 2 0 3 3 2 1 4 0 0 2 332
LYS 4 3L06.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 0 1 2 0 2 2 332
LYS 28 3L06.subA 9 0 2 2 0 1 1 1 2 5 4 0 1 1 2 1 1 0 0 2 332
LYS 31 3L06.subA 5 0 0 0 1 0 0 1 1 4 1 0 2 2 0 1 1 0 0 2 332
LYS 37 3L06.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 2 0 2 0 2 2 332
LYS 39 3L06.subA 5 1 1 2 0 0 0 0 0 3 2 1 0 1 1 0 1 0 0 2 332
LYS 74 3L06.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 0 2 332
LYS 116 3L06.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 1 1 0 1 2 332
LYS 122 3L06.subA 3 2 0 2 0 0 0 0 0 1 2 1 0 1 0 0 0 0 1 2 332
LYS 130 3L06.subA 5 1 1 2 0 0 0 0 0 3 2 1 1 2 0 0 1 0 1 2 332
LYS 134 3L06.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 1 1 1 0 0 2 332
LYS 171 3L06.subA 6 1 0 1 0 1 0 1 1 4 2 0 0 0 1 1 2 0 0 2 332
LYS 172 3L06.subA 6 2 0 2 0 0 0 0 0 4 2 0 0 0 3 1 0 0 0 2 332
LYS 182 3L06.subA 3 0 0 0 0 0 1 0 1 2 1 0 0 0 2 0 0 0 0 2 332
LYS 252 3L06.subA 12 0 0 0 0 0 1 0 1 11 1 2 1 3 3 1 3 0 3 2 332
LYS 267 3L06.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 1 1 0 0 0 2 332
LYS 282 3L06.subA 9 1 0 1 0 1 1 1 2 6 3 0 0 0 0 2 4 0 1 2 332
LYS 327 3L06.subA 12 0 0 0 1 0 1 1 2 10 2 0 3 3 2 1 4 0 1 2 332
LYS 2 3LA4.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 2 0 3 2 837
LYS 10 3LA4.subA 12 1 4 5 0 1 0 1 1 6 6 0 0 0 2 0 4 0 4 2 837
LYS 22 3LA4.subA 5 0 0 0 0 1 0 1 1 4 1 0 1 1 0 1 2 0 3 2 837
LYS 52 3LA4.subA 6 0 1 1 0 0 0 0 0 5 1 1 1 2 0 1 2 0 4 2 837
LYS 92 3LA4.subA 7 1 0 1 0 0 0 0 0 6 1 1 0 1 2 1 2 0 0 2 837
LYS 123 3LA4.subA 8 1 0 1 0 1 0 1 1 6 2 1 0 1 3 1 1 0 3 2 837
LYS 151 3LA4.subA 6 0 0 0 0 1 0 1 1 5 1 1 0 1 2 1 1 0 0 2 837
LYS 156 3LA4.subA 5 0 1 1 0 0 0 0 0 3 1 1 0 1 0 0 2 0 7 2 837
LYS 160 3LA4.subA 5 0 0 0 0 0 0 0 0 5 0 3 0 3 1 1 0 0 3 2 837
LYS 186 3LA4.subA 8 1 2 3 0 1 0 1 1 4 4 0 0 0 2 0 2 0 3 2 837
LYS 208 3LA4.subA 8 1 0 1 0 1 0 1 1 5 2 1 0 1 0 1 3 0 1 2 837
LYS 219 3LA4.subA 7 1 1 2 0 0 0 0 0 5 2 0 1 1 1 0 3 0 3 2 837
LYS 247 3LA4.subA 8 0 0 0 0 0 0 0 0 8 0 1 1 2 3 1 2 0 6 2 837
LYS 255 3LA4.subA 3 1 2 3 0 0 0 0 0 0 3 0 0 0 0 0 0 0 1 2 837
LYS 263 3LA4.subA 5 1 1 2 0 0 0 0 0 3 2 1 1 2 0 1 0 0 1 2 837
LYS 277 3LA4.subA 4 0 1 1 0 1 1 1 2 1 3 0 1 1 0 0 0 0 2 2 837
LYS 282 3LA4.subA 7 0 1 1 0 0 0 0 0 6 1 1 1 2 1 2 1 0 2 2 837
LYS 290 3LA4.subA 9 2 1 3 0 0 0 0 0 6 3 0 1 1 2 0 3 0 3 2 837
LYS 304 3LA4.subA 5 1 1 2 0 0 0 0 0 3 2 1 0 1 0 1 1 0 4 2 837
LYS 319 3LA4.subA 8 0 0 0 0 0 0 0 0 8 0 0 0 0 4 1 3 0 7 2 837
LYS 354 3LA4.subA 12 1 0 1 0 0 0 0 0 11 1 1 1 2 4 2 3 0 4 2 837
LYS 360 3LA4.subA 10 2 0 2 0 0 0 0 0 8 2 1 0 1 3 0 4 0 5 2 837
LYS 369 3LA4.subA 13 2 1 3 0 0 0 0 0 10 3 0 1 1 5 1 3 0 9 2 837
LYS 469 3LA4.subA 7 0 0 0 0 0 0 0 0 7 0 1 0 1 4 1 0 0 1 2 837
LYS 474 3LA4.subA 4 1 0 1 0 0 0 0 0 3 1 2 0 2 1 0 0 0 5 2 837
LYS 483 3LA4.subA 6 0 1 1 0 0 1 0 1 4 2 0 1 1 2 0 1 0 6 2 837
LYS 537 3LA4.subA 6 1 0 1 0 1 0 1 1 4 2 0 0 0 3 1 0 0 2 2 837
LYS 559 3LA4.subA 9 1 1 2 0 0 0 0 0 7 2 1 0 1 0 0 6 0 5 2 837
LYS 564 3LA4.subA 5 0 0 0 0 0 0 0 0 5 0 1 1 2 0 0 3 0 2 2 837
LYS 612 3LA4.subA 4 0 1 1 0 1 0 1 1 2 2 0 0 0 0 0 2 0 2 2 837
LYS 613 3LA4.subA 5 0 0 0 0 0 0 0 0 5 0 3 0 3 0 1 1 0 1 2 837
LYS 653 3LA4.subA 13 2 1 3 0 0 0 0 0 10 3 3 2 5 2 0 3 0 2 2 837
LYS 655 3LA4.subA 11 2 0 2 0 1 0 1 1 8 3 0 1 1 4 0 3 0 2 2 837
LYS 662 3LA4.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 0 1 0 1 2 837
LYS 679 3LA4.subA 14 1 0 1 0 1 0 1 1 12 2 2 2 4 2 2 4 0 0 2 837
LYS 701 3LA4.subA 12 0 1 1 0 0 0 0 0 11 1 2 1 3 2 1 5 0 5 2 837
LYS 709 3LA4.subA 6 0 1 1 0 0 0 0 0 5 1 1 0 1 1 2 1 0 7 2 837
LYS 716 3LA4.subA 6 0 0 0 0 0 0 0 0 6 0 2 0 2 1 3 0 0 3 2 837
LYS 722 3LA4.subA 13 1 0 1 0 0 0 0 0 12 1 2 1 3 3 0 6 0 0 2 837
LYS 745 3LA4.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 1 0 3 0 5 2 837
LYS 755 3LA4.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 3 2 837
LYS 767 3LA4.subA 6 1 0 1 0 0 1 0 1 4 2 1 0 1 2 0 1 0 5 2 837
LYS 782 3LA4.subA 13 2 0 2 0 1 0 1 1 10 3 2 1 3 3 1 3 0 3 2 837
LYS 792 3LA4.subA 2 0 0 0 0 1 0 1 1 1 1 0 0 0 0 0 1 0 4 2 837
LYS 795 3LA4.subA 15 1 0 1 0 0 0 0 0 14 1 2 1 3 5 0 6 0 1 2 837
LYS 799 3LA4.subA 8 1 0 1 0 0 0 0 0 7 1 0 1 1 1 1 4 0 8 2 837
LYS 817 3LA4.subA 6 0 0 0 0 0 0 0 0 6 0 2 0 2 2 0 2 0 2 2 837
LYS 821 3LA4.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 2 0 1 0 0 2 837
LYS 30 3LCE.subA 5 0 1 1 0 0 0 0 0 4 1 2 1 3 0 1 0 0 0 2 245
LYS 45 3LCE.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 1 1 1 0 1 2 245
LYS 49 3LCE.subA 4 0 0 0 0 0 0 0 0 4 0 4 0 4 0 0 0 0 3 2 245
LYS 61 3LCE.subA 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 4 2 245
LYS 84 3LCE.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 5 2 245
LYS 91 3LCE.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 6 2 245
LYS 95 3LCE.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 5 2 245
LYS 100 3LCE.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 6 2 245
LYS 134 3LCE.subA 5 0 0 0 1 1 0 2 2 3 2 0 1 1 0 1 1 0 10 2 245
LYS 137 3LCE.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 5 2 245
LYS 138 3LCE.subA 4 0 0 0 1 0 0 1 1 3 1 1 0 1 0 2 0 0 4 2 245
LYS 152 3LCE.subA 7 1 1 2 0 0 0 0 0 5 2 0 1 1 2 1 1 0 2 2 245
LYS 177 3LCE.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 0 0 2 0 2 2 245
LYS 182 3LCE.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 0 1 0 4 2 245
LYS 189 3LCE.subA 10 0 1 1 0 0 1 0 1 8 2 0 2 2 0 2 4 0 0 2 245
LYS 205 3LCE.subA 11 0 0 0 0 0 0 0 0 11 0 5 0 5 1 0 4 0 2 2 245
LYS 228 3LCE.subA 10 0 4 4 0 0 0 0 0 6 4 0 0 0 2 2 2 0 1 2 245
LYS 246 3LCE.subA 8 0 2 2 0 0 0 0 0 6 2 1 1 2 0 0 4 0 4 2 245
LYS 251 3LCE.subA 6 0 0 0 0 1 0 1 1 5 1 3 0 3 1 1 0 0 5 2 245
LYS 256 3LCE.subA 6 0 2 2 0 0 0 0 0 4 2 2 0 2 0 0 2 0 4 2 245
LYS 4 3M4J.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 0 1 2 0 1 2 332
LYS 28 3M4J.subA 9 0 2 2 0 1 1 1 2 5 4 0 1 1 2 1 1 0 5 2 332
LYS 31 3M4J.subA 5 0 0 0 1 0 0 1 1 4 1 0 2 2 0 1 1 0 5 2 332
LYS 37 3M4J.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 2 0 2 0 5 2 332
LYS 39 3M4J.subA 6 1 1 2 0 0 0 0 0 4 2 1 0 1 1 0 2 0 3 2 332
LYS 74 3M4J.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 0 2 332
LYS 116 3M4J.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 1 1 0 1 2 332
LYS 122 3M4J.subA 3 2 0 2 0 0 0 0 0 1 2 1 0 1 0 0 0 0 0 2 332
LYS 130 3M4J.subA 5 1 1 2 0 0 0 0 0 3 2 1 1 2 0 0 1 0 1 2 332
LYS 134 3M4J.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 1 1 1 0 0 2 332
LYS 171 3M4J.subA 7 2 0 2 0 1 0 1 1 4 3 0 0 0 1 1 2 0 0 2 332
LYS 172 3M4J.subA 6 2 0 2 0 0 0 0 0 4 2 0 0 0 3 1 0 0 1 2 332
LYS 182 3M4J.subA 4 0 0 0 1 0 1 1 2 2 2 0 0 0 2 0 0 0 0 2 332
LYS 252 3M4J.subA 12 0 0 0 0 0 1 0 1 11 1 2 1 3 3 1 3 0 3 2 332
LYS 267 3M4J.subA 6 1 2 3 1 0 0 1 1 2 4 0 0 0 1 1 0 0 0 2 332
LYS 282 3M4J.subA 9 1 0 1 0 1 1 1 2 6 3 0 0 0 0 2 4 0 1 2 332
LYS 327 3M4J.subA 11 0 0 0 1 0 1 1 2 9 2 0 2 2 2 1 4 0 2 2 332
LYS 40 3MBZ.subA 4 0 1 1 0 0 0 0 0 3 1 0 2 2 1 0 0 0 3 2 244
LYS 43 3MBZ.subA 6 1 1 2 0 0 0 0 0 4 2 1 1 2 1 0 1 0 4 2 244
LYS 58 3MBZ.subA 9 1 1 2 0 0 0 0 0 7 2 0 2 2 2 0 3 0 1 2 244
LYS 61 3MBZ.subA 2 0 0 0 0 0 0 0 0 2 0 0 1 1 1 0 0 0 1 2 244
LYS 75 3MBZ.subA 5 0 1 1 0 1 0 1 1 3 2 0 1 1 1 1 0 0 3 2 244
LYS 96 3MBZ.subA 8 0 0 0 0 1 1 1 2 6 2 1 1 2 2 0 2 0 1 2 244
LYS 104 3MBZ.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 4 2 244
LYS 108 3MBZ.subA 3 1 0 1 0 1 0 1 1 1 2 0 0 0 1 0 0 0 1 2 244
LYS 109 3MBZ.subA 2 0 1 1 0 1 0 1 1 0 2 0 0 0 0 0 0 0 3 2 244
LYS 117 3MBZ.subA 6 1 1 2 0 0 0 0 0 4 2 0 0 0 0 1 3 0 0 2 244
LYS 147 3MBZ.subA 5 0 2 2 0 1 0 1 1 2 3 0 1 1 0 0 1 0 1 2 244
LYS 150 3MBZ.subA 5 0 0 0 0 1 0 1 1 4 1 1 2 3 0 0 1 0 1 2 244
LYS 173 3MBZ.subA 7 0 1 1 0 0 0 0 0 6 1 0 1 1 1 1 3 0 4 2 244
LYS 194 3MBZ.subA 5 0 2 2 0 0 0 0 0 3 2 1 0 1 0 1 1 0 4 2 244
LYS 202 3MBZ.subA 9 0 1 1 1 0 1 1 2 6 3 0 2 2 1 1 2 0 2 2 244
LYS 203 3MBZ.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 0 2 0 4 2 244
LYS 208 3MBZ.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 0 1 3 0 0 2 244
LYS 214 3MBZ.subA 7 0 2 2 0 0 0 0 0 5 2 1 0 1 1 1 2 0 0 2 244
LYS 218 3MBZ.subA 10 0 0 0 0 0 0 0 0 10 0 5 0 5 1 0 3 0 2 2 244
LYS 242 3MBZ.subA 4 0 0 0 0 0 0 0 0 4 0 0 2 2 1 0 1 0 2 2 244
LYS 243 3MBZ.subA 8 0 1 1 1 0 1 1 2 5 3 0 0 0 2 1 2 0 2 2 244
LYS 253 3MBZ.subA 7 0 2 2 0 0 0 0 0 5 2 1 2 3 0 0 2 0 0 2 244
LYS 267 3MBZ.subA 6 0 2 2 0 0 0 0 0 4 2 1 0 1 0 2 1 0 2 2 244
LYS 8 3MJM.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 0 1 4 0 3 2 344
LYS 26 3MJM.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 8 2 344
LYS 131 3MJM.subA 4 0 2 2 0 1 0 1 1 1 3 0 0 0 0 0 1 0 1 2 344
LYS 172 3MJM.subA 6 0 0 0 0 1 0 1 1 5 1 0 0 0 2 1 2 0 4 2 344
LYS 181 3MJM.subA 6 2 0 2 0 0 1 0 1 3 3 1 0 1 1 0 1 0 3 2 344
LYS 226 3MJM.subA 11 0 2 2 0 1 1 1 2 7 4 2 0 2 1 1 3 0 6 2 344
LYS 259 3MJM.subA 12 0 1 1 0 2 2 2 4 7 5 0 0 0 6 1 0 0 3 2 344
LYS 346 3MJM.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 0 0 2 0 1 2 344
LYS 38 3MKV.subA 3 1 0 1 0 0 0 0 0 2 1 1 0 1 1 0 0 0 0 2 414
LYS 41 3MKV.subA 2 0 0 0 0 0 0 0 0 2 0 1 0 1 0 0 1 0 0 2 414
LYS 51 3MKV.subA 4 1 0 1 1 0 0 1 1 2 2 0 0 0 1 0 1 0 0 2 414
LYS 53 3MKV.subA 10 1 0 1 1 0 0 1 1 8 2 1 0 1 3 0 4 0 0 2 414
LYS 114 3MKV.subA 10 1 1 2 0 0 0 0 0 8 2 0 1 1 2 3 2 0 3 2 414
LYS 287 3MKV.subA 3 0 2 2 0 0 0 0 0 1 2 0 0 0 0 1 0 0 0 2 414
LYS 297 3MKV.subA 14 2 1 3 0 0 0 0 0 11 3 3 0 3 5 0 3 0 0 2 414
LYS 313 3MKV.subA 8 0 2 2 0 1 0 1 1 5 3 0 0 0 0 0 5 0 2 2 414
LYS 318 3MKV.subA 8 0 1 1 0 0 0 0 0 7 1 0 0 0 3 1 3 0 2 2 414
LYS 366 3MKV.subA 5 1 0 1 0 0 1 0 1 3 2 0 1 1 0 0 2 0 1 2 414
LYS 387 3MKV.subA 3 0 0 0 0 0 0 0 0 3 0 1 1 2 0 0 1 0 1 2 414
LYS 404 3MKV.subA 12 2 1 3 0 2 0 2 2 7 5 1 0 1 2 1 3 0 1 2 414
LYS 30 3MTW.subA 11 0 1 1 0 0 0 0 0 10 1 1 0 1 3 0 6 0 4 2 403
LYS 44 3MTW.subA 3 0 0 0 0 0 0 0 0 3 0 1 0 1 1 1 0 0 4 2 403
LYS 63 3MTW.subA 5 2 0 2 0 0 0 0 0 3 2 0 1 1 1 0 1 0 7 2 403
LYS 64 3MTW.subA 6 1 0 1 0 0 0 0 0 5 1 1 1 2 2 0 1 0 3 2 403
LYS 121 3MTW.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 3 0 0 0 8 2 403
LYS 122 3MTW.subA 5 0 0 0 0 0 1 0 1 4 1 1 1 2 2 0 0 0 8 2 403
LYS 184 3MTW.subA 5 2 0 2 0 0 0 0 0 3 2 0 3 3 0 0 0 0 5 2 403
LYS 197 3MTW.subA 5 1 1 2 0 1 0 1 1 2 3 1 0 1 1 0 0 0 8 2 403
LYS 203 3MTW.subA 6 0 0 0 0 1 0 1 1 5 1 0 0 0 2 0 3 0 7 2 403
LYS 204 3MTW.subA 3 0 0 0 0 1 0 1 1 2 1 1 0 1 0 1 0 0 4 2 403
LYS 236 3MTW.subA 6 1 1 2 0 0 0 0 0 4 2 0 0 0 2 1 1 0 4 2 403
LYS 248 3MTW.subA 8 1 1 2 0 0 0 0 0 6 2 0 0 0 2 0 4 0 11 2 403
LYS 282 3MTW.subA 5 1 1 2 0 0 0 0 0 3 2 0 1 1 0 0 2 0 3 2 403
LYS 287 3MTW.subA 8 0 0 0 0 1 0 1 1 7 1 0 1 1 4 0 2 0 5 2 403
LYS 306 3MTW.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 1 0 2 0 4 2 403
LYS 307 3MTW.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 1 0 0 0 2 2 403
LYS 317 3MTW.subA 14 3 1 4 0 2 0 2 2 8 6 1 2 3 3 1 1 0 3 2 403
LYS 331 3MTW.subA 7 2 1 3 0 1 0 1 1 3 4 0 1 1 1 0 1 0 8 2 403
LYS 334 3MTW.subA 3 0 0 0 0 1 0 1 1 2 1 0 0 0 1 0 1 0 9 2 403
LYS 338 3MTW.subA 8 0 1 1 0 0 0 0 0 7 1 1 0 1 3 1 2 0 7 2 403
LYS 355 3MTW.subA 5 1 1 2 0 0 0 0 0 3 2 0 1 1 1 1 0 0 9 2 403
LYS 386 3MTW.subA 3 1 0 1 0 1 0 1 1 1 2 1 0 1 0 0 0 0 6 2 403
LYS 415 3MTW.subA 3 0 1 1 0 0 0 0 0 2 1 1 0 1 1 0 0 0 5 2 403
LYS 421 3MTW.subA 11 1 0 1 0 1 0 1 1 9 2 1 0 1 4 0 4 0 6 2 403
LYS 427 3MTW.subA 7 0 1 1 0 0 0 0 0 6 1 0 0 0 2 1 3 0 5 2 403
LYS 53 3N2C.subA 10 2 1 3 0 1 0 1 1 6 4 1 0 1 1 1 3 0 0 2 408
LYS 132 3N2C.subA 13 0 1 1 0 2 0 2 2 10 3 1 0 1 8 0 1 0 0 2 408
LYS 182 3N2C.subA 6 0 2 2 0 1 0 1 1 3 3 0 1 1 2 0 0 0 0 2 408
LYS 259 3N2C.subA 5 0 2 2 0 0 0 0 0 3 2 0 0 0 2 0 1 0 0 2 408
LYS 280 3N2C.subA 4 1 0 1 0 0 1 0 1 2 2 0 0 0 2 0 0 0 0 2 408
LYS 294 3N2C.subA 11 1 1 2 0 0 0 0 0 9 2 2 0 2 4 0 3 0 0 2 408
LYS 301 3N2C.subA 6 0 1 1 0 0 0 0 0 5 1 0 2 2 1 0 2 0 0 2 408
LYS 315 3N2C.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 1 1 3 0 0 2 408
LYS 407 3N2C.subA 8 1 0 1 0 1 0 1 1 6 2 2 1 3 0 0 3 0 0 2 408
LYS 159 3NWR.subA 5 2 0 2 0 1 0 1 1 2 3 0 0 0 1 0 1 0 2 2 404
LYS 169 3NWR.subA 5 1 0 1 0 0 0 0 0 4 1 0 1 1 1 0 2 0 7 2 404
LYS 324 3NWR.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 2 1 0 0 2 2 404
LYS 19 3OJG.subA 9 0 0 0 0 0 1 0 1 8 1 1 1 2 2 2 2 0 1 2 323
LYS 57 3OJG.subA 7 0 3 3 0 1 0 1 1 3 4 0 0 0 1 0 2 0 2 2 323
LYS 59 3OJG.subA 7 0 1 1 0 1 0 1 1 5 2 1 0 1 1 0 3 0 1 2 323
LYS 140 3OJG.subA 9 0 3 3 0 0 0 0 0 6 3 1 0 1 3 0 2 0 6 2 323
LYS 150 3OJG.subA 5 0 1 1 0 0 0 0 0 4 1 2 1 3 1 0 0 0 4 2 323
LYS 158 3OJG.subA 7 0 2 2 0 1 0 1 1 4 3 1 0 1 0 1 2 0 7 2 323
LYS 169 3OJG.subA 7 0 1 1 0 1 1 1 2 4 3 0 1 1 3 0 0 0 3 2 323
LYS 200 3OJG.subA 6 1 0 1 1 0 0 1 1 4 2 0 0 0 2 1 1 0 7 2 323
LYS 201 3OJG.subA 11 1 0 1 1 0 0 1 1 9 2 1 1 2 3 1 3 0 3 2 323
LYS 218 3OJG.subA 5 1 0 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 3 2 323
LYS 259 3OJG.subA 5 0 1 1 0 0 0 0 0 4 1 0 3 3 0 1 0 0 7 2 323
LYS 287 3OJG.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 1 0 1 0 1 2 323
LYS 303 3OJG.subA 7 1 0 1 0 0 0 0 0 6 1 0 1 1 1 0 4 0 2 2 323
LYS 77 3OOD.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 9 2 329
LYS 82 3OOD.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 2 1 2 0 13 2 329
LYS 175 3OOD.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 7 2 329
LYS 185 3OOD.subA 9 0 2 2 0 1 0 1 1 6 3 1 0 1 1 0 4 0 4 2 329
LYS 285 3OOD.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 2 2 329
LYS 293 3OOD.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 0 1 3 0 11 2 329
LYS 339 3OOD.subA 6 0 1 1 1 0 0 1 1 4 2 0 0 0 1 1 2 0 1 2 329
LYS 77 3OQE.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 8 2 329
LYS 82 3OQE.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 2 1 2 0 13 2 329
LYS 175 3OQE.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 8 2 329
LYS 185 3OQE.subA 9 0 2 2 0 1 0 1 1 6 3 1 0 1 1 0 4 0 7 2 329
LYS 285 3OQE.subA 8 1 0 1 1 0 0 1 1 6 2 0 0 0 2 0 4 0 2 2 329
LYS 293 3OQE.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 0 1 2 0 11 2 329
LYS 339 3OQE.subA 6 0 1 1 1 0 0 1 1 4 2 0 0 0 1 1 2 0 4 2 329
LYS 19 3ORW.subA 8 0 0 0 0 0 1 0 1 7 1 1 1 2 2 1 2 0 0 2 324
LYS 57 3ORW.subA 7 0 3 3 0 1 0 1 1 3 4 0 0 0 1 0 2 0 0 2 324
LYS 59 3ORW.subA 8 0 1 1 0 1 0 1 1 6 2 1 0 1 2 0 3 0 0 2 324
LYS 140 3ORW.subA 9 0 3 3 0 0 0 0 0 6 3 1 0 1 3 0 2 0 0 2 324
LYS 150 3ORW.subA 5 0 1 1 0 0 0 0 0 4 1 2 1 3 1 0 0 0 0 2 324
LYS 158 3ORW.subA 7 0 2 2 0 1 0 1 1 4 3 1 0 1 0 1 2 0 0 2 324
LYS 169 3ORW.subA 6 0 1 1 0 1 1 1 2 3 3 0 1 1 2 0 0 0 0 2 324
LYS 200 3ORW.subA 4 1 0 1 1 0 0 1 1 2 2 0 0 0 2 0 0 0 0 2 324
LYS 201 3ORW.subA 12 1 0 1 1 0 0 1 1 10 2 1 1 2 4 1 3 0 0 2 324
LYS 218 3ORW.subA 5 1 0 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 0 2 324
LYS 259 3ORW.subA 5 0 1 1 0 0 0 0 0 4 1 0 3 3 0 1 0 0 0 2 324
LYS 287 3ORW.subA 5 0 1 1 0 0 0 0 0 4 1 0 1 1 1 0 2 0 0 2 324
LYS 303 3ORW.subA 7 1 0 1 0 0 0 0 0 6 1 0 1 1 1 0 4 0 0 2 324
LYS 4 3OVG.subA 7 1 1 2 0 0 0 0 0 5 2 0 1 1 1 2 1 0 5 2 351
LYS 17 3OVG.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 1 0 3 0 0 2 351
LYS 29 3OVG.subA 11 1 1 2 0 1 1 1 2 6 4 0 1 1 2 1 2 1 6 2 351
LYS 51 3OVG.subA 5 1 1 2 0 0 0 0 0 3 2 0 1 1 1 0 1 0 2 2 351
LYS 54 3OVG.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 1 1 3 0 1 2 351
LYS 78 3OVG.subA 7 1 1 2 0 0 0 0 0 5 2 1 0 1 0 0 4 0 1 2 351
LYS 87 3OVG.subA 6 0 0 0 0 0 0 0 0 6 0 0 2 2 2 0 2 0 3 2 351
LYS 101 3OVG.subA 10 0 1 1 1 0 2 1 3 6 4 1 1 2 3 1 0 0 3 2 351
LYS 103 3OVG.subA 4 0 1 1 1 0 0 1 1 2 2 0 0 0 1 1 0 0 1 2 351
LYS 107 3OVG.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 3 2 351
LYS 122 3OVG.subA 4 0 2 2 0 0 0 0 0 1 2 0 0 0 0 0 1 0 2 2 351
LYS 143 3OVG.subA 5 1 1 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 7 2 351
LYS 146 3OVG.subA 2 0 0 0 0 0 0 0 0 2 0 1 0 1 1 0 0 0 3 2 351
LYS 148 3OVG.subA 8 0 3 3 0 1 0 1 1 4 4 0 0 0 3 0 1 0 3 2 351
LYS 167 3OVG.subA 15 0 1 1 0 1 0 1 1 13 2 1 0 1 3 1 8 0 2 2 351
LYS 198 3OVG.subA 7 0 3 3 0 0 1 0 1 3 4 1 0 1 1 0 1 0 4 2 351
LYS 209 3OVG.subA 14 2 0 2 0 1 0 1 1 11 3 1 1 2 4 1 4 0 3 2 351
LYS 217 3OVG.subA 8 1 0 1 0 0 1 0 1 6 2 0 2 2 1 2 1 0 4 2 351
LYS 221 3OVG.subA 5 1 0 1 0 0 0 0 0 4 1 0 1 1 0 3 0 0 1 2 351
LYS 226 3OVG.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 0 2 2 0 2 2 351
LYS 229 3OVG.subA 5 0 2 2 1 0 0 1 1 2 3 0 0 0 0 0 2 0 1 2 351
LYS 244 3OVG.subA 5 0 0 0 0 2 0 2 2 3 2 0 0 0 0 2 1 0 2 2 351
LYS 256 3OVG.subA 8 1 1 2 0 0 0 0 0 6 2 0 0 0 2 1 3 0 1 2 351
LYS 261 3OVG.subA 7 1 1 2 1 0 0 1 1 4 3 0 0 0 1 1 2 0 1 2 351
LYS 265 3OVG.subA 7 1 0 1 0 1 1 1 2 4 3 0 3 3 0 0 1 0 3 2 351
LYS 286 3OVG.subA 7 0 0 0 1 0 0 1 1 6 1 1 0 1 2 1 2 0 1 2 351
LYS 288 3OVG.subA 6 0 0 0 1 0 0 1 1 5 1 1 1 2 2 0 1 0 1 2 351
LYS 305 3OVG.subA 7 0 0 0 1 0 0 1 1 6 1 0 1 1 1 0 4 0 4 2 351
LYS 311 3OVG.subA 3 0 1 1 1 0 0 1 1 1 2 1 0 1 0 0 0 0 1 2 351
LYS 323 3OVG.subA 7 0 1 1 0 1 0 1 1 5 2 0 1 1 1 0 3 0 6 2 351
LYS 331 3OVG.subA 7 0 2 2 0 1 0 1 1 4 3 1 0 1 1 0 2 0 5 2 351
LYS 338 3OVG.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 0 0 2 0 1 2 351
LYS 341 3OVG.subA 2 0 1 1 0 0 0 0 0 1 1 1 0 1 0 0 0 0 2 2 351
LYS 347 3OVG.subA 8 0 1 1 0 0 0 0 0 7 1 0 1 1 0 1 5 0 2 2 351
LYS 348 3OVG.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 1 2 351
LYS 2 3PNU.subA 7 0 1 1 0 0 0 0 0 6 1 1 0 1 1 0 4 0 0 2 338
LYS 4 3PNU.subA 5 1 0 1 0 0 0 0 0 4 1 0 2 2 0 0 2 0 0 2 338
LYS 50 3PNU.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 1 0 3 0 0 2 338
LYS 53 3PNU.subA 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 1 4 0 0 2 338
LYS 58 3PNU.subA 4 0 0 0 0 1 0 1 1 3 1 0 0 0 1 0 2 0 0 2 338
LYS 61 3PNU.subA 4 1 1 2 0 0 0 0 0 2 2 0 0 0 2 0 0 0 0 2 338
LYS 74 3PNU.subA 5 0 0 0 0 0 0 0 0 5 0 0 1 1 1 3 0 0 3 2 338
LYS 79 3PNU.subA 3 1 1 2 0 0 0 0 0 1 2 0 0 0 0 1 0 0 2 2 338
LYS 85 3PNU.subA 6 1 0 1 0 0 0 0 0 5 1 0 0 0 1 1 3 0 2 2 338
LYS 115 3PNU.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 1 1 2 0 1 2 338
LYS 147 3PNU.subA 5 1 0 1 0 0 0 0 0 4 1 0 1 1 1 1 1 0 4 2 338
LYS 151 3PNU.subA 7 1 1 2 0 0 1 0 1 4 3 0 0 0 0 1 3 0 0 2 338
LYS 154 3PNU.subA 4 0 1 1 0 0 1 0 1 2 2 0 0 0 1 1 0 0 0 2 338
LYS 160 3PNU.subA 8 0 0 0 0 0 0 0 0 8 0 0 2 2 1 2 3 0 3 2 338
LYS 169 3PNU.subA 6 1 2 3 0 0 0 0 0 3 3 2 0 2 1 0 0 0 4 2 338
LYS 176 3PNU.subA 6 1 1 2 0 0 0 0 0 4 2 0 0 0 1 1 2 0 2 2 338
LYS 201 3PNU.subA 3 0 0 0 0 0 0 0 0 3 0 0 1 1 1 0 1 0 2 2 338
LYS 209 3PNU.subA 10 0 1 1 0 1 1 1 2 7 3 1 0 1 3 1 2 0 4 2 338
LYS 213 3PNU.subA 11 1 2 3 0 1 0 1 1 7 4 2 0 2 2 0 3 0 7 2 338
LYS 218 3PNU.subA 10 1 1 2 0 1 1 1 2 6 4 1 0 1 1 1 3 0 1 2 338
LYS 231 3PNU.subA 9 0 2 2 0 0 0 0 0 7 2 0 1 1 0 3 3 0 1 2 338
LYS 243 3PNU.subA 4 0 0 0 0 0 1 0 1 3 1 0 1 1 1 0 1 0 0 2 338
LYS 270 3PNU.subA 8 0 2 2 1 0 0 1 1 5 3 0 2 2 1 2 0 0 0 2 338
LYS 280 3PNU.subA 7 1 1 2 0 0 0 0 0 5 2 0 3 3 1 1 0 0 1 2 338
LYS 288 3PNU.subA 5 1 1 2 0 0 0 0 0 3 2 0 1 1 1 0 1 0 3 2 338
LYS 293 3PNU.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 0 1 1 0 2 2 338
LYS 295 3PNU.subA 2 0 1 1 0 0 0 0 0 1 1 0 0 0 0 1 0 0 0 2 338
LYS 298 3PNU.subA 12 2 1 3 0 1 0 1 1 8 4 0 2 2 2 2 2 0 1 2 338
LYS 305 3PNU.subA 6 0 3 3 0 0 0 0 0 3 3 1 2 3 0 0 0 0 0 2 338
LYS 316 3PNU.subA 2 1 0 1 0 0 0 0 0 1 1 0 0 0 0 1 0 0 0 2 338
LYS 330 3PNU.subA 6 0 2 2 1 0 0 1 1 3 3 0 1 1 0 1 1 0 1 2 338
LYS 334 3PNU.subA 4 0 1 1 0 0 1 0 1 2 2 1 0 1 0 0 1 0 0 2 338
LYS 45 3PNZ.subA 8 2 2 4 0 0 0 0 0 4 4 0 2 2 0 0 2 0 4 2 329
LYS 47 3PNZ.subA 5 1 1 2 0 0 1 0 1 2 3 1 0 1 0 1 0 0 3 2 329
LYS 61 3PNZ.subA 7 0 0 0 0 0 0 0 0 7 0 1 1 2 3 1 1 0 5 2 329
LYS 83 3PNZ.subA 5 0 1 1 0 0 0 0 0 4 1 1 1 2 1 1 0 0 0 2 329
LYS 97 3PNZ.subA 8 0 1 1 0 0 1 0 1 6 2 2 1 3 0 2 1 1 6 2 329
LYS 104 3PNZ.subA 5 0 2 2 0 0 0 0 0 3 2 1 0 1 1 0 1 0 1 2 329
LYS 106 3PNZ.subA 8 0 2 2 0 2 0 2 2 4 4 0 0 0 2 0 2 0 4 2 329
LYS 110 3PNZ.subA 8 1 1 2 0 0 0 0 0 6 2 0 0 0 3 1 2 0 0 2 329
LYS 130 3PNZ.subA 5 1 1 2 0 0 0 0 0 3 2 2 0 2 0 0 1 0 1 2 329
LYS 149 3PNZ.subA 9 0 2 2 0 1 0 1 1 6 3 0 1 1 3 1 1 0 3 2 329
LYS 168 3PNZ.subA 12 0 3 3 0 0 0 0 0 9 3 4 1 5 1 0 3 0 5 2 329
LYS 181 3PNZ.subA 7 0 0 0 2 0 2 2 4 3 4 1 0 1 1 1 0 0 4 2 329
LYS 202 3PNZ.subA 5 0 1 1 0 0 0 0 0 4 1 0 1 1 0 0 3 0 7 2 329
LYS 226 3PNZ.subA 8 0 1 1 1 0 1 1 2 5 3 0 1 1 1 3 0 0 6 2 329
LYS 230 3PNZ.subA 5 0 0 0 1 0 0 1 1 4 1 1 1 2 1 1 0 0 5 2 329
LYS 242 3PNZ.subA 11 2 0 2 1 1 0 2 2 7 4 0 0 0 3 2 2 0 4 2 329
LYS 244 3PNZ.subA 5 0 0 0 1 1 0 2 2 3 2 0 0 0 0 2 1 0 2 2 329
LYS 276 3PNZ.subA 9 3 0 3 0 1 0 1 1 5 4 1 0 1 1 1 2 0 4 2 329
LYS 280 3PNZ.subA 9 2 1 3 1 0 1 1 2 4 5 0 0 0 3 1 0 0 5 2 329
LYS 293 3PNZ.subA 6 0 1 1 2 0 0 2 2 3 3 0 0 0 1 1 1 0 2 2 329
LYS 294 3PNZ.subA 4 0 1 1 1 0 0 1 1 2 2 0 0 0 0 2 0 0 4 2 329
LYS 306 3PNZ.subA 6 0 2 2 0 0 0 0 0 4 2 0 0 0 2 1 1 0 4 2 329
LYS 312 3PNZ.subA 5 1 2 3 1 0 0 1 1 1 4 0 0 0 0 0 1 0 1 2 329
LYS 315 3PNZ.subA 4 1 1 2 0 0 0 0 0 2 2 0 0 0 0 0 2 0 1 2 329
LYS 316 3PNZ.subA 7 2 1 3 1 0 0 1 1 3 4 0 1 1 0 1 1 0 4 2 329
LYS 329 3PNZ.subA 11 0 1 1 1 1 1 2 3 7 4 1 0 1 3 3 0 0 5 2 329
LYS 330 3PNZ.subA 5 0 0 0 1 0 0 1 1 4 1 2 0 2 0 2 0 0 3 2 329
LYS 10 3Q7V.subB 4 0 0 0 0 0 0 0 0 4 0 1 1 2 0 1 1 0 1 2 252
LYS 11 3Q7V.subB 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 1 0 0 1 2 252
LYS 22 3Q7V.subB 9 1 1 2 0 0 0 0 0 7 2 1 0 1 0 2 4 0 2 2 252
LYS 24 3Q7V.subB 2 0 0 0 0 0 0 0 0 2 0 1 0 1 0 0 1 0 2 2 252
LYS 49 3Q7V.subB 2 0 2 2 0 0 0 0 0 0 2 0 0 0 0 0 0 0 2 2 252
LYS 53 3Q7V.subB 6 0 1 1 0 2 0 2 2 3 3 0 1 1 0 2 0 0 2 2 252
LYS 87 3Q7V.subB 5 0 0 0 0 0 2 0 2 3 2 0 1 1 0 2 0 0 3 2 252
LYS 96 3Q7V.subB 7 1 1 2 0 0 0 0 0 5 2 1 2 3 1 1 0 0 1 2 252
LYS 121 3Q7V.subB 5 0 0 0 0 0 0 0 0 5 0 0 1 1 2 1 1 0 1 2 252
LYS 129 3Q7V.subB 5 0 0 0 1 0 0 1 1 4 1 1 1 2 1 0 1 0 1 2 252
LYS 136 3Q7V.subB 7 0 0 0 1 0 0 1 1 6 1 0 3 3 1 1 1 0 3 2 252
LYS 142 3Q7V.subB 3 0 0 0 0 0 0 0 0 3 0 2 0 2 0 1 0 0 1 2 252
LYS 151 3Q7V.subB 10 2 0 2 0 1 0 1 1 7 3 2 1 3 0 2 2 0 3 2 252
LYS 162 3Q7V.subB 10 0 1 1 0 0 0 0 0 9 1 0 2 2 0 4 3 0 1 2 252
LYS 173 3Q7V.subB 3 0 0 0 1 0 0 1 1 2 1 1 1 2 0 0 0 0 2 2 252
LYS 174 3Q7V.subB 3 1 0 1 0 0 0 0 0 2 1 1 0 1 1 0 0 0 2 2 252
LYS 176 3Q7V.subB 12 0 1 1 1 0 0 1 1 10 2 1 4 5 1 1 3 0 2 2 252
LYS 186 3Q7V.subB 4 0 0 0 0 0 0 0 0 4 0 0 1 1 0 1 2 0 2 2 252
LYS 187 3Q7V.subB 4 0 1 1 0 0 0 0 0 3 1 0 1 1 0 1 1 0 1 2 252
LYS 190 3Q7V.subB 5 0 1 1 0 0 0 0 0 4 1 0 2 2 0 1 1 0 2 2 252
LYS 196 3Q7V.subB 10 0 0 0 0 0 0 0 0 10 0 4 0 4 1 0 4 0 4 2 252
LYS 205 3Q7V.subB 4 0 0 0 0 0 0 0 0 4 0 0 1 1 1 1 1 0 1 2 252
LYS 221 3Q7V.subB 8 1 1 2 0 0 0 0 0 6 2 0 0 0 0 3 3 0 2 2 252
LYS 232 3Q7V.subB 6 0 0 0 0 0 0 0 0 6 0 1 2 3 2 1 0 0 2 2 252
LYS 236 3Q7V.subB 4 0 0 0 0 0 0 0 0 4 0 1 1 2 1 0 1 0 5 2 252
LYS 244 3Q7V.subB 5 0 1 1 0 0 0 0 0 4 1 0 0 0 0 0 4 0 0 2 252
LYS 247 3Q7V.subB 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 0 2 0 1 2 252
LYS 2 3QGA.subC 2 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 2 0 0 2 564
LYS 4 3QGA.subC 4 0 1 1 0 1 0 1 1 2 2 0 1 1 0 0 1 0 0 2 564
LYS 19 3QGA.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 2 1 1 0 0 2 564
LYS 43 3QGA.subC 9 0 2 2 1 0 0 1 1 6 3 2 1 3 0 2 1 0 0 2 564
LYS 48 3QGA.subC 6 0 0 0 1 0 0 1 1 5 1 1 0 1 3 1 0 0 0 2 564
LYS 82 3QGA.subC 13 1 0 1 0 0 1 0 1 11 2 1 1 2 4 3 2 0 0 2 564
LYS 88 3QGA.subC 8 1 1 2 0 0 1 0 1 5 3 1 2 3 0 0 2 0 2 2 564
LYS 91 3QGA.subC 6 0 0 0 0 0 1 0 1 5 1 0 1 1 2 0 2 0 5 2 564
LYS 97 3QGA.subC 15 2 2 4 0 0 0 0 0 11 4 2 1 3 4 1 3 0 0 2 564
LYS 101 3QGA.subC 5 2 0 2 0 0 0 0 0 3 2 0 2 2 0 0 1 0 7 2 564
LYS 197 3QGA.subC 8 0 0 0 0 0 0 0 0 8 0 0 2 2 3 1 1 0 0 2 564
LYS 211 3QGA.subC 6 0 1 1 0 0 0 0 0 5 1 0 0 0 1 1 3 0 4 2 564
LYS 325 3QGA.subC 4 1 0 1 0 1 0 1 1 2 2 0 0 0 0 0 2 0 0 2 564
LYS 381 3QGA.subC 11 2 1 3 0 0 0 0 0 8 3 2 2 4 3 1 0 0 0 2 564
LYS 383 3QGA.subC 10 3 0 3 0 2 0 2 2 5 5 0 1 1 2 1 1 0 1 2 564
LYS 384 3QGA.subC 5 3 1 4 0 0 0 0 0 1 4 0 0 0 1 0 0 0 0 2 564
LYS 393 3QGA.subC 3 0 2 2 0 0 0 0 0 1 2 0 0 0 1 0 0 0 3 2 564
LYS 402 3QGA.subC 11 1 1 2 0 1 0 1 1 8 3 0 1 1 1 3 3 0 3 2 564
LYS 407 3QGA.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 2 2 3 0 0 2 564
LYS 429 3QGA.subC 12 0 3 3 0 0 0 0 0 9 3 2 0 2 2 1 4 0 2 2 564
LYS 444 3QGA.subC 6 0 0 0 0 0 0 0 0 6 0 1 0 1 1 3 1 0 0 2 564
LYS 446 3QGA.subC 7 0 1 1 0 0 0 0 0 6 1 0 1 1 1 1 3 0 0 2 564
LYS 450 3QGA.subC 14 1 0 1 0 0 0 0 0 13 1 3 1 4 2 0 7 0 0 2 564
LYS 483 3QGA.subC 4 0 1 1 0 0 0 0 0 3 1 0 0 0 2 1 0 0 1 2 564
LYS 485 3QGA.subC 11 0 2 2 0 0 2 0 2 7 4 2 0 2 3 2 0 0 0 2 564
LYS 495 3QGA.subC 9 1 1 2 1 0 0 1 1 6 3 1 1 2 1 1 2 0 3 2 564
LYS 503 3QGA.subC 7 0 2 2 0 1 0 1 1 4 3 0 0 0 1 0 3 0 4 2 564
LYS 505 3QGA.subC 5 0 1 1 0 0 0 0 0 4 1 0 1 1 1 0 2 0 0 2 564
LYS 511 3QGA.subC 8 1 0 1 0 1 0 1 1 6 2 2 0 2 0 1 3 0 2 2 564
LYS 516 3QGA.subC 8 1 1 2 1 0 0 1 1 5 3 0 1 1 2 1 1 0 1 2 564
LYS 523 3QGA.subC 15 0 1 1 0 0 0 0 0 14 1 3 2 5 4 1 4 0 0 2 564
LYS 524 3QGA.subC 3 1 0 1 0 0 0 0 0 2 1 1 0 1 0 0 1 0 3 2 564
LYS 527 3QGA.subC 8 1 0 1 0 0 0 0 0 7 1 1 1 2 0 3 2 0 2 2 564
LYS 534 3QGA.subC 4 0 1 1 0 0 0 0 0 3 1 0 1 1 1 0 1 0 1 2 564
LYS 549 3QGA.subC 4 0 0 0 0 0 0 0 0 4 0 0 1 1 1 0 2 0 0 2 564
LYS 554 3QGA.subC 3 0 1 1 0 0 0 0 0 2 1 1 0 1 1 0 0 0 0 2 564
LYS 2 3QGK.subC 2 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 2 0 0 2 564
LYS 4 3QGK.subC 4 0 1 1 0 1 0 1 1 2 2 0 1 1 0 0 1 0 0 2 564
LYS 19 3QGK.subC 8 2 1 3 0 0 0 0 0 5 3 0 1 1 2 1 1 0 0 2 564
LYS 43 3QGK.subC 9 0 2 2 1 0 0 1 1 6 3 2 1 3 0 2 1 0 0 2 564
LYS 48 3QGK.subC 6 0 0 0 1 0 0 1 1 5 1 1 0 1 3 1 0 0 0 2 564
LYS 82 3QGK.subC 13 1 0 1 0 0 1 0 1 11 2 1 1 2 4 3 2 0 0 2 564
LYS 88 3QGK.subC 8 1 1 2 0 0 1 0 1 5 3 1 2 3 0 0 2 0 0 2 564
LYS 91 3QGK.subC 7 0 0 0 0 0 1 0 1 6 1 0 1 1 3 0 2 0 0 2 564
LYS 97 3QGK.subC 15 2 2 4 0 0 0 0 0 11 4 2 1 3 4 1 3 0 0 2 564
LYS 101 3QGK.subC 5 2 0 2 0 0 0 0 0 3 2 0 2 2 0 0 1 0 0 2 564
LYS 197 3QGK.subC 8 0 0 0 0 0 0 0 0 8 0 0 2 2 3 1 1 0 0 2 564
LYS 211 3QGK.subC 6 0 1 1 0 0 0 0 0 5 1 0 0 0 1 1 3 0 0 2 564
LYS 325 3QGK.subC 3 1 0 1 0 1 0 1 1 1 2 0 0 0 0 0 1 0 0 2 564
LYS 381 3QGK.subC 11 2 1 3 0 0 0 0 0 8 3 2 2 4 3 1 0 0 0 2 564
LYS 383 3QGK.subC 10 3 0 3 0 2 0 2 2 5 5 0 1 1 2 1 1 0 0 2 564
LYS 384 3QGK.subC 5 3 1 4 0 0 0 0 0 1 4 0 0 0 1 0 0 0 0 2 564
LYS 393 3QGK.subC 3 0 2 2 0 0 0 0 0 1 2 0 0 0 1 0 0 0 0 2 564
LYS 402 3QGK.subC 11 1 1 2 0 1 0 1 1 8 3 0 1 1 1 3 3 0 0 2 564
LYS 407 3QGK.subC 13 1 0 1 0 1 1 1 2 10 3 3 1 4 2 2 2 0 0 2 564
LYS 429 3QGK.subC 12 0 3 3 0 0 0 0 0 9 3 2 0 2 2 1 4 0 0 2 564
LYS 444 3QGK.subC 6 0 0 0 0 0 0 0 0 6 0 1 0 1 1 3 1 0 0 2 564
LYS 446 3QGK.subC 7 0 1 1 0 0 0 0 0 6 1 0 1 1 1 1 3 0 0 2 564
LYS 450 3QGK.subC 14 1 0 1 0 0 0 0 0 13 1 3 1 4 2 0 7 0 0 2 564
LYS 483 3QGK.subC 4 0 1 1 0 0 0 0 0 3 1 0 0 0 2 1 0 0 0 2 564
LYS 485 3QGK.subC 11 0 2 2 0 0 2 0 2 7 4 2 0 2 3 2 0 0 0 2 564
LYS 495 3QGK.subC 9 1 1 2 1 0 0 1 1 6 3 1 1 2 1 1 2 0 0 2 564
LYS 503 3QGK.subC 7 0 2 2 0 1 0 1 1 4 3 0 0 0 1 0 3 0 0 2 564
LYS 505 3QGK.subC 5 0 1 1 0 0 0 0 0 4 1 0 1 1 1 0 2 0 0 2 564
LYS 511 3QGK.subC 8 1 0 1 0 1 0 1 1 6 2 2 0 2 0 1 3 0 0 2 564
LYS 516 3QGK.subC 8 1 1 2 1 0 0 1 1 5 3 0 1 1 2 1 1 0 0 2 564
LYS 523 3QGK.subC 15 0 1 1 0 0 0 0 0 14 1 3 2 5 4 1 4 0 0 2 564
LYS 524 3QGK.subC 3 1 0 1 0 0 0 0 0 2 1 1 0 1 0 0 1 0 0 2 564
LYS 527 3QGK.subC 8 1 0 1 0 0 0 0 0 7 1 1 1 2 0 3 2 0 0 2 564
LYS 534 3QGK.subC 4 0 1 1 0 0 0 0 0 3 1 0 1 1 1 0 1 0 0 2 564
LYS 549 3QGK.subC 4 0 0 0 0 0 0 0 0 4 0 0 1 1 1 0 2 0 0 2 564
LYS 554 3QGK.subC 3 0 1 1 0 0 0 0 0 2 1 1 0 1 1 0 0 0 0 2 564
LYS 30 3QNB.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 1 0 0 1 2 242
LYS 45 3QNB.subA 6 0 1 1 0 0 0 0 0 5 1 2 0 2 1 1 1 0 2 2 242
LYS 49 3QNB.subA 4 0 0 0 0 0 0 0 0 4 0 4 0 4 0 0 0 0 1 2 242
LYS 61 3QNB.subA 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 8 2 242
LYS 84 3QNB.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 1 2 242
LYS 91 3QNB.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 1 2 242
LYS 95 3QNB.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 2 1 0 0 0 2 242
LYS 100 3QNB.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 0 2 242
LYS 134 3QNB.subA 6 0 0 0 2 1 0 3 3 3 3 0 1 1 0 1 1 0 3 2 242
LYS 137 3QNB.subA 4 0 0 0 1 0 0 1 1 3 1 0 2 2 0 0 1 0 5 2 242
LYS 138 3QNB.subA 4 0 0 0 1 0 0 1 1 3 1 1 0 1 0 2 0 0 3 2 242
LYS 152 3QNB.subA 7 1 1 2 0 0 0 0 0 5 2 0 1 1 2 1 1 0 0 2 242
LYS 177 3QNB.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 0 0 2 0 3 2 242
LYS 182 3QNB.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 0 1 0 2 2 242
LYS 189 3QNB.subA 10 0 1 1 0 0 1 0 1 8 2 0 2 2 0 2 4 0 0 2 242
LYS 205 3QNB.subA 10 0 0 0 0 0 0 0 0 10 0 5 0 5 0 0 4 0 4 2 242
LYS 228 3QNB.subA 11 0 4 4 0 0 0 0 0 7 4 0 0 0 2 2 3 0 1 2 242
LYS 246 3QNB.subA 8 1 1 2 0 0 0 0 0 6 2 1 1 2 0 0 4 0 2 2 242
LYS 251 3QNB.subA 6 0 0 0 0 1 0 1 1 5 1 3 0 3 1 1 0 0 2 2 242
LYS 256 3QNB.subA 6 0 2 2 0 0 0 0 0 4 2 2 0 2 0 0 2 0 2 2 242
LYS 30 3QNC.subA 5 0 1 1 0 0 0 0 0 4 1 2 1 3 0 1 0 0 0 2 242
LYS 45 3QNC.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 1 1 1 0 6 2 242
LYS 49 3QNC.subA 3 0 0 0 0 0 0 0 0 3 0 3 0 3 0 0 0 0 2 2 242
LYS 61 3QNC.subA 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 9 2 242
LYS 84 3QNC.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 10 2 242
LYS 91 3QNC.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 0 2 1 0 8 2 242
LYS 95 3QNC.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 3 2 242
LYS 100 3QNC.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 0 0 1 0 2 2 242
LYS 134 3QNC.subA 5 0 0 0 1 1 0 2 2 3 2 0 1 1 0 1 1 0 6 2 242
LYS 137 3QNC.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 0 1 0 6 2 242
LYS 138 3QNC.subA 4 0 0 0 1 0 0 1 1 3 1 1 0 1 0 2 0 0 6 2 242
LYS 152 3QNC.subA 8 1 1 2 0 1 0 1 1 5 3 0 1 1 2 1 1 0 3 2 242
LYS 177 3QNC.subA 4 0 0 0 0 0 0 0 0 4 0 1 1 2 0 0 2 0 6 2 242
LYS 182 3QNC.subA 4 0 1 1 0 0 0 0 0 3 1 1 1 2 0 0 1 0 5 2 242
LYS 189 3QNC.subA 10 0 1 1 0 0 1 0 1 8 2 0 2 2 0 2 4 0 0 2 242
LYS 205 3QNC.subA 11 0 0 0 0 0 0 0 0 11 0 5 0 5 1 0 4 0 4 2 242
LYS 218 3QNC.subA 10 1 1 2 0 0 0 0 0 8 2 2 1 3 1 1 3 0 2 2 242
LYS 228 3QNC.subA 9 0 4 4 0 0 0 0 0 5 4 0 0 0 1 2 2 0 1 2 242
LYS 246 3QNC.subA 8 1 1 2 0 0 0 0 0 6 2 1 1 2 0 0 4 0 2 2 242
LYS 251 3QNC.subA 6 0 0 0 0 1 0 1 1 5 1 3 0 3 1 1 0 0 10 2 242
LYS 256 3QNC.subA 6 0 2 2 0 0 0 0 0 4 2 2 0 2 0 0 2 0 7 2 242
LYS 2 3S46.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 0 1 0 3 2 367
LYS 10 3S46.subA 9 0 2 2 0 0 0 0 0 7 2 1 0 1 3 1 2 0 4 2 367
LYS 35 3S46.subA 12 2 0 2 0 0 0 0 0 10 2 1 0 1 1 0 8 0 2 2 367
LYS 53 3S46.subA 6 0 0 0 0 0 0 0 0 6 0 0 0 0 4 0 2 0 6 2 367
LYS 80 3S46.subA 9 3 0 3 0 0 0 0 0 6 3 1 1 2 2 1 1 0 4 2 367
LYS 97 3S46.subA 9 1 1 2 0 0 0 0 0 7 2 1 0 1 3 0 3 0 0 2 367
LYS 117 3S46.subA 7 1 2 3 0 0 0 0 0 4 3 0 0 0 1 0 3 0 1 2 367
LYS 185 3S46.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 0 2 1 0 0 2 367
LYS 192 3S46.subA 3 0 2 2 0 0 0 0 0 1 2 0 0 0 0 0 1 0 0 2 367
LYS 253 3S46.subA 8 0 0 0 0 0 1 0 1 7 1 3 0 3 0 0 4 0 3 2 367
LYS 319 3S46.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 1 1 0 5 2 367
LYS 326 3S46.subA 6 0 1 1 0 0 0 0 0 5 1 3 1 4 0 0 1 0 6 2 367
LYS 336 3S46.subA 7 1 1 2 0 0 0 0 0 5 2 0 2 2 2 0 1 0 3 2 367
LYS 77 3SO7.subA 4 0 1 1 0 1 0 1 1 2 2 1 0 1 1 0 0 0 4 2 329
LYS 82 3SO7.subA 8 1 1 2 0 1 0 1 1 5 3 0 0 0 2 1 2 0 10 2 329
LYS 175 3SO7.subA 4 0 1 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 4 2 329
LYS 185 3SO7.subA 8 0 2 2 0 1 0 1 1 5 3 0 0 0 1 0 4 0 3 2 329
LYS 285 3SO7.subA 9 1 0 1 1 0 0 1 1 7 2 0 0 0 3 0 4 0 3 2 329
LYS 293 3SO7.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 0 1 3 0 8 2 329
LYS 339 3SO7.subA 6 0 1 1 1 0 0 1 1 4 2 0 0 0 1 1 2 0 3 2 329
LYS 5 3TW6.subA 8 1 0 1 1 0 0 1 1 6 2 3 0 3 2 0 1 0 0 2 1007
LYS 28 3TW6.subA 5 1 0 1 1 0 0 1 1 3 2 2 0 2 0 0 1 0 1 2 1007
LYS 38 3TW6.subA 7 1 2 3 0 1 0 1 1 3 4 0 0 0 1 1 1 0 0 2 1007
LYS 45 3TW6.subA 7 0 0 0 0 0 1 0 1 6 1 0 0 0 2 2 2 0 0 2 1007
LYS 79 3TW6.subA 4 0 0 0 0 1 0 1 1 3 1 0 0 0 1 0 2 0 0 2 1007
LYS 105 3TW6.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 2 0 0 0 0 2 1007
LYS 114 3TW6.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 1 2 1007
LYS 124 3TW6.subA 2 0 0 0 0 0 0 0 0 2 0 0 1 1 0 0 1 0 0 2 1007
LYS 209 3TW6.subA 2 0 0 0 0 0 0 0 0 2 0 0 0 0 1 0 1 0 0 2 1007
LYS 245 3TW6.subA 9 0 3 3 0 0 1 0 1 5 4 0 2 2 1 0 2 0 0 2 1007
LYS 269 3TW6.subA 6 0 0 0 0 0 0 0 0 6 0 0 0 0 2 1 3 0 0 2 1007
LYS 292 3TW6.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 1 2 0 0 0 2 1007
LYS 319 3TW6.subA 7 1 1 2 0 0 1 0 1 4 3 0 0 0 2 0 2 0 1 2 1007
LYS 404 3TW6.subA 12 1 1 2 0 1 1 1 2 8 4 1 1 2 1 1 4 0 0 2 1007
LYS 446 3TW6.subA 7 1 0 1 0 0 1 0 1 5 2 1 0 1 1 3 0 0 0 2 1007
LYS 468 3TW6.subA 2 0 0 0 0 1 0 1 1 1 1 0 0 0 0 0 1 0 0 2 1007
LYS 475 3TW6.subA 10 2 0 2 0 1 0 1 1 7 3 3 1 4 1 1 1 0 0 2 1007
LYS 492 3TW6.subA 6 1 0 1 0 0 1 0 1 4 2 0 1 1 3 0 0 0 2 2 1007
LYS 496 3TW6.subA 4 0 0 0 0 1 0 1 1 3 1 0 0 0 2 1 0 0 1 2 1007
LYS 515 3TW6.subA 4 1 0 1 0 0 0 0 0 3 1 0 1 1 1 0 1 0 0 2 1007
LYS 519 3TW6.subA 10 2 0 2 0 1 0 1 1 7 3 1 2 3 2 0 2 0 0 2 1007
LYS 528 3TW6.subA 4 0 1 1 1 0 1 1 2 1 3 0 0 0 1 0 0 0 3 2 1007
LYS 529 3TW6.subA 5 0 1 1 1 0 0 1 1 3 2 0 0 0 1 1 1 0 2 2 1007
LYS 538 3TW6.subA 3 0 2 2 0 1 0 1 1 0 3 0 0 0 0 0 0 0 2 2 1007
LYS 637 3TW6.subA 6 1 0 1 0 0 0 0 0 5 1 0 1 1 1 1 2 0 3 2 1007
LYS 645 3TW6.subA 6 1 0 1 0 1 0 1 1 4 2 0 1 1 2 1 0 0 0 2 1007
LYS 675 3TW6.subA 8 1 1 2 0 0 0 0 0 6 2 0 1 1 1 1 3 0 1 2 1007
LYS 694 3TW6.subA 8 1 0 1 0 1 0 1 1 6 2 1 0 1 3 1 1 0 1 2 1007
LYS 698 3TW6.subA 5 1 0 1 0 0 0 0 0 4 1 1 1 2 0 1 1 0 3 2 1007
LYS 709 3TW6.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 1 0 1 0 0 2 1007
LYS 725 3TW6.subA 8 1 0 1 0 0 0 0 0 7 1 0 0 0 3 1 3 0 5 2 1007
LYS 730 3TW6.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 4 0 1 0 4 2 1007
LYS 734 3TW6.subA 4 0 1 1 0 0 0 0 0 3 1 0 0 0 1 1 1 0 4 2 1007
LYS 829 3TW6.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 0 2 1007
LYS 849 3TW6.subA 4 0 1 1 0 0 0 0 0 3 1 1 0 1 0 1 1 0 0 2 1007
LYS 880 3TW6.subA 14 1 1 2 0 1 1 1 2 10 4 3 1 4 1 1 4 0 3 2 1007
LYS 886 3TW6.subA 6 1 0 1 0 0 0 0 0 5 1 2 1 3 1 0 1 0 1 2 1007
LYS 923 3TW6.subA 10 1 1 2 0 0 0 0 0 8 2 1 0 1 2 1 4 0 4 2 1007
LYS 939 3TW6.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 1 0 0 0 0 2 1007
LYS 940 3TW6.subA 12 1 1 2 0 1 0 1 1 9 3 1 0 1 3 1 4 0 0 2 1007
LYS 943 3TW6.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 0 1 0 0 2 1007
LYS 946 3TW6.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 2 2 1007
LYS 957 3TW6.subA 2 0 1 1 0 0 0 0 0 1 1 0 0 0 0 0 1 0 2 2 1007
LYS 966 3TW6.subA 5 1 1 2 0 1 0 1 1 2 3 0 0 0 1 0 1 0 1 2 1007
LYS 970 3TW6.subA 2 0 1 1 0 0 0 0 0 1 1 0 0 0 0 0 1 0 0 2 1007
LYS 971 3TW6.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 0 0 3 0 0 2 1007
LYS 989 3TW6.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 0 2 1007
LYS 1028 3TW6.subA 2 0 1 1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 0 2 1007
LYS 1030 3TW6.subA 4 0 1 1 0 0 0 0 0 3 1 1 0 1 1 0 1 0 0 2 1007
LYS 1062 3TW6.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 0 0 4 0 0 2 1007
LYS 5 3TW7.subA 10 1 1 2 1 0 0 1 1 7 3 3 0 3 2 0 2 0 0 2 1004
LYS 28 3TW7.subA 6 1 0 1 1 0 0 1 1 4 2 2 0 2 1 0 1 0 0 2 1004
LYS 38 3TW7.subA 6 1 1 2 0 1 0 1 1 3 3 0 0 0 1 1 1 0 0 2 1004
LYS 45 3TW7.subA 8 0 0 0 0 1 1 1 2 6 2 0 0 0 2 2 2 0 0 2 1004
LYS 79 3TW7.subA 5 0 0 0 0 1 0 1 1 4 1 0 0 0 2 0 2 0 0 2 1004
LYS 105 3TW7.subA 5 1 0 1 0 0 0 0 0 4 1 0 1 1 3 0 0 0 0 2 1004
LYS 114 3TW7.subA 5 1 0 1 0 0 0 0 0 4 1 1 0 1 2 0 1 0 0 2 1004
LYS 124 3TW7.subA 9 0 0 0 0 1 0 1 1 8 1 1 1 2 3 1 2 0 0 2 1004
LYS 166 3TW7.subA 7 0 1 1 0 0 0 0 0 6 1 0 0 0 1 1 4 0 0 2 1004
LYS 209 3TW7.subA 6 1 1 2 0 0 0 0 0 4 2 0 0 0 2 0 2 0 0 2 1004
LYS 245 3TW7.subA 10 0 4 4 0 0 1 0 1 5 5 0 2 2 1 0 2 0 0 2 1004
LYS 269 3TW7.subA 6 0 0 0 0 0 0 0 0 6 0 0 0 0 1 1 4 0 0 2 1004
LYS 292 3TW7.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 1 2 0 0 0 2 1004
LYS 319 3TW7.subA 8 1 1 2 0 0 1 0 1 5 3 0 0 0 2 0 3 0 0 2 1004
LYS 404 3TW7.subA 13 1 1 2 0 1 1 1 2 9 4 2 1 3 1 0 5 0 0 2 1004
LYS 446 3TW7.subA 7 0 0 0 0 0 1 0 1 6 1 2 0 2 1 3 0 0 0 2 1004
LYS 475 3TW7.subA 8 1 0 1 0 0 0 0 0 7 1 3 0 3 2 1 1 0 0 2 1004
LYS 492 3TW7.subA 6 1 0 1 0 0 1 0 1 4 2 0 1 1 3 0 0 0 0 2 1004
LYS 496 3TW7.subA 4 0 0 0 0 1 0 1 1 3 1 0 0 0 2 1 0 0 1 2 1004
LYS 515 3TW7.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 0 2 1004
LYS 519 3TW7.subA 10 2 0 2 0 1 0 1 1 7 3 1 2 3 2 0 2 0 0 2 1004
LYS 528 3TW7.subA 4 0 1 1 1 0 1 1 2 1 3 0 0 0 1 0 0 0 0 2 1004
LYS 529 3TW7.subA 5 0 1 1 1 0 0 1 1 3 2 0 0 0 1 1 1 0 0 2 1004
LYS 538 3TW7.subA 3 0 2 2 0 1 0 1 1 0 3 0 0 0 0 0 0 0 0 2 1004
LYS 637 3TW7.subA 8 2 1 3 0 0 0 0 0 5 3 0 1 1 1 1 2 0 0 2 1004
LYS 645 3TW7.subA 5 0 0 0 0 1 0 1 1 4 1 0 1 1 2 1 0 0 0 2 1004
LYS 675 3TW7.subA 8 1 1 2 0 0 0 0 0 6 2 0 1 1 1 1 3 0 0 2 1004
LYS 694 3TW7.subA 7 1 0 1 0 1 0 1 1 5 2 1 0 1 2 1 1 0 0 2 1004
LYS 698 3TW7.subA 5 1 0 1 0 0 0 0 0 4 1 1 1 2 0 1 1 0 0 2 1004
LYS 709 3TW7.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 1 0 1 0 0 2 1004
LYS 725 3TW7.subA 7 1 0 1 0 0 0 0 0 6 1 0 0 0 3 1 2 0 1 2 1004
LYS 730 3TW7.subA 7 0 0 0 0 0 0 0 0 7 0 0 0 0 6 0 1 0 0 2 1004
LYS 734 3TW7.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 2 1 1 0 0 2 1004
LYS 829 3TW7.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 0 2 1004
LYS 849 3TW7.subA 5 0 1 1 0 0 0 0 0 4 1 1 0 1 0 2 1 0 0 2 1004
LYS 880 3TW7.subA 12 1 1 2 0 0 1 0 1 9 3 3 1 4 0 1 4 0 0 2 1004
LYS 886 3TW7.subA 5 1 0 1 0 0 0 0 0 4 1 1 1 2 1 0 1 0 0 2 1004
LYS 923 3TW7.subA 9 1 1 2 0 0 0 0 0 7 2 1 0 1 2 1 3 0 0 2 1004
LYS 939 3TW7.subA 3 0 1 1 0 0 0 0 0 2 1 0 1 1 1 0 0 0 0 2 1004
LYS 940 3TW7.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 0 1 0 0 2 1004
LYS 943 3TW7.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 0 1 0 0 2 1004
LYS 946 3TW7.subA 3 0 1 1 0 0 0 0 0 2 1 0 0 0 2 0 0 0 0 2 1004
LYS 957 3TW7.subA 2 0 1 1 0 0 0 0 0 1 1 0 0 0 0 0 1 0 0 2 1004
LYS 966 3TW7.subA 4 1 0 1 0 1 0 1 1 2 2 0 0 0 1 0 1 0 0 2 1004
LYS 970 3TW7.subA 2 0 1 1 0 0 0 0 0 1 1 0 0 0 0 0 1 0 0 2 1004
LYS 971 3TW7.subA 4 0 0 0 0 0 0 0 0 4 0 1 0 1 0 0 3 0 0 2 1004
LYS 989 3TW7.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 0 2 1004
LYS 1028 3TW7.subA 2 0 1 1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 0 2 1004
LYS 1030 3TW7.subA 3 0 1 1 0 0 0 0 0 2 1 1 0 1 1 0 0 0 0 2 1004
LYS 1062 3TW7.subA 6 1 0 1 0 0 0 0 0 5 1 1 0 1 0 0 4 0 0 2 1004
LYS 6 3UAG.subA 7 2 0 2 0 0 0 0 0 5 2 0 2 2 1 1 1 0 3 2 430
LYS 45 3UAG.subA 4 1 0 1 0 0 0 0 0 3 1 0 0 0 1 0 2 0 3 2 430
LYS 115 3UAG.subA 14 0 1 1 0 0 0 0 0 13 1 5 2 7 2 0 3 1 4 2 430
LYS 127 3UAG.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 3 0 2 0 7 2 430
LYS 206 3UAG.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 0 2 0 2 2 430
LYS 251 3UAG.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 3 2 430
LYS 254 3UAG.subA 8 0 1 1 0 1 1 1 2 5 3 0 1 1 0 1 3 0 1 2 430
LYS 259 3UAG.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 4 2 430
LYS 262 3UAG.subA 7 0 1 1 0 0 0 0 0 6 1 3 0 3 0 1 2 0 2 2 430
LYS 291 3UAG.subA 6 0 0 0 0 0 0 0 0 6 0 2 0 2 3 0 1 0 2 2 430
LYS 319 3UAG.subA 10 1 1 2 0 1 1 1 2 6 4 3 1 4 1 1 0 0 6 2 430
LYS 348 3UAG.subA 11 3 0 3 0 0 1 0 1 7 4 2 1 3 3 0 1 1 11 2 430
LYS 420 3UAG.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 1 0 0 0 2 430
LYS 434 3UAG.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 4 2 430
LYS 2 3UBP.subC 2 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 2 0 1 2 570
LYS 33 3UBP.subC 4 1 1 2 0 0 0 0 0 2 2 0 0 0 0 2 0 0 10 2 570
LYS 48 3UBP.subC 7 0 0 0 0 0 0 0 0 7 0 0 1 1 3 1 2 0 2 2 570
LYS 84 3UBP.subC 12 1 0 1 0 0 0 0 0 11 1 1 1 2 4 3 2 0 4 2 570
LYS 90 3UBP.subC 11 2 1 3 0 1 0 1 1 7 4 1 1 2 0 1 4 0 8 2 570
LYS 99 3UBP.subC 14 3 1 4 0 0 0 0 0 10 4 3 0 3 3 1 3 0 5 2 570
LYS 127 3UBP.subC 11 1 2 3 0 1 0 1 1 7 4 0 0 0 3 1 3 0 5 2 570
LYS 169 3UBP.subC 8 0 1 1 0 0 0 0 0 7 1 1 0 1 5 0 1 0 7 2 570
LYS 182 3UBP.subC 5 0 1 1 1 0 0 1 1 3 2 0 1 1 0 1 1 0 3 2 570
LYS 185 3UBP.subC 6 0 1 1 1 0 0 1 1 4 2 1 0 1 0 0 3 0 4 2 570
LYS 199 3UBP.subC 8 0 0 0 0 0 1 0 1 7 1 0 1 1 3 1 1 0 0 2 570
LYS 326 3UBP.subC 6 1 0 1 0 0 1 0 1 4 2 0 2 2 0 0 2 0 2 2 570
LYS 383 3UBP.subC 12 2 1 3 0 0 0 0 0 9 3 1 2 3 3 1 2 0 2 2 570
LYS 385 3UBP.subC 11 3 0 3 0 1 0 1 1 7 4 0 0 0 4 1 2 0 4 2 570
LYS 386 3UBP.subC 4 1 1 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 7 2 570
LYS 395 3UBP.subC 3 0 2 2 0 0 0 0 0 1 2 0 1 1 0 0 0 0 4 2 570
LYS 404 3UBP.subC 10 1 1 2 0 1 0 1 1 7 3 0 1 1 1 3 2 0 5 2 570
LYS 409 3UBP.subC 14 1 0 1 0 1 1 1 2 11 3 3 1 4 1 2 4 0 0 2 570
LYS 431 3UBP.subC 12 0 3 3 0 0 1 0 1 8 4 1 0 1 3 1 3 0 8 2 570
LYS 441 3UBP.subC 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 1 0 0 4 2 570
LYS 446 3UBP.subC 7 0 0 0 0 0 0 0 0 7 0 0 1 1 2 3 1 0 2 2 570
LYS 452 3UBP.subC 14 1 0 1 0 0 0 0 0 13 1 2 2 4 3 1 5 0 2 2 570
LYS 497 3UBP.subC 8 1 0 1 1 0 0 1 1 6 2 3 1 4 1 0 1 0 11 2 570
LYS 507 3UBP.subC 5 0 1 1 0 0 0 0 0 4 1 0 0 0 2 0 2 0 11 2 570
LYS 511 3UBP.subC 5 0 1 1 0 1 1 1 2 2 3 0 0 0 1 0 1 0 11 2 570
LYS 518 3UBP.subC 7 1 0 1 1 0 1 1 2 4 3 1 1 2 1 0 1 0 11 2 570
LYS 525 3UBP.subC 16 1 0 1 0 0 1 0 1 14 2 2 1 3 6 0 5 0 3 2 570
LYS 526 3UBP.subC 2 1 0 1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 8 2 570
LYS 529 3UBP.subC 8 1 0 1 0 0 0 0 0 7 1 1 1 2 0 2 3 0 9 2 570
LYS 547 3UBP.subC 7 1 2 3 0 0 0 0 0 4 3 0 0 0 1 0 3 0 4 2 570
LYS 559 3UBP.subC 3 0 1 1 0 0 0 0 0 2 1 0 0 0 0 0 2 0 6 2 570
LYS 6 4UAG.subA 9 2 1 3 0 0 0 0 0 6 3 0 2 2 1 1 2 0 2 2 429
LYS 45 4UAG.subA 5 1 0 1 0 0 0 0 0 4 1 0 0 0 2 0 2 0 0 2 429
LYS 115 4UAG.subA 14 0 1 1 0 0 0 0 0 13 1 5 2 7 2 0 3 0 4 2 429
LYS 127 4UAG.subA 6 0 1 1 0 0 0 0 0 5 1 0 0 0 3 0 2 0 6 2 429
LYS 206 4UAG.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 0 2 0 3 2 429
LYS 251 4UAG.subA 3 1 0 1 0 0 0 0 0 2 1 0 0 0 1 0 1 0 0 2 429
LYS 254 4UAG.subA 6 0 1 1 0 1 0 1 1 4 2 0 1 1 0 1 2 0 1 2 429
LYS 259 4UAG.subA 4 0 2 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 5 2 429
LYS 262 4UAG.subA 7 0 1 1 0 0 0 0 0 6 1 3 0 3 0 1 2 0 3 2 429
LYS 291 4UAG.subA 6 0 0 0 0 0 0 0 0 6 0 2 0 2 3 0 1 0 2 2 429
LYS 319 4UAG.subA 9 1 0 1 0 1 1 1 2 6 3 2 0 2 2 1 1 0 5 2 429
LYS 348 4UAG.subA 12 3 0 3 0 0 1 0 1 8 4 2 1 3 4 0 1 0 8 2 429
LYS 420 4UAG.subA 3 0 0 0 0 0 0 0 0 3 0 0 2 2 0 1 0 0 0 2 429
LYS 434 4UAG.subA 6 0 1 1 0 1 0 1 1 4 2 0 0 0 2 0 2 0 1 2 429
LYS 2 4UBP.subC 2 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 2 0 3 2 570
LYS 33 4UBP.subC 4 1 1 2 0 0 0 0 0 2 2 0 0 0 0 2 0 0 8 2 570
LYS 48 4UBP.subC 7 0 0 0 0 0 0 0 0 7 0 0 1 1 3 1 2 0 7 2 570
LYS 84 4UBP.subC 12 1 0 1 0 0 0 0 0 11 1 1 1 2 4 3 2 0 4 2 570
LYS 90 4UBP.subC 11 2 1 3 0 1 0 1 1 7 4 1 1 2 0 1 4 0 10 2 570
LYS 99 4UBP.subC 14 3 1 4 0 0 0 0 0 10 4 3 0 3 3 1 3 0 5 2 570
LYS 127 4UBP.subC 11 1 2 3 0 1 0 1 1 7 4 0 0 0 3 1 3 0 4 2 570
LYS 169 4UBP.subC 7 0 1 1 0 0 0 0 0 6 1 1 0 1 4 0 1 0 7 2 570
LYS 182 4UBP.subC 5 0 1 1 1 0 0 1 1 3 2 0 1 1 0 1 1 0 2 2 570
LYS 185 4UBP.subC 6 0 1 1 1 0 0 1 1 4 2 1 0 1 0 0 3 0 4 2 570
LYS 199 4UBP.subC 8 0 0 0 0 0 1 0 1 7 1 0 1 1 3 1 1 0 0 2 570
LYS 326 4UBP.subC 6 1 0 1 0 0 1 0 1 4 2 0 2 2 0 0 2 0 0 2 570
LYS 383 4UBP.subC 12 2 1 3 0 0 0 0 0 9 3 1 2 3 3 1 2 0 2 2 570
LYS 385 4UBP.subC 11 3 0 3 0 1 0 1 1 7 4 0 0 0 4 1 2 0 3 2 570
LYS 386 4UBP.subC 4 1 1 2 0 0 0 0 0 2 2 0 1 1 0 0 1 0 9 2 570
LYS 395 4UBP.subC 7 0 3 3 1 0 0 1 1 3 4 1 2 3 0 0 0 0 3 2 570
LYS 404 4UBP.subC 12 1 1 2 1 1 0 2 2 8 4 0 1 1 1 3 3 0 3 2 570
LYS 409 4UBP.subC 13 1 0 1 0 1 1 1 2 10 3 3 1 4 1 2 3 0 0 2 570
LYS 431 4UBP.subC 12 0 3 3 0 0 1 0 1 8 4 1 0 1 3 1 3 0 6 2 570
LYS 441 4UBP.subC 3 0 1 1 0 0 0 0 0 2 1 0 0 0 1 1 0 0 4 2 570
LYS 446 4UBP.subC 7 0 0 0 0 0 0 0 0 7 0 0 1 1 2 3 1 0 2 2 570
LYS 452 4UBP.subC 13 1 0 1 0 0 0 0 0 12 1 2 2 4 3 1 4 0 3 2 570
LYS 497 4UBP.subC 8 1 0 1 1 0 0 1 1 6 2 3 1 4 1 0 1 0 10 2 570
LYS 507 4UBP.subC 5 0 1 1 0 0 0 0 0 4 1 0 0 0 2 0 2 0 10 2 570
LYS 511 4UBP.subC 5 0 1 1 0 1 1 1 2 2 3 0 0 0 1 0 1 0 7 2 570
LYS 518 4UBP.subC 7 1 0 1 1 0 1 1 2 4 3 1 1 2 1 0 1 0 10 2 570
LYS 525 4UBP.subC 16 1 0 1 0 0 1 0 1 14 2 2 1 3 6 0 5 0 2 2 570
LYS 526 4UBP.subC 4 1 0 1 0 0 0 0 0 3 1 0 1 1 1 0 1 0 3 2 570
LYS 529 4UBP.subC 8 1 0 1 0 0 0 0 0 7 1 1 1 2 0 2 3 0 8 2 570
LYS 547 4UBP.subC 6 1 1 2 0 0 0 0 0 4 2 0 0 0 1 0 3 0 3 2 570
LYS 559 4UBP.subC 3 0 1 1 0 0 0 0 0 2 1 0 0 0 0 0 2 0 6 2 570
LYS 14 8RUC.subA 4 0 0 0 1 0 0 1 1 3 1 0 0 0 2 1 0 0 2 2 467
LYS 18 8RUC.subA 7 1 0 1 1 0 0 1 1 5 2 2 0 2 1 0 2 0 6 2 467
LYS 21 8RUC.subA 7 1 2 3 0 0 0 0 0 4 3 0 0 0 1 2 1 0 0 2 467
LYS 81 8RUC.subA 10 1 0 1 0 1 0 1 1 8 2 4 0 4 1 2 1 0 1 2 467
LYS 128 8RUC.subA 3 0 0 0 0 0 0 0 0 3 0 0 0 0 2 1 0 0 0 2 467
LYS 146 8RUC.subA 5 0 0 0 0 0 0 0 0 5 0 1 0 1 1 0 3 0 3 2 467
LYS 161 8RUC.subA 6 1 1 2 0 0 1 0 1 3 3 0 0 0 1 0 2 0 2 2 467
LYS 164 8RUC.subA 9 1 1 2 0 2 0 2 2 5 4 0 1 1 0 1 3 0 3 2 467
LYS 175 8RUC.subA 8 1 0 1 1 0 0 1 1 6 2 1 0 1 3 0 2 1 1 2 467
LYS 177 8RUC.subA 6 1 1 2 1 0 0 1 1 3 3 0 1 1 1 0 1 1 3 2 467
LYS 183 8RUC.subA 5 0 0 0 0 1 0 1 1 4 1 1 1 2 1 1 0 0 1 2 467
LYS 227 8RUC.subA 5 0 1 1 0 0 0 0 0 4 1 0 0 0 2 2 0 0 3 2 467
LYS 236 8RUC.subA 12 1 1 2 0 1 0 1 1 9 3 1 1 2 3 0 4 0 2 2 467
LYS 252 8RUC.subA 5 1 1 2 0 1 0 1 1 2 3 0 0 0 0 1 1 0 7 2 467
LYS 305 8RUC.subA 8 1 0 1 0 2 0 2 2 5 3 1 2 3 0 0 2 0 5 2 467
LYS 316 8RUC.subA 11 1 0 1 0 2 0 2 2 8 3 0 1 1 2 1 4 0 1 2 467
LYS 334 8RUC.subA 5 0 0 0 0 0 0 0 0 5 0 0 0 0 3 1 1 0 5 2 467
LYS 356 8RUC.subA 7 1 1 2 0 1 0 1 1 4 3 1 0 1 0 2 1 0 9 2 467
LYS 450 8RUC.subA 5 0 1 1 0 1 0 1 1 3 2 1 0 1 1 1 0 0 4 2 467
LYS 463 8RUC.subA 6 0 2 2 0 1 0 1 1 3 3 0 0 0 1 2 0 0 6 2 467
LYS 466 8RUC.subA 4 0 2 2 0 0 0 0 0 2 2 0 0 0 0 1 1 0 2 2 467
";

	my $sizecontrol = 0;
	my $c = 0;
	
	my %hashofarrays = ();

	my @data = split('\n',$lystraining);

	foreach my $line (@data)
	{
		chomp($line);
		my @temp = split(" ",$line);
		my $size = @temp;
		#Size control for the first one
		if ($c == 0)
		{
			$sizecontrol = $size;
		}
		
		#All the others - -- - - - -- ---
		if($size != $sizecontrol)
		{
			print "\n\nFATAL ERROR PROBLEM\n";
			print $command." error 518\n";
			print "Please, contact us at prelyscar\@gmail.com\n\n";
			exit;
		}
		#- - - - - - -- - -end of control
		
		$hashofarrays{$c} = [@temp];
		
		$c++;
	}
	return (\%hashofarrays);
}

