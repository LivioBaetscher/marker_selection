#!/usr/bin/perl

## Project: PhD
## Script purpose: This script can be used in combination or as followup script
##  of the script 'select_markers.sh'. It organizes processes and can filter the
##  BLAST results collected from 'select_markers.sh'.
## Date: 01.02.2021
## Author: Livio Bätscher

#**************************** Variable input *****************************#

## The programm needs the directory with the results of a certain gen set
##  as input. Furthermore, a outputdirectory needs to be declaired. Additonally,
##  a e value for filtering as well as a number of transcriptomes in which
##  the gene needs to be present can be declaired. The user has also the chance
##  to choose a minimal length that a sequence needs to have.

#**************************** Main Programm *****************************#
# Import the packages
use strict;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use Storable 'dclone';
use Data::Dumper;

my $nrTrSelection = 1; #Set the default for the number of transcriptomes
my $lengthSelection = 1; #Set the default for the minimal sequence length

# Process the input arguments
my ($genSetDir, $outDir, $eValSelection, $pidentSelection) = ();
GetOptions(
	"in|i=s" => \$genSetDir,
	"out|o=s"=> \$outDir,
	"eVal:s" => \$eValSelection,
	"len|l:s" => \$lengthSelection,
	"ident:s" => \$pidentSelection,
	"Tr:s" => \$nrTrSelection
);

#Error handling for missing input directory
unless(defined $genSetDir and -s $genSetDir){
	die "###########################################################################\n".
	"Input directory is missing or incorrect.\n".
	"###########################################################################\n";
}

#Error handling for missing output file
unless(defined $outDir){
	die "###########################################################################\n".
	"Output directory is missing.\n".
	"###########################################################################\n";
}

#Error handling for the eValue
if ($eValSelection) {
	unless (looks_like_number($eValSelection) and $eValSelection < 1 and $eValSelection > 0) {
		die "###########################################################################\n".
		"E value needs to be a positive number below 1.\n".
		"###########################################################################\n";
	}
}

#Error handling for the percent identity
if ($pidentSelection) {
	unless (looks_like_number($pidentSelection) and $pidentSelection < 100 and $pidentSelection > 0) {
		die "###########################################################################\n".
		"The percent identity needs to be a positive number between 0 and 100.\n".
		"###########################################################################\n";
	}
}

#Load all the files in the directory
opendir my $dir, $genSetDir or die "Cannot open directory: $!.";
my @filesArray = readdir $dir;
closedir $dir;

#Change the directory for simplisity
chdir $genSetDir;

#Remove the '.' and '..' file entries (first tow entries)
splice @filesArray, 0, 2;

#Error handling for the number of transcriptomes
if ($nrTrSelection) {
	unless (looks_like_number($nrTrSelection) and $nrTrSelection >= 1 and $nrTrSelection <= scalar @filesArray) {
		die "###########################################################################\n".
		"Number of transcriptome needs to be a natural number with a maximum that\nis equal to the number of files in the input directory.\n".
		"###########################################################################\n";
	}
}

#Get the BLAST results from the files stored in a hash
my %fileHash = create_file_hash(@filesArray);

#Converte the file hash to a gene hash
my %geneHash = convert_to_gene_hash(%fileHash);

#generate a result hash that is identical with the gene hash
my %resultHash = %geneHash;

#Print a print a short summary
printf "%-25s | %-10s | %-13s | %-12s\n", "Before filtering";
print_summary(%resultHash);

#Handle the selection by a certain eValue
if ($pidentSelection) {
	print "Select for BLAST hits with an percent identity that is above or equal: "
	 . $pidentSelection . "\n";
	#Call the selection function
	%resultHash = do_ident_selection(%resultHash);
}

#Handle the selection by a certain eValue
if ($eValSelection) {
	print "Select for BLAST hits with an eValue that is below or equal: "
	 . $eValSelection . "\n";
	#Call the selection function
	%resultHash = do_eValue_selection(%resultHash);
}

#Handle the selection by length transcriptomes
print "AND\nSelect for genes with the minimum length of: " . $lengthSelection .
	 "\n\n";
%resultHash = do_length_selection(%resultHash); #Call the selection function

#Handle the selection by number of transcriptomes
print "AND\nSelect for genes present in at least: " . $nrTrSelection .
	 " transcriptomes\n";
%resultHash = do_transcriptome_selection(%resultHash); #Call the selection function

#Print a print a short summary
printf "%-25s | %-10s | %-13s | %-12s\n", "After filtering";
print_summary(%resultHash);

#Generate a fasta hash containing all 12 species
my %fastaHash = create_fasta_hash($genSetDir);

#Create the output directory and move into it
unless(mkdir $outDir) {
	die "Unable to create $outDir\n";
}
chdir $outDir;

#Print the sequences into a file
print_sequences(\%resultHash, \%fastaHash);

#**************************** Subroutines *****************************#
## This subroutine reads the ".txt" file and splits the lines.
sub read_lines {
	my ($inFile) = @_; #Pass the argument

 	#Try to open the file and read the content
	open my $in, "<:encoding(utf8)", $inFile or die "$inFile: $!";
	local $/ = undef;
	my $content = <$in>;
	close $in;

	#Split the content into lines
  return split(/\n/, $content);
}


## This subroutine reads a file array and returns a hash object with the file as
##  key and the BLAST hits as array elements
sub create_file_hash {
  my (@files) = @_; #Pass the array

  my %fileHash; #Create a empty hash
  #Loop through all the files
  foreach my $file (@files) {
    print "Processing file: " . $file . "\n"; #Printstatement 1

    #Load the blast hits as array
    my @blastHits = read_lines($file);

    #Split each Blast line into a array
    foreach my $hit (@blastHits) {
      my @hit = split(/\t/, $hit);
      $hit = \@hit;
    }

    #Get the name of the genetic_mat
    my $geneticMatName = $file;
    $geneticMatName =~ s/blast_//;
    $geneticMatName =~ s/_cont//;
    $geneticMatName =~ s/_genomic//;
		$geneticMatName =~ s/_\Q$genSetDir//;
		$geneticMatName =~ s/.txt//;

    $fileHash{$geneticMatName} = \@blastHits;
  }

  #Retun the hash
  return %fileHash;
}

## This subroutine reads a file hash and converte it into a gene hash. Where
##  the genes are key. In this hash is another hash with the species as key.
##  In a first iteration only the longest elment above the e value of 0.001
##  is added. In a secand iteration genes don't overlap are added and with the
##  gen key + a unique didgit number. So that each gen key remains unique.
sub convert_to_gene_hash {
  my (%inputHash) = @_; #Pass the hash
  print "\nGenerate a summary hash\n\n";
  my %genes; #Create a empty hash

  #Loop through all files
  foreach my $species (keys %inputHash) {
    #Loop through each element of the BLAST hits
    foreach my $hit (@{$inputHash{$species}}) {
      #Restore the variable to a proper array
      my @hitArray = @{$hit};

      #Create variables needed to evaluate the a gene already in the hash
      my $genKey = $hitArray[1];
			my $newLength = $hitArray[3];
      my $newEVal = $hitArray[10];
      splice @hitArray, 1, 1; #Remove redundant information

      #Check if the key combination of gene and species already exists
      if (exists $genes{$genKey}{$species}) {
				#Yes, so check if a replacement is needed

				my $currentLength = $genes{$genKey}{$species}[2]; #Get current length

        #Check if the new evalue is smaller than the existing one
        if (($newEVal <= 0.001)  and ($newLength > $currentLength)) {
          #Replace the old element
          $genes{$genKey}{$species} = \@hitArray;
        }
      } else {
				#The combination of keys does not exist, add the element
        $genes{$genKey}{$species} = \@hitArray;
      }
    }
  }

	#Loop again through all files and search for genes that are detected but witht
	# a gap in between and add them too.
  foreach my $species (keys %inputHash) {
		#Loop through each element of the BLAST hits
    foreach my $hit (@{$inputHash{$species}}) {
      #Restore the variable to a proper array
      my @hitArray = @{$hit};

      #Create variables needed to evaluate the BLAST hits a secand time
			my $newContig = $hitArray[0];
			$newContig =~ s/:.*//;
			my $genKey = $hitArray[1];
			my $newStart = $hitArray[6];
			my $newEnd = $hitArray[7];
      splice @hitArray, 1, 1; #Remove redundant information

			#Create the variables that are currently hold these gen and species key
			my $currentContig = $genes{$genKey}{$species}[0];
			$currentContig =~ s/:.*//;
			my $currentStart = $genes{$genKey}{$species}[5];
			my $currentEnd = $genes{$genKey}{$species}[6];

			#If another BLAST hit does not overlap with the first hit add the hit
			if (($newEnd <= $currentStart) or ($currentEnd <= $newStart)) {
				#Add a number at the end of a gen name (must be unique)
				my $number = 1;

				while (exists $genes{$genKey}{$species}) {
					$genKey =~ s/-\Q$number//;
					$number = $number + 1;
					$genKey = $genKey . "-" . $number;
				}

				$genes{$genKey}{$species} = \@hitArray; #Add the element
			}
		}
	}
	#print Dumper sort \%genes;
	#Return the hash
  return %genes;
}


## This subroutine takes the result hash aplies filtering according to the
##  percent identity (input -ident). A BLAST hit needs to at lest the similarity
##  definedin the input to pass the filtering.
sub do_ident_selection {
	my (%filterHash) = @_;
	#Loop through all genes
	foreach my $gene (keys %filterHash) {
		#Loop through all species
		foreach my $species (keys %{$filterHash{$gene}}) {
			#Restore as a a proper array
			my @currentArray = @{$filterHash{$gene}{$species}};
			my $pident = $currentArray[1]; #Get the current percent identit

			#Remove the element if the filter criteria is not fulfilled
			if ($pidentSelection > $pident) {
				delete $filterHash{$gene}{$species};
			}
		}
	}

	#Return the hash
	return %filterHash;
}


## This subroutine takes the result hash aplies filtering according to the
##  eValue (input -eVal). A BLAST hit needs to be in at lest the eValue defined
##  in the input to pass the filtering.
sub do_eValue_selection {
	my (%filterHash) = @_;
	#Loop through all genes
	foreach my $gene (keys %filterHash) {
		#Loop through all species
		foreach my $species (keys %{$filterHash{$gene}}) {
			#Restore as a a proper array
			my @currentArray = @{$filterHash{$gene}{$species}};
			my $eVal = $currentArray[9]; #Get the current eValue

			#Remove the element if the filter criteria is not fulfilled
			if ($eValSelection < $eVal) {
				delete $filterHash{$gene}{$species};
			}
		}
	}

	#Return the hash
	return %filterHash;
}


## This subroutine takes the result hash aplies filtering according to the
##  number of transcriptomes (input -Tr). A gen needs to be in at lest the
##  number of transcriptomes defined in the input to pass the filtering.
sub do_transcriptome_selection {
	my (%filterHash) = @_;
	#Loop through all genes
	foreach my $gene (keys %filterHash) {
		#In how many genes is the value present
		my $nrTr = scalar keys %{$filterHash{$gene}};

		#Remove the element if the filter criteria is not fulfilled
		if ($nrTr < $nrTrSelection) {
			delete $filterHash{$gene};
		}
	}

	#Return the hash
	return %filterHash;
}


## This subroutine takes the result hash aplies filtering according to the
##  length of a hit. A BLAST hit needs to be in at lest the length defined
##  in the input to pass the filtering.
sub do_length_selection {
	my (%filterHash) = @_;
	#Loop through all genes
	foreach my $gene (keys %filterHash) {
		#Loop through all species
		foreach my $species (keys %{$filterHash{$gene}}) {
			#Restore as a a proper array
			my @currentArray = @{$filterHash{$gene}{$species}};
			#Get the length
			my $length = $currentArray[2];

			#Remove the element if the filter criteria is not fulfilled
			if ($length < $lengthSelection) {
				delete $filterHash{$gene}{$species};
			}
		}
	}

	#Return the hash
	return %filterHash;
}


## This subroutine reads the gebne hash returns a short print statement in a
##  good looking table.
sub print_summary {
	my (%printHash) = @_;

	printf "%-25s | %-10s | %-13s | %-12s\n",
					 "Gen:", "Presences:", "Mean e value:", "Mean length:";
	printf "#####################################################################\n";

	#Loop through each gene
	foreach my $gene (sort keys %printHash) {
		#Get the number of Transcriptome hits in each gene
		my $trNr =  scalar keys %{$printHash{$gene}}; #After filtering

		#Create two empty variables
		my $meaneVal;
		my $meanLen;

		foreach my $species (sort keys %{$printHash{$gene}}) {
			#Restore as a a proper array
			my @currentArray = @{$printHash{$gene}{$species}};

			#Sum up the elements
			$meaneVal = $meaneVal + $currentArray[9];
			my $length = $currentArray[2];
			$meanLen = $meanLen + $length;
		}

		#Error handling
		if ($trNr > 0) {
			#Divide the elements
			$meaneVal = $meaneVal / $trNr;
			$meanLen =  $meanLen / $trNr;
		}

		printf "%-25s | %-10.0f | %-13e | %-12.2f\n",
						 $gene, $trNr, $meaneVal, $meanLen;

	}
	print "\n";
}


## This subroutine loads all the 12 fasta files in the folder 'genetic_mat'
##  and stores them in a hash, it takes the initial directory as input, in order
##  to be able to move back to this directory.
sub create_fasta_hash {
	my ($oldDir) = @_;
	#Set the directory witht the fastafiles
	my $fastaDir = "../../../genetic_mat/";

	#Load all the files in the directory
	opendir my $dir, $fastaDir or die "Cannot open directory: $!.";
	my @fastaArray = readdir $dir;
	closedir $dir;

	#Remove the '.' and '..' file entries (first tow entries)
	splice @fastaArray, 0, 2;

	#Change the directory for simplisity
	chdir $fastaDir;
	my %fastaSpecies; #Create a empty hash

	#Loop through the fasta array
	foreach my $fastaFile (@fastaArray) {
		#Generate a fasta hash for each file
		my %currentFastaHash = read_fasta_as_hash($fastaFile);
		print "Load file: " . $fastaFile . "\n"; #Print the current state

		#Chang the name, so that it matches the result hash
		my $species = $fastaFile;
		$species =~ s/_cont//;
    $species =~ s/_genomic//;
		$species =~s/.fa//;

		#Add the element to the fasta hash
		$fastaSpecies{$species} = \%currentFastaHash;
	}

	chdir "../potential_markers/results/" . $genSetDir; #Move back to the inital directory

	#Return the hash
	return %fastaSpecies;
}


## This subroutine splits fasta files into correct hashes from filenames, as
##  input this subroutine, gets the target file name.
sub read_fasta_as_hash{
  my ($file) = @_; #Pass the argument

  #Instanciate empty variables
  my $id = '';
  my %sequences;

  #Try to open the file and read the content
  open my $in, "<:encoding(utf8)", $file or die "$file: $!";

  #If the file has lines
  while( my $line = <$in> ){
        chomp $line; #Split it
        #Create a id if it starts with '>'
        if ( $line =~ /^(>.*)$/ ){
            $id = $1;
						$id =~ s/>//;
						$id =~ s/\s.+//;
        # Else add the line to the content
        } else {
             $sequences{$id} .= $line;
        }
   }
   close $in;

   #Return the hash
   return %sequences;
}


## This subroutine takes the filtered result hash and the fasta hash and
##  and calls the sequences printing. Therefore the sequences get stored in
##  files according to thier gen name. All the finished files will then be
##  aligned by using MAFT.
sub print_sequences {
	my ($finalHash, $fastaHash) = @_;
	my %finalHash = %{$finalHash};
	my %fastaHash = %{$fastaHash};

	print "\nPrint into FASTA file: " . $genSetDir . $outDir . "\n\n";

	#Loop through all the genes
	foreach my $gene (keys %finalHash) {
		#Create a file name
		my $file = $gene . ".fa";

		#Loop through all the species
		foreach my $species (keys %{$finalHash{$gene}}) {
			#Restore as a a proper array
			my @blastHit = @{$finalHash{$gene}{$species}};
			#Instanciate some variables
			my $qseqid = $blastHit[0];
			my $qstart  = $blastHit[5];
			my $qend = $blastHit[6];

			#Call the sequence
			my $seq = $fastaHash{$species}{$qseqid};

			$seq = substr($seq, $qstart, ($qend-$qstart)); #Cut out the match fragment
			my $id = get_id_name($gene, $species, $qseqid); #Create a id
			add_to_fasta_file($file, $id, $seq);
		}

		#If the file exists
		if (-e $file) {
		print "Creating a alignment for the file: " . $file . "\n";

			#Create a file name for the alignment
			my $alingmentFile = $gene . "_aligned.fa";

			#Align the FASTA files
			system ("mafft --maxiterate 1000 --localpair --adjustdirection --quiet $file > $alingmentFile") == 0
				or die "Error: mafft failed. Make sure mafft is installed and available.";
		}
	}
}

## This subroutine generates a sequence name so that they  match the
##  requirements of Arbor Bioscience.
sub get_id_name {
	my ($gene, $species, $qseqid) = @_;

	my $newId = $gene . "_" . $species . "_" . $qseqid;
	#Replace the species id by a proper species name
	$newId =~ s/ath://;
	$newId =~ s/sly://;
	$newId =~ s/SLY://;

	#Replace the species id by a proper species name
	$newId =~ s/_PRJNA238546_SRR7346504_/-P-sikkimensis-/;
	$newId =~ s/_PRJNA238546_SRR6830996_/-P-maximowiczii-/;
	$newId =~ s/_GBRY01_1_/-P-vulgaris-/;
	$newId =~ s/_PRJNA238546_SRR2039591_/-P-chrysochlora-/;
	$newId =~ s/_PRJNA238546_SRR3355026_/-P-forbesii-/;
	$newId =~ s/_PRJNA238546_SRR3355043_/-P-veris-/;
	$newId =~ s/_PRJNA238546_SRR640158_/-P-wilsonii-/;
	$newId =~ s/_PRJNA238546_SRR3307913_/-P-sinensis-/;
	$newId =~ s/_PRJNA238546_SRR5377219_/-P-ovalifolia-/;
	$newId =~ s/_PRJNA238546_SRR629689_/-P-poissonii-/;
	$newId =~ s/_PRJNA238546_SRR866502_/-P-obconica-/;
	$newId =~ s/_GCA_000788445_1_ASM78844v1_/-P-veris-/;
	$newId =~ s/_CMFF_.*/-L-polyphyllus/;
	$newId =~ s/_CYVA_.*/-C-racemosa/;
	$newId =~ s/_GBVZ_.*/-T-thalictroides/;
	$newId =~ s/_TTRG_.*/-L-angustifolius/;
	$newId =~ s/_UPOG_.*/-A-pulsatilla/;
	$newId =~ s/_VGHH_.*/-H-canadensis/;
	$newId =~ s/_ZUHO_.*/-A-hupehensis/;

	#Remove unnecesarry token
	$newId =~ s/.1:/-/;
	$newId =~ s/.1$//;

	#Sorten the onetowtree gene names
	$newId =~ s/Androsace_/A-/;
	$newId =~ s/Cortusa_/C-/;
	$newId =~ s/Dionysia_/D-/;
	$newId =~ s/Primula_/P-/;
	$newId =~ s/_Cluster_/-/;


	#Remove unnecesarry token in the s locus gen name
	$newId =~ s/-T1//;

	return $newId;
}


## This subroutine creates for each gene a fasta file and adds the sequences in
##  FASTA format.
sub add_to_fasta_file {
	my ($fileName, $seqId, $fastaSeq) = @_;

	#Open a file with the gen name
	open my $outFasta, ">>", $fileName or die "Error: opening process failed at the file $!";

	#Print the sequence into the file according to the FASTA-format
	print $outFasta ">" . $seqId . "\n";
	print $outFasta $fastaSeq . "\n";

	close $outFasta; #Close the file
}
