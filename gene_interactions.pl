#!/usr/bin/perl
use warnings; use strict; use autodie;
use Getopt::Long;
#
#@author John 
# johnwilliams@post.harvard.edu
# j.williams@har.mrc.ac.uk
#
# if you encounter any problems, contact John.
#

=pod
					if errors of the type:
		Use of uninitialized value $first_sig in split at gene_interactions.pl line xxx.
		Use of uninitialized value in split at gene_interactions.pl line xxx+2.
		Use of uninitialized value in split at gene_interactions.pl line xxx+3.
		Use of uninitialized value in pattern match (m//) at gene_interactions.pl line xxx+19.
		Use of uninitialized value in pattern match (m//) at gene_interactions.pl line xxx+36.

Appear, do not worry: there is a bin size with no interactions for you given gene. You can confirm this by looking at the interaction files generated in the gene folder.

=cut

my $chrs;

my $help = 0;

my $circos = 0;

my $gene;

my $genome;

GetOptions (
"gene=s" => \$gene,
"help=s" => \$help,
"circos" => \$circos,
"genome=s" => \$genome
);


sub helper { die(qq/
		This script must be run from the project directory of your Hi-C experiment. For example, if you ran in directory "project", cd into project then run this script. 
		The script creates functional annotation of interactions, tables which can be copied to Excel for display, and circos visualizations
		The Circos .conf files and R scripts for prommoter annotation remain, you can change the R script and re-run it, or change the Circos --conf file and then re-run circos on that file
		_______________________________________________
		
    USAGE : perl gene_interactions.pl <arguments> 
			ex: perl gene_interactions.pl --gene Klf14 --circos
    ARGUMENTS :  
                    REQUIRED 
           
					-- gene: the name of the gene whose promoter interactions you would like to annotate and visualize
					-- genome: either mm9 or mm10, for correct annotation and circos plotting
					
					OPTIONAL
					
					-- circos: create circos diagrams
					-- help: display this message
					

               \n/); 
} 
 
if ( $help) {
 &helper();
 } 

unless ($gene) {
 print "\n MISSING ARGUMENTS : Specify a gene name\n";
 &helper();
 } 


unless ($genome) {
 print "\n MISSING ARGUMENTS : Specify either mm10 or mm9 genome\n";
 &helper();
 } 




&get_gene($gene);

chdir("./$gene");

# make write_table() script to make human readible table and to annotate with GO annotation
print "Functionally annotating interactions with $gene\n";
&write_table($gene);
print "Done with $gene functional annotation\n";

print "Genome is $genome\n";

if ($circos) {
	print "writing interactions\n";
	&writeinter($gene);
	print "done writing interactions\n";
	print "writing circos config files\n";
	&writeconf($gene,$genome);
	print "done writing circos files\n";
	print "making circos diagrams\n";
	&make_diagrams($gene);
	print "done making circos diagrams\n";
}

## 
# begin subroutines
##
sub get_gene{

	my $myGene = $_[0];
	
	# write bash script to run, since couldn't get permissions correct with chmod or mkdir in perl
	
	open (EXTRACT, ">extract.sh");
		print EXTRACT "mkdir -m 777 ./$myGene\n";
		print EXTRACT "grep -w '$myGene' ./interactions/sigints100k/sigints100k.annotation.txt >> $myGene/100k.txt\n";
		print EXTRACT "grep -w '$myGene' ./interactions/sigints1M/sigints1M.annotation.txt >> $myGene/1M.txt\n";
		print EXTRACT "grep -w '$myGene' ./interactions/sigints20k/sigints20k.annotation.txt >> $myGene/20k.txt\n";
		print EXTRACT "grep -w '$myGene' ./interactions/sigints40k/sigints40k.annotation.txt >> $myGene/40k.txt\n";
		print EXTRACT "grep -w '$myGene' ./interactions/sigints70k/sigints70k.annotation.txt >> $myGene/70k.txt\n";
		print EXTRACT "grep -w '$myGene' ./interactions/sigints250k/sigints250k.annotation.txt >> $myGene/250k.txt\n";
		print EXTRACT "mkdir -m 777 ./$myGene/Tables\n";
		print EXTRACT "mkdir -m 777 ./$myGene/Graphs\n";
	close(EXTRACT);
	
	system(q(chmod 777 ./extract.sh));
	system(q(bash extract.sh));
	# check each file size; if all return empty then gene is not in set; die;
	my $k20k = -s "$myGene/20k.txt";
	my $k40k = -s "$myGene/40k.txt";
	my $k70k = -s "$myGene/70k.txt";
	my $k100k = -s "$myGene/100k.txt";
	my $k250k = -s "$myGene/250k.txt";				
	my $k1M = -s "$myGene/1M.txt";
	my $sum  = eval join "+", ($k20k,$k40k, $k70k, $k100k, $k250k, $k1M);
	if ($sum == 0) {die "No interactions were found for your gene.\nYou either mistyped it, or no significant p < 0.05 interactions were found in that region.\n\n";}
}

sub writeinter{
	my $myGene = $_[0];
	my %chr_hash;
	# get genome
	#my $myGenome = $_[1];
	
	
	#opendir(DIR, ".");
	#my @sig_files = grep(/\.txt$/,readdir(DIR));
	#closedir(DIR);
	my @sig_files = ('20k.txt', '40k.txt', '70k.txt', '100k.txt', '250k.txt', '1M.txt');
	#print "@sig_files\n\n";
	
	foreach my $eachfile (@sig_files) {
		open (CIRC_INTS, ">$eachfile.circos.interactions.txt");
	 	open (SIG, "$eachfile");
		chomp (my @sig = (<SIG>));
		close(SIG);

	# set counter in bare block
	{
		my $interact_count = 0;
		foreach (@sig) {
			my @sig_fields = split(/,/, $_);	
			$interact_count++;
			my $thickness = $sig_fields[11];
			# regex to parse gene locations
			$sig_fields[0] =~ /(\w+):(\d+)-(\d+)/;
			# save coordiantes for location one
			my ($chr1, $start1, $end1) = ($1,$2,$3);
			# capture and save coordinates for location two
			$sig_fields[3] =~ /(\w+):(\d+)-(\d+)/;
			my ($chr2, $start2, $end2) = ($1,$2,$3);
			# fix mm vs chr problem
			$chr1 =~ s/chr/mm/g;
			$chr2 =~ s/chr/mm/g;
			# print locations and thicknesses to use with circos 
			print CIRC_INTS "interaction$interact_count\t$chr1\t$start1\t$end1\tthickness=$thickness\n";
			print CIRC_INTS "interaction$interact_count\t$chr2\t$start2\t$end2\tthickness=$thickness\n";		
			
			# chr hash stuff
			
			unless (exists $chr_hash{$chr1} ) {
				$chr_hash{$chr1} = ();
			}
		  unless (exists $chr_hash{$chr2} ) {
				$chr_hash{$chr2} = ();
			}
		}
	}
	
	$chrs = join(";", keys %chr_hash);
	$chrs .= ";";

	# fix chrs to mm
	$chrs =~ s/chr/mm/g;
	close(CIRC_INTS);
	# end writing circos file 
	
	}

}

sub writeconf{
	# open master circos file
	my $myGene = $_[0];
	# get genome
	my $myGenome = $_[1];
	open (BASE, " /NGS/Software/HiC/HiC_MASTER_circos.conf.txt");   ############################################## change to diectory on grid
	chomp (my @base = (<BASE>));
	close(BASE);
	
	# save stuff for cis writing
	# open 100k.txt file to extract chr and ranges
	open (SIG, "100k.txt");
	chomp (my @sig = (<SIG>));
	close(SIG);


	# variable for storing chr location for gene
	my $master_chr;
	# variable for storing start and end of coordinates on chr $master_chr
	my $master_start;
	my $master_end;

	# shift off first to capture chr 
	my $first_sig = shift @sig;
	# split first sig fields
	my @first_sig_fields = split(/,/, $first_sig);
	# check to see which anchor contains gene
	if ($first_sig_fields[2] =~ /$myGene/) {
		# split chr/coords from anchor 1
		$first_sig_fields[0] =~ /(\w+):(\d+)-(\d+)/;
		# save coordiantes for location one
		($master_chr, $master_start, $master_end) = ($1,$2,$3);
	} else {
		$first_sig_fields[3] =~ /(\w+):(\d+)-(\d+)/;
		# save coordiantes for location one
		($master_chr, $master_start, $master_end) = ($1,$2,$3);
	}

	foreach(@sig) {
			my @sig_fields = split(/,/, $_);
			# if anchor 1 contains gene
			if ($sig_fields[0] =~ /$master_chr/) {
				# split chr/coords from anchor 1
				$sig_fields[0] =~ /\w+:(\d+)-(\d+)/;
				# save coordiantes for location one, check if they are greater or less than previous
				# want to keep $1 if less than $master_start
				# want to keep $2 if greater than $master_end
				if ($1 < $master_start) {
					$master_start = $1;
				}
				if ($2 > $master_end) {
					$master_end = $2;
				}
			# else, anchor 2 must contain gene
			} 
			if ($sig_fields[3] =~ /$master_chr/) {
				$sig_fields[3] =~ /\w+:(\d+)-(\d+)/;
				# save coordiantes for location one
				if ($1 < $master_start) {
					$master_start = $1;
				}
				if ($2 > $master_end) {
					$master_end = $2;
				}
			}
		}



		# 
		# fix chr problem
		$master_chr =~ s/chr/mm/g;

		# calculate zoom lengths
		# divide by unit length
		$master_start = ($master_start / 250000);
		$master_end = ($master_end / 250000);
		# round start down with int
		$master_start = int($master_start);
		# round end up with int + 1
		$master_end = (int($master_end) + 1);
		# cat u to each
		$master_start = $master_start.'u';
		$master_end = $master_end.'u';




	
			# open TRANS circos file
			open (NEW_CONF_TRANS, ">$myGene.trans.circos.conf");
			foreach my $line (@base) {
				# if match regex for delim 1, change and print then next
				# ££change1	file = ./interChrom.circos.png ########## change for each file / image location
				if ($line =~ /^££change1/) {
					$line = "\tfile = ./$myGene.trans.circos.png";
					print NEW_CONF_TRANS "$line\n";
					next;
				}
				# if match regex for delim 2, change and print then next
				# ££change2 karyotype = interChrom.circos.karyotype.txt     #### change depending on genome 
				if ($line =~ /^££change2/) {
					# take organism name for karotype file
					#$line = "karyotype = /NGS/Software/HiC/karyotype.mouse.mm10.txt";    ############################################### change to directory on grid
					if ($myGenome eq "mm10") {
					$line = "karyotype = /NGS/Software/HiC/karyotype.mouse.mm10.txt";
					}
					if ($myGenome eq "mm9") {
					$line = "karyotype = /NGS/Software/HiC/karyotype.mouse.mm9.txt";
					}
					print NEW_CONF_TRANS "$line\n";
					next;
				}
				# if match regex for delim 3, change and print then next
			#	££change3 chromosomes = chr1;chr2;chr3;chr4;chr5;chr6;chr7;chr8;chr9;chr10;chr11;chr12;chr13;chr14;chr15;chr16;chr17;chr18;chr19;chr20;chr21;chr22;chrX;chrY;chrM    #### change depending on chroms involved
				if ($line =~ /^££change3/) {
					$line = "chromosomes = $chrs";
					print NEW_CONF_TRANS "$line\n";
					next;
				}
				# if match regex for delim 4, change and print then next
				# ££change4		file = interChrom.circos.interactions.txt ## change for each file, input interactions file 
				if ($line =~ /^££change4/) {
				# CHANGE TO 1OOK WHEN DONE WITH TESTS
					$line = "\t\tfile = 100k.txt.circos.interactions.txt";
					print NEW_CONF_TRANS "$line\n";
					next;
				}
				if ($line =~ /^££change5/) {
					$line = "\t\tfile = 1M.txt.circos.interactions.txt";
					print NEW_CONF_TRANS "$line\n";
					next;
				}
				if ($line =~ /^££change6/) {
					$line = "\t\tfile = 250k.txt.circos.interactions.txt";
					print NEW_CONF_TRANS "$line\n";
					next;
				}
				if ($line =~ /^££change7/) {
					$line = "\t\tfile = 20k.txt.circos.interactions.txt";
					print NEW_CONF_TRANS "$line\n";
					next;
				}
				if ($line =~ /^££change8/) {
					$line = "\t\tfile = 40k.txt.circos.interactions.txt";
					print NEW_CONF_TRANS "$line\n";
					next;
				}
				if ($line =~ /^££change9/) {
					$line = "\t\tfile = 70k.txt.circos.interactions.txt";
					print NEW_CONF_TRANS "$line\n";
					next;
				}
						if ($line =~ /^££changeBBB/) {
					$line = "";
					print NEW_CONF_TRANS "$line\n";
					next;
				}
				if ($line =~ /^££changeAAA/) {
					$line = "";
					print NEW_CONF_TRANS "$line\n";
					next;
				}
				if ($line =~ /^££changeCCC/) { # change file = /NGS/Software/HiC/rpkm_hist.txt
					if ($myGenome eq "mm10") {
					$line = "\t\tfile = /NGS/Software/HiC/rpkm_hist.txt";
					}
					if ($myGenome eq "mm9") {
					$line = "\t\tfile = /NGS/Software/HiC/rpkm_hist_mm9.txt";
					}
					print NEW_CONF_TRANS "$line\n";
					next;
				}				
				 else {
				 	print NEW_CONF_TRANS "$line\n";
				 }	
			}
			close(NEW_CONF_TRANS);
	
			# open CIS circos file
	
			open (NEW_CONF_CIS, ">$myGene.cis.circos.conf");
			open (BASE1, " /NGS/Software/HiC/HiC_MASTER_circos.conf.txt");   ############################################### change to directory on grid
			chomp (my @base1 = (<BASE1>));
			close(BASE1);
			foreach my $line1 (@base1) {
				# if match regex for delim 1, change and print then next
				# ££change1	file = ./interChrom.circos.png ########## change for each file / image location
				if ($line1 =~ /^££change1/) {
					$line1 = "\tfile = ./$myGene.cis.circos.png";
					print NEW_CONF_CIS "$line1\n";
					next;
				}
				# if match regex for delim 2, change and print then next
				# ££change2 karyotype = interChrom.circos.karyotype.txt     #### change depending on genome 
				if ($line1 =~ /^££change2/) {
					# take organism name for karotype file
					if ($myGenome eq "mm10") {
					$line1 = "karyotype = /NGS/Software/HiC/karyotype.mouse.mm10.txt";
					}
					if ($myGenome eq "mm9") {
					$line1 = "karyotype = /NGS/Software/HiC/karyotype.mouse.mm9.txt";
					}
					#$line1 = "karyotype = /NGS/Software/HiC/karyotype.mouse.mm10.txt";    ############################################### change to directory on grid
					print NEW_CONF_CIS "$line1\n";
					next;
				}
				# if match regex for delim 3, change and print then next
			#	££change3 chromosomes = chr1;chr2;chr3;chr4;chr5;chr6;chr7;chr8;chr9;chr10;chr11;chr12;chr13;chr14;chr15;chr16;chr17;chr18;chr19;chr20;chr21;chr22;chrX;chrY;chrM    #### change depending on chroms involved
				if ($line1 =~ /^££change3/) {
					$line1 = "chromosomes = $master_chr"; # changes
					print NEW_CONF_CIS "$line1\n";
					next;
				}
				# if match regex for delim 4, change and print then next
				# ££change4		file = interChrom.circos.interactions.txt ## change for each file, input interactions file 
				if ($line1 =~ /^££change4/) {
				# CHANGE TO 1OOK WHEN DONE WITH TESTS
					$line1 = "\t\tfile = 100k.txt.circos.interactions.txt";
					print NEW_CONF_CIS "$line1\n";
					next;
				}
				if ($line1 =~ /^££change5/) {
					$line1 = "\t\tfile = 1M.txt.circos.interactions.txt";
					print NEW_CONF_CIS "$line1\n";
					next;
				}
				if ($line1 =~ /^££change6/) {
					$line1 = "\t\tfile = 250k.txt.circos.interactions.txt";
					print NEW_CONF_CIS "$line1\n";
					next;
				}
				if ($line1 =~ /^££change7/) {
					$line1 = "\t\tfile = 20k.txt.circos.interactions.txt";
					print NEW_CONF_CIS "$line1\n";
					next;
				}
				if ($line1 =~ /^££change8/) {
					$line1 = "\t\tfile = 40k.txt.circos.interactions.txt";
					print NEW_CONF_CIS "$line1\n";
					next;
				}
				if ($line1 =~ /^££change9/) {
					$line1 = "\t\tfile = 70k.txt.circos.interactions.txt";
					print NEW_CONF_CIS "$line1\n";
					next;
				}
				if ($line1 =~ /^££changeBBB/) {
					$line1 = "";
					print NEW_CONF_CIS "$line1\n";
					print NEW_CONF_CIS "<plot>\n";
					print NEW_CONF_CIS "type             = text\n";
					print NEW_CONF_CIS "color            = red\n";
					if ($myGenome eq "mm10") {
					print NEW_CONF_CIS "file             = /NGS/Software/HiC/MM10_label2.txt\n";
					}
					if ($myGenome eq "mm9") {
					print NEW_CONF_CIS "file             = /NGS/Software/HiC/MM9_label2.txt\n";
					}
					#print NEW_CONF_CIS "file             = /NGS/Software/HiC/MM10_label2.txt\n";      ############################################### change to directory on grid
					print NEW_CONF_CIS "r0 = 1r\n";
					print NEW_CONF_CIS "r1 = 1.10r\n";
					print NEW_CONF_CIS "show_links     = yes\n";
					print NEW_CONF_CIS "link_dims      = 4p,4p,8p,4p,4p\n";
					print NEW_CONF_CIS "link_thickness = 2p\n";
					print NEW_CONF_CIS "link_color     = dred\n";
					print NEW_CONF_CIS "label_size   = 28p\n";
					print NEW_CONF_CIS "label_font   = condensed\n";
					print NEW_CONF_CIS "padding  = 2p\n";
					print NEW_CONF_CIS "rpadding = 2p\n";
					print NEW_CONF_CIS "</plot>\n";
					next;
				}
				if ($line1 =~ /^££changeAAA/) {
					$line1 = "";
					print NEW_CONF_CIS "$line1\n";
					print NEW_CONF_CIS "<zooms>\n";
					print NEW_CONF_CIS "<zoom>\n";
					print NEW_CONF_CIS "chr = $master_chr\n";   # changes
					print NEW_CONF_CIS "start = $master_start\n";   # changes
					print NEW_CONF_CIS "end   = $master_end\n";   # changes
		 	    print NEW_CONF_CIS "scale = 1000\n";
		 	    print NEW_CONF_CIS "</zoom>\n";
					print NEW_CONF_CIS "</zoom>\n";
					next;
				}
				if ($line1 =~ /^££changeCCC/) { # change file = /NGS/Software/HiC/rpkm_hist.txt
					if ($myGenome eq "mm10") {
					$line1 = "\t\tfile = /NGS/Software/HiC/rpkm_hist.txt";
					}
					if ($myGenome eq "mm9") {
					$line1 = "\t\tfile = /NGS/Software/HiC/rpkm_hist_mm9.txt";
					}
					print NEW_CONF_CIS "$line1\n";
					next;
				}
				 else {
				 	print NEW_CONF_CIS "$line1\n";
				 }	
			}
			close(NEW_CONF_CIS);	
	
	
}

sub make_diagrams {
	my $myGene = $_[0];
	system("/NGS/Software/circos-0.64/bin/circos -conf $myGene.trans.circos.conf 2>&1 | tee circos_trans_log.txt");
	system("/NGS/Software/circos-0.64/bin/circos -conf $myGene.cis.circos.conf 2>&1 | tee circos_cis_log.txt");
	system(q(find . -name "*circos.png" -type f -exec mv {} ./Graphs/ \;));
	system(q(find . -name "*circos.svg" -type f -exec mv {} ./Graphs/ \;));
}

sub write_table {

# write table of gene interactions in human readible form, with one copy of the bin containing gene
# make list of genes for each interaction file
# write R script for GO annotation with data frames printed
# execute R scripts

		# get gene name
		my $myGene = $_[0];
		# open directory (location is in gene name folder)
	#	opendir(DIR, ".");
		# grep gene interaction files
		#my @sig_files = grep(/\.txt$/,readdir(DIR));
		#closedir(DIR);
		#print "@sig_files\n\n";
		my @sig_files = ('20k.txt', '40k.txt', '70k.txt', '100k.txt', '250k.txt', '1M.txt');
		#print "@sig_list\n\n";

		# for each interaction bin size, going to make an R script and a better formatted table
		foreach my $eachfile (@sig_files) {
			open (GO_ANNO, ">$eachfile.GO_Anno.R");
		 	open (SIG, "$eachfile");
			chomp (my @sig = (<SIG>));
			close(SIG);
			open (TABLE, ">$eachfile.Interaction.Table.txt");
		
			# gene hash to keep gene names
			my %gene_hash;
			
			
			# shift @sig to extract 
			my $first_sig = shift @sig;
			
			
			# add gene to hash as below
			# get list of genes in each anchor (interaction partner)
			my @first_sig_fields = split(/,/, $first_sig);
			# store genes for each 	
			my @firstgene1 = split(/\s/, $first_sig_fields[2]);
			my @firstgene2 = split(/\s/, $first_sig_fields[5]);
			# add gene names from first anchor to hash
			foreach (@firstgene1) {
				unless (exists $gene_hash{$_}) {
					$gene_hash{$_} = ();
				}
			}
			# add gene names from second anchor to hash
			foreach (@firstgene2) {
				unless (exists $gene_hash{$_}) {
					$gene_hash{$_} = ();
				}
			}
			
			# for $first_sig_fields, extract $myGene bin info
			# if first anchor contains gene of interest
			if ($first_sig_fields[2] =~ /$myGene/) {
				print TABLE qq{
				Interactions involving $myGene:
	$myGene containing bin:
	$first_sig_fields[0]
	Entrez Gene IDs:
	$first_sig_fields[1]
	Genes Identified:
	$first_sig_fields[2]

				Information from Interacting Bins:

Location\tInteraction_Counts\tLog10-P-value\tFDR\tZ-score\tGene_Names\tGene_Entrez_IDs
$first_sig_fields[3]\t$first_sig_fields[6]\t$first_sig_fields[8]\t$first_sig_fields[9]\t$first_sig_fields[10]\t$first_sig_fields[5]\t$first_sig_fields[4]\n
}
			}
			# if second container contains gene of interest
			if ($first_sig_fields[5] =~ /$myGene/) {
				print TABLE qq{
				Interactions involving $myGene:
	$myGene containing bin:
	$first_sig_fields[3]
	Entrez Gene IDs:
	$first_sig_fields[4]
	Genes Identified:
	$first_sig_fields[5]

				Information from Interacting Bins:

Location\tInteraction_Counts\tLog10-P-value\tFDR\tZ-score\tGene_Names\tGene_Entrez_IDs
$first_sig_fields[0]\t$first_sig_fields[6]\t$first_sig_fields[8]\t$first_sig_fields[9]\t$first_sig_fields[10]\t$first_sig_fields[2]\t$first_sig_fields[1]\n
}
			}
			
			# for each line in a bin interaction/gene file
			foreach (@sig) {
				# chr6:30900000-30999999,116732 619665,Tsga13 Klf14,chr6:30500000-30599999,72649 75647 232680 71791,Tmem209 Ssmem1 Cpa2 Cpa4,483,0.0000000000,-109.292646,0.000000,2.191928,10
		
				# get list of genes in each anchor (interaction partner)
				my @sig_fields = split(/,/, $_);
				# store genes for each 	
				my @gene1 = split(/\s/, $sig_fields[2]);
				my @gene2 = split(/\s/, $sig_fields[5]);
				# write gene interaction details to Table
				# if first anchor contains gene of interest
				if ($sig_fields[2] =~ /$myGene/) {
					# print interaction data for second anchor and test results
					print TABLE "$sig_fields[3]\t$sig_fields[6]\t$sig_fields[8]\t$sig_fields[9]\t$sig_fields[10]\t$sig_fields[5]\t$sig_fields[4]\n";
				}
				# if second container contains gene of interest
				if ($sig_fields[5] =~ /$myGene/) {
					# print interaction data for first anchor and test results
					print TABLE "$sig_fields[0]\t$sig_fields[6]\t$sig_fields[8]\t$sig_fields[9]\t$sig_fields[10]\t$sig_fields[2]\t$sig_fields[1]\n";
				}
					
			# add gene names from first anchor to hash
				foreach (@gene1) {
					unless (exists $gene_hash{$_}) {
						$gene_hash{$_} = ();
					}
				}
				# add gene names from second anchor to hash
				foreach (@gene2) {
					unless (exists $gene_hash{$_}) {
						$gene_hash{$_} = ();
					}
				}
		
			}
			
		my $geneID  = join("\",\"", keys %gene_hash);
		$geneID .= "\"";
		$geneID = "\"".$geneID;

		
			
	
			# write R script to do g:Profile annotation
	
			my $Annostuff = qq{
				dir.create("./$eachfile.GO_Annotations")
				setwd("./$eachfile.GO_Annotations")
				library(gProfileR)
				# get list of genes
				genes <- c($geneID)
				# annotate with gProfiler using all available databases, max p-val of 005, fdr correction
				# none uses no hierarchical filtering, moderate uses moderate (best parent) hierarchical filtering, and strong uses strong (best parent group) hierarchical filtering
				# see gprofiler website for more info 
				none <- gprofiler(query = genes, organism = "mmusculus", max_p_value = 0.05, significant=TRUE, correction_method = "fdr", src_filter = c("GO", "KEGG", "REAC", "TF", "MI", "CORUM", "HP"), hier_filtering = "none")
				moderate <- gprofiler(query = genes, organism = "mmusculus", max_p_value = 0.05, significant=TRUE, correction_method = "fdr", src_filter = c("GO", "KEGG", "REAC", "TF", "MI", "CORUM", "HP"), hier_filtering = "moderate")
				strong <- gprofiler(query = genes, organism = "mmusculus", max_p_value = 0.05, significant=TRUE, correction_method = "fdr", src_filter = c("GO", "KEGG", "REAC", "TF", "MI", "CORUM", "HP"), hier_filtering = "strong")
				write.table(x = strong, file = "$eachfile.Gene_strong_GO.txt",sep = "\t")
				write.table(x = moderate, file = "$eachfile.Gene_moderate_GO.txt",sep = "\t")
				write.table(x = none, file = "$eachfile.Gene_none_GO.txt",sep = "\t")
				q("no")
			};
			
			print GO_ANNO $Annostuff;
			close(GO_ANNO);
			# execute R script to make functional annotation
			system(qq(chmod 777 $eachfile.GO_Anno.R));
			system(qq(/NGS/Software/HiC/R-patched/bin/R CMD BATCH $eachfile.GO_Anno.R));

	}
	
	system(q(find . -name "*.Table.txt" -type f -exec mv {} ./Tables/ \;));	

}
