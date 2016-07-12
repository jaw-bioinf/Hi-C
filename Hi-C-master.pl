#!/usr/bin/perl
#
# @author: John
# johnwilliams@post.harvard.edu
# j.williams@har.mrc.ac.uk
# HiC analysis master script

# to run: perl Hi-C-master.pl --project_dir directory --threads 10 --enzyme1 hindIII --fastq1 test_dataset1.fastq.gz --fastq2 test_dataset.fastq.gz --genome mm10
# to run with annotation only: perl HiC-master.pl --project_dir directory --threads 10 --enzyme1 hindIII --fastq test_dataset1.fastq.gz --fastq2 test_dataset.fastq.gz --genome mm10 --AnnotateOnly
# to run on grid, use master.sh and give the appropriate options. ie:
# 	bash master.pl --project_dir myDirectory --threads 10 --enzyme1 Hind3 --fastq1 myfastq_1.fastq.gz --fastq2 myfastq_2.fastq.gz --genome mm10
# to run with test data perl Hi-C-master.pl --project_dir test --threads 8 --enzyme1 hindIII --fastq1 test_dataset1.fastq.gz --fastq2 test_dataset2.fastq.gz --genome GRCh37
#
# if you encounter any problems, contact John. 
#


#
use strict; use warnings; use autodie;
use Getopt::Long;
use Cwd;


# get working directory
my $cwd = cwd();
open (TIMING, ">$cwd/timing.txt") or die("$!");
my $StartScriptTime = localtime;
print TIMING "Script started on $StartScriptTime\n";

# create parameter hash
my %hicup_param = (
help => '',
enzyme1 => '',
enzyme2 => '',
project_dir => '',
threads => '',
fastq1 => '',
fastq2 => '',
genome => '',
AnnotateOnly => ''
);

# use GetOptions to grab params passed by user
my %config_result = GetOptions(
    "help"        => \$hicup_param{help},
    "project_dir=s"    => \$hicup_param{project_dir},
    "enzyme1=s"       => \$hicup_param{enzyme1},
    "enzyme2=s"       => \$hicup_param{enzyme2},
    "threads=i"   => \$hicup_param{threads},
    "fastq1=s"    => \$hicup_param{fastq1},
    "fastq2=s"    => \$hicup_param{fastq2},
    "genome=s"    => \$hicup_param{genome},
    "AnnotateOnly"	=> \$hicup_param{AnnotateOnly}
);

# help module

sub helper { die(qq/
		This script may or may not be run in a directory containing fastq files. Fastq files can be gzipped, zipped, or uncompressed. 
		
		The scirpt will digest the genome, pre- and post-filter reads, align reads, create several hicup libraries with different bin sizes, and calculate significant reads for each bin size. Graphs from these stages will be included in the plots_tables directory.
		
		The script then takes interaction files and annotates them based on promoter regions. These annotated files, called <bin size>.annotation.txt, are located in subdirectories of the interactions directory. The annotation is done using the most recent UCSC RefGene info available via the GenomicFeatures and org.Mm.eg.db libraries.
		
		To extract interactions involving a gene of interest, run the separate gene_interactions.pl <still in development>.
		
		_______________________________________________
		
    USAGE : perl Hi-C-master <arguments> 
			ex: perl Hi-C-master --project_dir hic_project --threads 8 --enzyme1 hindIII --fastq1 test_dataset1.fastq.gz --fastq2 test_dataset2.fastq.gz --genome GRCh37
    ARGUMENTS :  
                    REQUIRED 
                    -- project_dir: the name of the project folder you would like to use
					-- enzyme1: the name of the enzyme used to digest, followed by sonication
					-- enzyme2: the name of the second enzyme if a double-digest and NOT sonication was used 
					-- threads: the number of threads used when applicable 
					-- fastq1: the name of one of the paired-end fastq files 
					-- fastq2: the name of the other paired-end fastq file 
					-- genome: the genome used for pre-processing and alignment, specify GRCh37, mm9, or mm10 
					-- help: display this message
					-- AnnotateOnly: Only run the interaction annotation, must be run after alignment, filtering and homer interactions have been made
					
					enzyme1 and enzyme2 can be of the following two forms:
						if HindIII or NCoI, they can be named:
						--enzyme1 HindIII
						else please give the enzyme and cut site:
						--enzyme1 A^AGCTT
						

               \n/); 
} 
 
if ( $hicup_param{help}) {
 &helper();
 } 

unless ($hicup_param{project_dir}) {
 print "\n MISSING ARGUMENTS : Give all the required options\n";
 &helper();
 } 
unless ($hicup_param{enzyme1}) { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &helper();} 
unless ($hicup_param{threads}) { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &helper();} 
unless ($hicup_param{fastq1}) { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &helper();} 
unless ($hicup_param{fastq2}) { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &helper();} 
unless ($hicup_param{genome}) { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &helper();} 



# pass relative project directory back to hash with full path 
my $prodir = $hicup_param{project_dir};
my $prodir1 = "$cwd/$prodir";
$hicup_param{project_dir} = "$prodir1";

# choose genome 
if ($hicup_param{genome} =~ "GRCh37") {
	$hicup_param{genome} = "/NGS/users/John/HIC-Roadmap/HICUP_TEST/reAlignRep1/sigInts100k_each/testing/GRCh37/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index";
	$hicup_param{org} = "GRCh37";
}

if ($hicup_param{genome} =~ "mm10") {
	$hicup_param{genome} = "/NGS/musRefs_10";
	$hicup_param{org} = "mm10";
}

# get genome mm9 location to do mm9 alignment

 if ($hicup_param{genome} =~ "mm9") {
 	$hicup_param{genome} = "/NGS/Software/HiC/mm9/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index";
 	$hicup_param{org} = "mm9";
 }

die "Could not parse options" unless (%config_result);

# check for program locations, store in hash, 
my $pathsref = &check_progs();



my %pathhash = %$pathsref;



  # check for directory
if (-d "$hicup_param{project_dir}") {
  } else {mkdir $hicup_param{project_dir};
  }
# cd to directory for rest of script exicutions
chdir("$hicup_param{project_dir}");


my $tohicup = \%hicup_param; 


 
if ( $hicup_param{AnnotateOnly}) {
&AnnotateOnly();
} else {
&RunFull();
}


sub RunFull {

	# run truncater step
	&hicup_trunc($pathsref,$tohicup);
	# run mapping step
	&hicup_map($pathsref,$tohicup);
	# run digester step
	&hicup_digest($pathsref,$tohicup);
	# run filter step
	&hicup_filter($pathsref,$tohicup);
	# run deduplication step
	&hicup_deduplicator($pathsref,$tohicup);
	# run hicup to homer conversion
	&hicup_2_homer($pathsref,$tohicup);
	# make homer tags
	&homer_tags($pathsref,$tohicup);
	# make homer background Hi-C models and compute interaction frequencies
	&homer_interactions($pathsref,$tohicup); 
	my $orgtime = localtime;
	print TIMING "File organizaiton step started on $orgtime, next step is end of script.\n";

	# chmod project directory 
	system("chmod -R 777 $hicup_param{project_dir}");

	# move svg files and summary text files to plots_tables directory
	chdir ("$hicup_param{project_dir}");
	mkdir("./plots_tables");
	system(q{find . -name "*.svg" -type f -exec mv {} ./plots_tables/ \;});
	system(q{find . -name "*summary*" -type f -exec mv {} ./plots_tables/ \;});

	# move FASTQ, SAM, and BAM files to common directory
	mkdir("./BIG_FILES");
	system(q{find . -name "*.sam" -type f -exec mv {} ./BIG_FILES/ \;});
	system(q{find . -name "*.bam" -type f -exec mv {} ./BIG_FILES/ \;});
	system(q{find . -name "*.fastq.gz" -type f -exec mv {} ./BIG_FILES/ \;});
	system(q{find . -name "*.sam.homer" -type f -exec mv {} ./BIG_FILES/ \;});

	# annotate interactions and run futher stats/plots
	&genomicinteractions_annotations($pathsref,$tohicup); 

	my $EndScriptTime = localtime;
	print TIMING "Script Ended on $EndScriptTime\n";

	close(TIMING);
	
	# move files to clean up things 
	# move timing file 
	chdir ("$hicup_param{project_dir}");
	chdir ("..");
	rename "timing.txt", "$hicup_param{project_dir}/timing.txt";


  # annotate interactions and run futher stats/plots
	&genomicinteractions_annotations($pathsref,$tohicup); 


}



# if --AnnotateOnly selected, skip all the above and only do annotate interactions 


sub AnnotateOnly {

  # annotate interactions and run futher stats/plots
	&genomicinteractions_annotations($pathsref,$tohicup); 


}
  
sub check_progs {
  my $start = localtime;
 # print TIMING "check_progs() started at $start\n";

       my %paths;
       $paths{samtools} = "/NGS/Software/samtools-0.1.19/samtools";
       $paths{R} = "/NGS/Software/HiC/R-patched/bin/R";
       $paths{bowtie2} = "/NGS/Software/bowtie/bowtie2";
       $paths{homer} = "/NGS/Software/Homer/bin/homer";
       $paths{hicup} = "/NGS/users/John/HIC-Roadmap/HICUP_TEST/hicup_v0.5.7/hicup";
       $paths{circos} = "/NGS/Software/circos-0.64/bin/circos";
		   $paths{hicup2homer} = "/NGS/users/John/HIC-Roadmap/HICUP_TEST/hicup_v0.5.7/Conversion/hicup2homer";
			 $paths{hicup_truncater} = "/NGS/users/John/HIC-Roadmap/HICUP_TEST/hicup_v0.5.7/hicup_truncater";
       $paths{hicup_mapper} = "/NGS/users/John/HIC-Roadmap/HICUP_TEST/hicup_v0.5.7/hicup_mapper";
       $paths{hicup_digester} = "/NGS/users/John/HIC-Roadmap/HICUP_TEST/hicup_v0.5.7/hicup_digester";
       $paths{hicup_deduplicator} = "/NGS/users/John/HIC-Roadmap/HICUP_TEST/hicup_v0.5.7/hicup_deduplicator";
       $paths{hicup_filter} = "/NGS/users/John/HIC-Roadmap/HICUP_TEST/hicup_v0.5.7/hicup_filter";
		   $paths{makeTagDirectory} = "/NGS/Software/Homer/bin/makeTagDirectory";
		   $paths{analyzeHiC} = "/NGS/Software/Homer/bin/analyzeHiC";

		
      #print Dumper(\%paths);

        
      my $pathsref = \%paths; 
          my $end = localtime;
  #  print TIMING "check_progs() finished at $end\n";
 #   print TIMING "**************************************************\n"; 

      return $pathsref
}

sub hicup_trunc {
  my $start = localtime;
  print TIMING "hicup_trunc() started at $start\n";

  my $pathvars = $_[0];
  my %pathhash = %$pathvars;
  
  my $configvars = $_[1];
  my %confighash = %$configvars;

  
  # begin truncation step. since writing a config file for hicup seems to work better than directly passing arguments, and gives a concrete record of arguments passed for future comparisons, we'll write a config file
  open(TRUNC,">trunc.conf") or die("$!");
  print TRUNC "#Directory for output\n";
  print TRUNC "Outdir:$confighash{project_dir}\n";
  print TRUNC "# Number of threads\n";
  print TRUNC "Threads:$confighash{threads}\n";
  print TRUNC "#Suppress progress updates (0: off, 1: on)\n";
  print TRUNC "Quiet:0\n";
  print TRUNC "#Zip output (0: off, 1: on)\n";
  print TRUNC "Zip:1\n";
  
  print TRUNC "#Fill-in of sticky ends was not performed prior to ligation\n#Only select if fill-in was NOT performed\nNoFill:0\n";
  print TRUNC "#Restriction site that generates the Hi-C junction.\n#Insert a caret ('^') to represent the position where the restriction enzyme cuts the DNA strand\n";
  die("You must give enzyme1. Make sure it is the correct enzyme used in your HiC protocol.") unless $confighash{enzyme1};
  # make it so enzymes can be given as HindIII, NcoI, or as a sequence with ^ cut place
  if ($confighash{enzyme1} =~ /HindIII|Hind3/i) {
    print TRUNC "re1:A^AGCTT\n";
  } elsif ($confighash{enzyme1} =~ /NcoI|Nco1/i) {
    print TRUNC "re1:C^CATGG\n";
  } else { print TRUNC "re1:$confighash{enzyme1}\n"};
  if ($confighash{enzyme2}) {
    if ($confighash{enzyme2} =~ /HindIII|Hind3/i) {
      print TRUNC "re2:A^AGCTT\n";
    } elsif ($confighash{enzyme2} =~ /NcoI|Nco1/i) {
      print TRUNC "re2:C^CATGG\n";
    } else {print TRUNC "re2:$confighash{enzyme2}\n"};
  }
  print TRUNC "#FASTQ sequence files, can be .zip, .gz, or uncompressed\n";
  
  print TRUNC "$cwd/$confighash{fastq1}\n";
  print TRUNC "$cwd/$confighash{fastq2}\n";

  close(TRUNC);

  # store program for hicup truncater
  #my $hicup_trunc = "$pathhash{hicup_truncater} --config trunc.conf";
  system (qq($pathhash{hicup_truncater} --config trunc.conf));
  # print done
  print "HiCUP Truncation done! Preparing for mapping steps.\n";
  
  my $end = localtime;
  print TIMING "hicup_trunc() finished at $end\n";
  print TIMING "**************************************************\n"; 

}


sub hicup_map {
  my $start = localtime;
  print TIMING "hicup_map() started at $start\n";

  my $pathvars = $_[0];
  my %pathhash = %$pathvars;
  
  my $configvars = $_[1];
  my %confighash = %$configvars;

  # store current (project) directory
  my $wd = cwd();
  # find truncated fastq files
  opendir(DIR, ".");
  my @files = grep(/trunc\.fastq/,readdir(DIR));
  closedir(DIR);


  
  
  # begin mapping step.
  # make output directory
  mkdir "$wd/map";
  open(MAP,">map.conf") or die("$!");
  print MAP "#Directory for output\n";
  print MAP "Outdir:$wd/map\n";
  print MAP "# Number of threads\n";
  print MAP "Threads:$confighash{threads}\n";
  print MAP "#Suppress progress updates (0: off, 1: on)\n";
  print MAP "Quiet:0\n";
  print MAP "#Zip output (0: off, 1: on)\n";
  print MAP "Zip:1\n";  
  print MAP "#Path to the alignment program Bowtie (include the executable Bowtie filename)\n";
  if ($pathhash{bowtie2}) {
          print MAP "Bowtie2:$pathhash{bowtie2}\n";} else {
          print MAP "Bowtie:$pathhash{bowtie}\n";
        }

  print MAP "#Path to the reference genome indices\n#This should include the basename of the index:\n#[path to folder]/[name of folder containing indices]/[index basename]\n";
  # choose correct directory 
  if ($confighash{org} =~ "GRCh37") {
  print MAP "Index:$confighash{genome}/genome\n";
  }
  
  if ($confighash{org} =~ "mm10") {
  print MAP "Index:$confighash{genome}/ref_mm10\n";
  }
  
  if ($confighash{org} =~ "mm9") {
  print MAP "Index:$confighash{genome}/genome\n";
  }
  
  print MAP "#FASTQ sequence files, can be .zip, .gz, or uncompressed\n";
  foreach my $file (@files) {
    print MAP "$wd/$file\n";
  }
  

  close(MAP);

  # store program for hicup mapper
 # my $hicup_map = "$pathhash{hicup_mapper} --config map.conf";
  system (qq($pathhash{hicup_mapper} --config map.conf));
    print "HiCUP Mapping done! Preparing for digesting steps.\n";
      my $end = localtime;
  print TIMING "hicup_map() finished at $end\n";
  print TIMING "**************************************************\n"; 

}

sub hicup_digest {
  my $start = localtime;
  print TIMING "hicup_digest() started at $start\n";

  my $pathvars = $_[0];
  my %pathhash = %$pathvars;
  
  my $configvars = $_[1];
  my %confighash = %$configvars;

  # store current (project) directory
  my $wd = cwd();
  # find .fa genome files
  opendir(DIR, "$confighash{genome}");
  my @files = grep(/\.fa/,readdir(DIR));
  closedir(DIR);
  
  
  # begin digestion step.
  # make output directory
  mkdir "$wd/digest";
  open(DIGEST,">digest.conf") or die("$!");
  print DIGEST "#Directory for output\n";
  print DIGEST "Outdir:$wd/digest\n";
  print DIGEST "#Zip output (0: off, 1: on)\n";
  print DIGEST "Zip:1\n";  
  print DIGEST "#Restriction site that generates the Hi-C junction.\n#Insert a caret ('^') to represent the position where the restriction enzyme cuts the DNA strand\n";
  die("You must give enzyme1. Make sure it is the correct enzyme used in your HiC protocol.") unless $confighash{enzyme1};
  # make it so enzymes can be given as HindIII, NcoI, or as a sequence with ^ cut place
  if ($confighash{enzyme1} =~ /HindIII|Hind3/i) {
    if ($confighash{enzyme2} =~ /NcoI|Nco1/i) {
      print DIGEST "re1:A&AGCTT:C^CATGG\n";
    } else{ print DIGEST "re1:A^AGCTT\n";}
  } elsif ($confighash{enzyme1} =~ /NcoI|Nco1/i) {
    print DIGEST "re1:C^CATGG\n";
  } else { print DIGEST "re1:$confighash{enzyme1}\n"};
  
  
  print DIGEST "#Name of the genome to be digested\nGenome:Genome\n";
  
  # create dummy file 
  open (file1, ">./digest/file.fa");
  print file1 ">test\naaa";
  close(file1);
  
    

  # if genome is GRCh37
  
  if ($confighash{org} =~ "GRCh37") {
	  print DIGEST "#Sequence files to be digested, must end in .fa\n";
	  foreach my $file (@files) {
		print DIGEST "$confighash{genome}/$file\n";
		print DIGEST "$wd/digest/file.fa"
	  }
  }
  
  # if genome is mm10
  if ($confighash{org} =~ "mm10") {
	  print DIGEST "#Sequence files to be digested, must end in .fa\n";
		print DIGEST "$confighash{genome}/ref_mm10.fa\n";
		print DIGEST "$wd/digest/file.fa"
  }
  
  #if genome is mm9
   if ($confighash{org} =~ "mm9") {
	   print DIGEST "#Sequence files to be digested, must end in .fa\n";
		foreach my $file (@files) {
			print DIGEST "$confighash{genome}/$file\n";
			print DIGEST "$wd/digest/file.fa"
	  }
   }
  
  
  close(DIGEST);

  # store program for hicup digester
  #my $hicup_dig = "$pathhash{hicup_digester} --config digest.conf";
  system (qq($pathhash{hicup_digester} --config digest.conf));
    print "HiCUP Genome Digesting done! Preparing for filtering steps.\n";
    
  my $end = localtime;
  print TIMING "hicup_digest() finished at $end\n";
  print TIMING "**************************************************\n"; 

}

sub hicup_filter {
  my $start = localtime;
  print TIMING "hicup_filter() started at $start\n";

  my $pathvars = $_[0];
  my %pathhash = %$pathvars;
  
  my $configvars = $_[1];
  my %confighash = %$configvars;

  # store current (project) directory
  my $wd = cwd();

  # find digest fastq file
  opendir(DIR1, "$wd/digest") or die("$!");
  my @files = grep(/^Digest/,readdir(DIR1));
  closedir(DIR1);
  
  # find paired files from mapping output
  opendir(DIR2, "$wd/map") or die("$!");
  my @files1 = grep(/\.pair/,readdir(DIR2));
  closedir(DIR2);
  
  
  # begin filter step.
  # make output directory
  mkdir "$wd/filter";
  open(FILTER,">filter.conf") or die("$!");
  print FILTER "#Directory for output\n";
  print FILTER "Outdir:$wd/digest\n";
  print FILTER "##Reference genome generated by the 'hicup_digester' perl script\n";
 foreach my $file (@files) {
    print FILTER "Digest:$wd/digest/$file\n";

  }
  
  print FILTER "# Number of threads\n";
  print FILTER "Threads:$confighash{threads}\n";
  print FILTER "#Suppress progress updates (0: off, 1: on)\n";
  print FILTER "Quiet:0\n";
  # print long and short thresholds, these settings should be good for most uses
  print FILTER "#Minimum allowable di-tag size (optional)\nShortest: 100\n#Maximum allowable di-tag size (optional)\nLongest:  800\n";
  
  
  
  
  # make dummy file to get around hicup bug
  open (file3, ">./filter/dummy.sam");
  print file3 ">test\naaa";
  close(file3);
  
  print FILTER "#List '.pair' files to be categorised\n";
  foreach my $file1 (@files1) {
    print FILTER "$wd/map/$file1\n";
    print FILTER "$wd/filter/dummy.sam\n";
  }
  

  close(FILTER);

  # store program for hicup filter
 # my $hicup_filter = "$pathhash{hicup_filter} --config filter.conf";
  system (qq($pathhash{hicup_filter} --config filter.conf));
    print "HiCUP Filtering done! Preparing for deduplication in next step.\n";
  my $end = localtime;
  print TIMING "hicup_filter() finished at $end\n";
  print TIMING "**************************************************\n"; 

}

sub hicup_deduplicator {
  my $start = localtime;
  print TIMING "hicup_deduplicator() started at $start\n";

  my $pathvars = $_[0];
  my %pathhash = %$pathvars;
  
  my $configvars = $_[1];
  my %confighash = %$configvars;
  
  my $wd1 = cwd();

  
  # delete stupid duplicate file
  unlink glob "./digest/dummy.filt.sam";
  
  # find filtered sam/bam file
  opendir(DIRE, "$wd1/digest") or die("$!");
  my @filez = grep(/\.filt\./,readdir(DIRE));
  closedir(DIRE);
  
  
  # write deduplication config file
  open(DEDUP,">dedup.conf") or die("$!");
  print DEDUP "#Directory for output\n";
  print DEDUP "Outdir:$wd1\n";
 
  print DEDUP "# Number of threads\n";
  print DEDUP "Threads:$confighash{threads}\n";
  print DEDUP "#Suppress progress updates (0: off, 1: on)\n";
  print DEDUP "Quiet:0\n";

  
  
  # make dummy file to get around hicup bug
  open (file3, ">./dummy1.filt.sam");
  print file3 ">test\naaa";
  close(file3);
  
  print DEDUP "#List '.SAM/BAM' files to be de-duplicated\n";
  print DEDUP "$wd1/dummy1.filt.sam\n";
  print DEDUP "$wd1/digest/@filez\n";
  
  close(DEDUP);

  
 # store program for hicup filter
  #my $hicup_dedup = "$pathhash{hicup_deduplicator} --config dedup.conf";
  system (qq($pathhash{hicup_deduplicator} --config dedup.conf));
    print "HiCUP Deduplication done! Preparing for building Homer-style HiC summary table in next step.\n";

	
  my $end = localtime;
  print TIMING "hicup_deduplicator() finished at $end\n";
  print TIMING "**************************************************\n"; 
}

sub hicup_2_homer {
  my $start = localtime;
  print TIMING "hicup_2_homer() started at $start\n";

  my $pathvars = $_[0];
  my %pathhash = %$pathvars;
  
  my $configvars = $_[1];
  my %confighash = %$configvars;
  
   my $wd2= cwd();

  
  #delete stupid duplicate file
   unlink glob "./dummy1.*";
  
  #find filtered sam/bam file
   opendir(DIRE, "$wd2") or die("$!");
   my @filez = grep(/\.dedup\./,readdir(DIRE));
   closedir(DIRE);
  
 
  # store program for homer summary file generation
  #my $hicup_2_hom = "$pathhash{hicup2homer} $wd2/@filez";
  system (qq($pathhash{hicup2homer} $wd2/@filez));
    print "Homer file generation complete, ready to send interactions to Homer! Good bye HiCUP, and thanks for the help!.\n\n";
  my $end = localtime;
  print TIMING "hicup_2_homer() finished at $end\n";
  print TIMING "**************************************************\n"; 

}

sub homer_tags {
  my $start = localtime;
  print TIMING "homer_tags() started at $start\n";

  my $pathvars = $_[0];
  my %pathhash = %$pathvars;
  
  my $configvars = $_[1];
  my %confighash = %$configvars;

  # create directory for homer tags and move there
  my $twd = cwd();
  if (-d "$twd/tags") {
    } else {mkdir "$twd/tags" or die("$!");
    }
  chdir("$twd/tags");
  # find hic-summary file created by HiCUP
  opendir(DIRE, "$twd") or die("$!");
  my @filez = grep(/\.homer$/,readdir(DIRE));
  closedir(DIRE);
  
  # store program for tag directory generation
  #my $tag_maker = "$pathhash{makeTagDirectory} $twd/tags -format HiCsummary $twd/@filez";

  system(qq($pathhash{makeTagDirectory} $twd/tags -format HiCsummary $twd/@filez));
  print "Tag directory made!\n";
    my $end = localtime;
  print TIMING "homer_tags() finished at $end\n";
  print TIMING "**************************************************\n"; 

}

sub homer_interactions {
  my $start = localtime;
  print TIMING "homer_interactions() started at $start\n";

  my $pathvars = $_[0];
  my %pathhash = %$pathvars;
  
  my $configvars = $_[1];
  my %confighash = %$configvars;
  
  my $awd = cwd();
  print "$awd\n\n";
  
  # make folder for storing interaction files & cd there 
  if (-d "../interactions") {
    } else {mkdir "../interactions" or die("$!");
    }
  chdir("../interactions");
  
  # store programs to create background models and analyze interactions. Background models will be created automatically unless they alredy exist.
 # my $bin100k = "$pathhash{analyzeHiC} $awd -res 100000 -interactions sigints100k.txt -nomatrix -cpu $confighash{threads} -pvalue 0.05 2>&1 | tee 100k_log.txt";
 # my $bin250k = "$pathhash{analyzeHiC} $awd -res 250000 -interactions sigints250k.txt -nomatrix -cpu $confighash{threads} -pvalue 0.05 2>&1 | tee 250k_log.txt";
 # my $bin1M   = "$pathhash{analyzeHiC} $awd -res 1000000 -interactions sigints1M.txt -nomatrix -cpu $confighash{threads} -pvalue 0.05 2>&1 | tee 1M_log.txt";
 # my $bin40k  = "$pathhash{analyzeHiC} $awd -res 40000 -interactions sigints40k.txt -nomatrix -cpu $confighash{threads} -pvalue 0.05 2>&1 | tee 40k_log.txt";
 # my $bin70k  = "$pathhash{analyzeHiC} $awd -res 70000 -interactions sigints70k.txt -nomatrix -cpu $confighash{threads} -pvalue 0.05 2>&1 | tee 70k_log.txt";
 # my $bin20k  = "$pathhash{analyzeHiC} $awd -res 20000 -interactions sigints20k.txt -nomatrix -cpu $confighash{threads} -pvalue 0.05 2>&1 | tee 20k_log.txt";
  system(qq($pathhash{analyzeHiC} $awd -res 100000 -interactions sigints100k.txt -nomatrix -cpu $confighash{threads} -pvalue 0.05 2>&1 | tee 100k_log.txt));
  system(qq($pathhash{analyzeHiC} $awd -res 250000 -interactions sigints250k.txt -nomatrix -cpu $confighash{threads} -pvalue 0.05 2>&1 | tee 250k_log.txt));
  system(qq($pathhash{analyzeHiC} $awd -res 1000000 -interactions sigints1M.txt -nomatrix -cpu $confighash{threads} -pvalue 0.05 2>&1 | tee 1M_log.txt));
  system(qq($pathhash{analyzeHiC} $awd -res 40000 -interactions sigints40k.txt -nomatrix -cpu $confighash{threads} -pvalue 0.05 2>&1 | tee 40k_log.txt));
  system(qq($pathhash{analyzeHiC} $awd -res 70000 -interactions sigints70k.txt -nomatrix -cpu $confighash{threads} -pvalue 0.05 2>&1 | tee 70k_log.txt));
  system(qq($pathhash{analyzeHiC} $awd -res 20000 -interactions sigints20k.txt -nomatrix -cpu $confighash{threads} -pvalue 0.05 2>&1 | tee 20k_log.txt));

  print "HiC interactions analyzed, no matrices generated.\n";
  print "Changing FDR notation to make R stats easier.\n";
  
  # find sig interaction files file
	opendir(DIRa, ".");
	my @fileza = grep(/^sigints.*.txt/,readdir(DIRa));
	closedir(DIRa);
	# add fdr column to sigints files 
	foreach my $file (@fileza) {
		open (OLD, $file);
		chomp (my @OLD = (<OLD>));
		my $new = "new.txt";
		open (NEW, ">$new");
	
		# copy fdr column with new name simplified to fdr (without unkown number of tests mentioned)
		foreach my $line (@OLD) {
			my @fields = split(/\t/, $line);
			my $newcol = $fields[18];
			push(@fields, $newcol);
			if ($line =~ /^InteractionID/) {
				$fields[18] = "fdr";
			}
			$line = join ("\t", @fields);
			print NEW "$line\n";
		}
		# end doing stuff
			close (OLD);
			close (NEW);
			rename ($file, "$file.orig");
			rename ($new, $file);
	}
  
  # end timing 
  my $end = localtime;
  print TIMING "homer_interactions() finished at $end\n";
  print TIMING "**************************************************\n"; 

}

sub genomicinteractions_annotations {
	my $pathvars = $_[0];
  my %pathhash = %$pathvars;
  
  my $configvars = $_[1];
  my %confighash = %$configvars;

	chdir ("$hicup_param{project_dir}");
	chdir ("./interactions");
	my $lwd = cwd();

	my @files = ("sigints20k", "sigints40k", "sigints70k", "sigints250k", "sigints100k", "sigints1M");

	#print @files;
	
	

	foreach  my $sigIntsFile (@files) {
		open (SCRIPT, ">$sigIntsFile.RSCRIPT.R");
		
		# write if loop for genome == mm10 here surrounding my $Rthing 
  if ($confighash{org} =~ "mm10") {
		my $Rthing = q{
		mytime <- Sys.time()
		cat(sprintf("Annotation of £££ started on: %s\n", mytime),file="./annotation_timing.txt", append=TRUE)

		library(GenomicFeatures)
		mm10.refseq.db <- makeTxDbFromUCSC(genome="mm10", table="refGene")
		# get genes and transcripts out, keep transcripts with gene names
		refseq.genes = genes(mm10.refseq.db)
		refseq.transcripts = transcriptsBy(mm10.refseq.db, by="gene")
		refseq.transcripts = refseq.transcripts[ names(refseq.transcripts) %in% unlist(refseq.genes$gene_id) ]
		# get promoter IDs, upstream 1.5kb, downstream 1.5kb ####### this can change
		mm10_refseq_promoters <- promoters(refseq.transcripts, upstream = 1500, downstream = 1500)
		# annotate with promoter ids
		mm10_annotation.features <- list(promoter = mm10_refseq_promoters)

		# create library for converting gene entrez IDs to gene symbols
		library(org.Mm.eg.db)
		x <- org.Mm.egSYMBOL
		mapped_genes <- mappedkeys(x)
		xx <- as.list(x[mapped_genes])
		 
		library(GenomicInteractions)
		library(GenomicRanges)
		# read data in from summary output from Homer's analysis

		sigD <- makeGenomicInteractionsFromFile("£££.txt", type="homer", experiment_name = "myHi-Cexpt", description="myHi-Cexpt")
		 
		# make p-values for each
		sigD$p.value <- exp(sigD$LogP)
	
		# make subset with fdr <=0.05 for cutoff
		sigD_subset <- sigD[sigD$fdr <= 0.05]
	
		# annotate FDR subset interactions
			annotateInteractions(sigD_subset, mm10_annotation.features)


		 # make directory for output charts and tables
		 dir.create("./£££")
		 setwd("./£££/")
		 # create file for logging sigD
		 sigDlog <- file("£££.log.txt", open="at")
		 # open connection to sigD log
		 sink(sigDlog)
		 # write what doing
		 cat("Anchor One width -- £££\n")
		 # check bins for each anchor for each expt
		 summary(width(anchorOne(sigD)))
		 cat("Anchor Two width -- £££\n")
		 # check bins for each anchor for each expt
		 summary(width(anchorTwo(sigD)))
		 # check counts
		 cat("mean interaction counts\n")
		 mean(interactionCounts(sigD))
		 # number of sig interactions
		 cat("number of interactions q < 0.05\n")
		 sum(sigD$fdr < 0.05)
		 # close connection to sigD log
		 sink()
		 # plot densities for FDR and p-vals
		 svg("£££.p-valuePlot.svg")
		 plot(density(sigD$p.value))
		 dev.off()
		 svg("£££.fdrPlot.svg")
		 plot(density(sigD$fdr))
		 dev.off()
		 # get number of fdr < 0.05
		 #sum(hg18_data$fdr < 0.05)
		 # plot interaction types
		 svg("£££.interactionTypes.svg")
		 plotInteractionAnnotations(sigD_subset)
		 dev.off()
			# plot cis and trans percentage
		 svg("£££.cisTransAll.svg")
		 plotCisTrans(sigD)
		 dev.off()
		 svg("£££.cisTransStrict.svg")
		 plotCisTrans(sigD_subset)
		 dev.off()
		 # plot read counts to check if subsetting by strict FDR drops low counts
		 svg("£££.countsAll.svg")
		 plotCounts(sigD, cut=50)
		 dev.off()
		 svg("sigD_countsStrict.svg")
		 plotCounts(sigD_subset, cut=50)
		 dev.off()
		 
		cat("location1,geneIDs1,geneSymbols1,location2,geneIDs2,geneSymbols2,counts,p-value,logP-val,fdr,z-score,circos-thickness\n", file="£££.annotation.txt", append=TRUE)
		my_annot_function <- function() {
			for (i in 1:length(sigD_subset)) {
				# store current chr number
				chrN <- as.character(sigD_subset[i,]@anchor_one@seqnames@values)
				# store second chr number
				chrN2 <- as.character(sigD_subset[i,]@anchor_two@seqnames@values)
				# access anchor 1 start location
				s1 <- sigD_subset[i,]@anchor_one@ranges@start
				# access anchor 1 end location
				e1 <- sigD_subset[i,]@anchor_one@ranges@width-1 + sigD_subset[i,]@anchor_one@ranges@start
				# access anchor 2 start
				s2 <- sigD_subset[i,]@anchor_two@ranges@start
				# access anchor 2 end
				e2 <- sigD_subset[i,]@anchor_two@ranges@width-1 + sigD_subset[i,]@anchor_two@ranges@start
				# access Entrez Gene IDs for anchor 1
				gID1 <- unlist(sigD_subset[i,]@anchor_one$promoter.id)
				# access Entrez Gene IDs for anchor 2
				gID2 <- unlist(sigD_subset[i,]@anchor_two$promoter.id)
				# access counts (number of interactions)
				count1 <- sigD_subset[i,]@counts
				# access p-val
				pv1 <-  sigD_subset[i,]$p.value
				#access logP value
				lpv1 <- sigD_subset[i,]$LogP
				# access fdr
				fdr1 <- sigD_subset[i,]$fdr
				# access Z-score
				zscore1 <- sigD_subset[i,]$Z.score
				# access circos thickness
				circ1 <- sigD_subset[i,]$Circos.Thickness
				# map gene entrez IDs to gene symbols
				name1 <- unlist(xx[gID1])
				name2 <- unlist(xx[gID2])
				cat(sprintf("%s:%s-%s,%s,%s,%s:%s-%s,%s,%s,%d,%.10f,%s,%f,%f,%d\n", chrN, s1, e1, paste(gID1, collapse=" "), paste(name1, collapse=" "), chrN2, s2, e2, paste(gID2, collapse=" "),  paste(name2, collapse=" "),count1, pv1, lpv1, fdr1, zscore1, circ1), file="£££.annotation.txt", append=TRUE)
			}
		}
		my_annot_function()
		mytime2 <- Sys.time()
		cat(sprintf("Annotation of £££ finished on: %s\n", mytime2),file="../annotation_timing.txt", append=TRUE)
		q("no")
		};
		
				print SCRIPT $Rthing;

		
	};


		# copy Rthing code above, wrap in if mm9, change mmm10 references to mm9 
		# for mm9 annotation 
  if ($confighash{org} =~ "mm9") {
		my $Rthing = q{
		mytime <- Sys.time()
		cat(sprintf("Annotation of £££ started on: %s\n", mytime),file="./annotation_timing.txt", append=TRUE)

		library(GenomicFeatures)
		mm9.refseq.db <- makeTxDbFromUCSC(genome="mm9", table="refGene")
		# get genes and transcripts out, keep transcripts with gene names
		refseq.genes = genes(mm9.refseq.db)
		refseq.transcripts = transcriptsBy(mm9.refseq.db, by="gene")
		refseq.transcripts = refseq.transcripts[ names(refseq.transcripts) %in% unlist(refseq.genes$gene_id) ]
		# get promoter IDs, upstream 1.5kb, downstream 1.5kb ####### this can change
		mm9_refseq_promoters <- promoters(refseq.transcripts, upstream = 1500, downstream = 1500)
		# annotate with promoter ids
		mm9_annotation.features <- list(promoter = mm9_refseq_promoters)

		# create library for converting gene entrez IDs to gene symbols
		library(org.Mm.eg.db)
		x <- org.Mm.egSYMBOL
		mapped_genes <- mappedkeys(x)
		xx <- as.list(x[mapped_genes])
		 
		library(GenomicInteractions)
		library(GenomicRanges)
		# read data in from summary output from Homer's analysis

		sigD <- makeGenomicInteractionsFromFile("£££.txt", type="homer", experiment_name = "myHi-Cexpt", description="myHi-Cexpt")
		 
		# make p-values for each
		sigD$p.value <- exp(sigD$LogP)
	
		# make subset with fdr <=0.05 for cutoff
		sigD_subset <- sigD[sigD$fdr <= 0.05]
	
		# annotate FDR subset interactions
			annotateInteractions(sigD_subset, mm9_annotation.features)


		 # make directory for output charts and tables
		 dir.create("./£££")
		 setwd("./£££/")
		 # create file for logging sigD
		 sigDlog <- file("£££.log.txt", open="at")
		 # open connection to sigD log
		 sink(sigDlog)
		 # write what doing
		 cat("Anchor One width -- £££\n")
		 # check bins for each anchor for each expt
		 summary(width(anchorOne(sigD)))
		 cat("Anchor Two width -- £££\n")
		 # check bins for each anchor for each expt
		 summary(width(anchorTwo(sigD)))
		 # check counts
		 cat("mean interaction counts\n")
		 mean(interactionCounts(sigD))
		 # number of sig interactions
		 cat("number of interactions q < 0.05\n")
		 sum(sigD$fdr < 0.05)
		 # close connection to sigD log
		 sink()
		 # plot densities for FDR and p-vals
		 svg("£££.p-valuePlot.svg")
		 plot(density(sigD$p.value))
		 dev.off()
		 svg("£££.fdrPlot.svg")
		 plot(density(sigD$fdr))
		 dev.off()
		 # get number of fdr < 0.05
		 #sum(hg18_data$fdr < 0.05)
		 # plot interaction types
		 svg("£££.interactionTypes.svg")
		 plotInteractionAnnotations(sigD_subset)
		 dev.off()
			# plot cis and trans percentage
		 svg("£££.cisTransAll.svg")
		 plotCisTrans(sigD)
		 dev.off()
		 svg("£££.cisTransStrict.svg")
		 plotCisTrans(sigD_subset)
		 dev.off()
		 # plot read counts to check if subsetting by strict FDR drops low counts
		 svg("£££.countsAll.svg")
		 plotCounts(sigD, cut=50)
		 dev.off()
		 svg("sigD_countsStrict.svg")
		 plotCounts(sigD_subset, cut=50)
		 dev.off()
		 
		cat("location1,geneIDs1,geneSymbols1,location2,geneIDs2,geneSymbols2,counts,p-value,logP-val,fdr,z-score,circos-thickness\n", file="£££.annotation.txt", append=TRUE)
		my_annot_function <- function() {
			for (i in 1:length(sigD_subset)) {
				# store current chr number
				chrN <- as.character(sigD_subset[i,]@anchor_one@seqnames@values)
				# store second chr number
				chrN2 <- as.character(sigD_subset[i,]@anchor_two@seqnames@values)
				# access anchor 1 start location
				s1 <- sigD_subset[i,]@anchor_one@ranges@start
				# access anchor 1 end location
				e1 <- sigD_subset[i,]@anchor_one@ranges@width-1 + sigD_subset[i,]@anchor_one@ranges@start
				# access anchor 2 start
				s2 <- sigD_subset[i,]@anchor_two@ranges@start
				# access anchor 2 end
				e2 <- sigD_subset[i,]@anchor_two@ranges@width-1 + sigD_subset[i,]@anchor_two@ranges@start
				# access Entrez Gene IDs for anchor 1
				gID1 <- unlist(sigD_subset[i,]@anchor_one$promoter.id)
				# access Entrez Gene IDs for anchor 2
				gID2 <- unlist(sigD_subset[i,]@anchor_two$promoter.id)
				# access counts (number of interactions)
				count1 <- sigD_subset[i,]@counts
				# access p-val
				pv1 <-  sigD_subset[i,]$p.value
				#access logP value
				lpv1 <- sigD_subset[i,]$LogP
				# access fdr
				fdr1 <- sigD_subset[i,]$fdr
				# access Z-score
				zscore1 <- sigD_subset[i,]$Z.score
				# access circos thickness
				circ1 <- sigD_subset[i,]$Circos.Thickness
				# map gene entrez IDs to gene symbols
				name1 <- unlist(xx[gID1])
				name2 <- unlist(xx[gID2])
				cat(sprintf("%s:%s-%s,%s,%s,%s:%s-%s,%s,%s,%d,%.10f,%s,%f,%f,%d\n", chrN, s1, e1, paste(gID1, collapse=" "), paste(name1, collapse=" "), chrN2, s2, e2, paste(gID2, collapse=" "),  paste(name2, collapse=" "),count1, pv1, lpv1, fdr1, zscore1, circ1), file="£££.annotation.txt", append=TRUE)
			}
		}
		my_annot_function()
		mytime2 <- Sys.time()
		cat(sprintf("Annotation of £££ finished on: %s\n", mytime2),file="../annotation_timing.txt", append=TRUE)
		q("no")
		};
				print SCRIPT $Rthing;

	};		
		
		
		
		

		close(SCRIPT);

		open (FixR, "$sigIntsFile.RSCRIPT.R");
		open (NewR, ">$sigIntsFile.FixedR.R");
		chomp (my @sub_script = (<FixR>));
		foreach (@sub_script) {
			$_ =~ s/£££/$sigIntsFile/;
			print NewR "$_\n";
		}
		close(NewR); close(FixR);
	}

	## execute all scripts by sending to background,,, on grid send each as sepearate process

	unlink("sigints100k.RSCRIPT.R", "sigints40k.RSCRIPT.R", "sigints70k.RSCRIPT.R", "sigints20k.RSCRIPT.R", "sigints1M.RSCRIPT.R", "sigints250k.RSCRIPT.R");

	system(q(chmod 777 *.R));


	# divide number of threads to split annotation into 6 parts
	#my $cores = $hicup_param{threads};
	#my $each = int($cores/6);

	# run graphs and then write annotation in R
	system(qq(qsub -cwd -j y -b yes -N Anno100k -P NGS -o $lwd -q small.q /NGS/Software/HiC/R-patched/bin/R CMD BATCH sigints100k.FixedR.R));
	system(qq(qsub -cwd -j y -b yes -N Anno40k -P NGS -o $lwd -q small.q /NGS/Software/HiC/R-patched/bin/R CMD BATCH sigints40k.FixedR.R));
	system(qq(qsub -cwd -j y -b yes -N Anno70k -P NGS -o $lwd -q small.q /NGS/Software/HiC/R-patched/bin/R CMD BATCH sigints70k.FixedR.R));
	system(qq(qsub -cwd -j y -b yes -N Anno20k -P NGS -o $lwd -q small.q /NGS/Software/HiC/R-patched/bin/R CMD BATCH sigints20k.FixedR.R));
	system(qq(qsub -cwd -j y -b yes -N Anno250k -P NGS -o $lwd -q small.q /NGS/Software/HiC/R-patched/bin/R CMD BATCH sigints250k.FixedR.R));
	system(qq(qsub -cwd -j y -b yes -N Anno1M -P NGS -o $lwd -q small.q /NGS/Software/HiC/R-patched/bin/R CMD BATCH sigints1M.FixedR.R));	
	
	
}
