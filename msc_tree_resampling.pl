#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use File::Which;


our $gt_file;
our $gt_dir;
our $sampling;
our $reps;
our $jkprob;
our $method;
our $astral_jar;
our $mpest_file;
our $mpest_ctl;
our $mpest_num = 1;
our $Rscript = "Rscript";
our $star_r_code = "r_code/star.r";
our $njst_r_code = "r_code/njst.r";
our $outgroup;
our $output_dir;


my $usage = 
"\nUsage: perl $0 [arguments]
   
   ARGUMENTS
   
   Gene trees must be specified with one of the following options (but
   not both) [required]: 

         --gt_file      - a single file containing one or more trees

         --gt_dir       - a directory containing one or more files each 
                          containing a single tree 

   Specify number of sampling replicates [required]:
   
         --reps         - Number of bootstrap/jackknife pseudoreplicates


   Specify resampling technique [required]:

         --sampling     - choices: bs or jk
                            
                          bs: Bootstrap
                            
                          jk: Jackknife. Must specify jkprob (see below).


   Specify deletion resampling probability [required if --sampling=jk]
   
         --jkprob       - Probability (between 0 and 1) of deleting a 
                          given tree in jackknife analysis or, to use a 
                          deletion probability of 1/e following Farris 
                          et al. (1996), specify --jkprob=Farris


   Specify method for species tree inference [required]:
   
         --method       - choices: astral, mpest, star, njst, or none
                          
                          astral: ASTRAL. Must specify astral_jar (see below).
                          
                          mpest: MP-EST. Must specify mpest_file and mpest_ctl 
                          		 (see below).  
                          
                          star: STAR. Must specify outgroup and must specify 
                                Rscript and star_r_code if not default 
                                (see below).
                          
                          njst: NJst. Must specify Rscript and njst_r_code if
                                not default (see below).
                          
                          none: No species tree inference. Only output
                                bootstrap/jackknife resamplings.


   Astral jar file name [required if --method=astral]
   
   	     --astral_jar   - Astral jar file (including path)


   MP-EST executable file name [required if --method=mpest]
   
   	     --mpest_file   - MP-EST executable file (including path)


   MP-EST control file template [required if --method=mpest]
   
   	     --mpest_ctl    - MP-EST control file template (including path)
                          See example file in sample data.


   Number of MP-EST searches to perform [Only relevant if --method=mpest]
   
   	     --mpest_num    - default: 1


   Path and filename for the Rscript application installed on local machine
     [required if --method=star or --method=njst and not already in PATH]
   
   	     --Rscript      - default: $Rscript

   	     
   Outgroup species:
   
   	     --outgroup     - name of outgroup taxon (must be a single species, 
   	                      not a clade). [required if --method=star]


   Path and filename for the star.r code on local machine
     [required if --method=star and not in local directory]
   
   	     --star_r_code  - default: $star_r_code

   	     
   Path and filename for the njst.r code on local machine
     [required if --method=njst and not in local directory]
   
   	     --njst_r_code  - default: $njst_r_code

   	     
   Directory for output files (will be created and should not 
   exist already) [required]:

   	     --output_dir   - name of directory for output files

   
   EXAMPLES 

         perl $0 
            --gt_file=sample_data/gene_trees.tre  
            --output_dir=sample_output 
            --reps=100 
            --sampling=jk 
            --jkprob=0.5 
            --method=astral 
            --astral_jar=/PATH/TO/ASTRAL/astral.4.10.5.jar

        perl $0 
            --gt_file=sample_data/gene_trees.tre  
            --output_dir=sample_output 
            --reps=100 
            --sampling=jk 
            --jkprob=0.5 
            --method=mpest 
            --mpest_file=/PATH/TO/MPEST/mpest 
            --mpest_ctl=sample_data/gene_trees_mpest.ctl.txt
            --mpest_num=5

        perl $0 
            --gt_file=sample_data/gene_trees.tre  
            --output_dir=sample_output 
            --reps=100 
            --sampling=jk 
            --jkprob=0.5 
            --method=star 
            --outgroup=xantusia_vigilis

        perl $0 
            --gt_file=sample_data/gene_trees.tre  
            --output_dir=sample_output 
            --reps=100 
            --sampling=jk 
            --jkprob=0.5 
            --method=njst

\n\n";

GetOptions(
    'gt_file=s'  => \$gt_file,
    'gt_dir=s'  => \$gt_dir,
    'reps=s'  => \$reps,
    'sampling=s'  => \$sampling,
    'jkprob=s'  => \$jkprob,
    'method=s'  => \$method,
    'astral_jar=s'  => \$astral_jar,
    'mpest_file=s'  => \$mpest_file,
    'mpest_ctl=s'  => \$mpest_ctl,
    'mpest_num=s'  => \$mpest_num,
    'Rscript=s'  => \$Rscript,
    'outgroup=s'  => \$outgroup,
    'star_r_code=s'  => \$star_r_code,
    'njst_r_code=s'  => \$njst_r_code,
    'output_dir=s'  => \$output_dir,
);


#Check that required command line arguments were provided

$gt_file or $gt_dir or die ("\n$usage\nERROR: Must specify gene trees with either --gt_file or --gt_dir\n\n");

$gt_file and $gt_dir and die ("\n$usage\nERROR: --gt_file and --gt_dir cannot both be used at the same time\n\n");

if ($gt_file){
	-d $gt_file and die ("\n$usage\nERROR: the --gt_file option was used but $gt_file is a directory. You may want to use --gt_dir instead.\n\n");
}else{
	-d $gt_dir or die ("\n$usage\nERROR: the --gt_dir option was used but $gt_dir is NOT a directory. You may want to use --gt_file instead.\n\n");
}

$output_dir or die ("\n$usage\nERROR: Must specify a directory name for output files with --output_dir.\n\n");

-e $output_dir and die ("\nERROR: $output_dir already exists. Please delete output directory before running.\n\n");

$method or die ("\n$usage\nERROR: Must specify a method for species tree inference with --method.\n\n");

unless ($method eq 'astral' || $method eq 'mpest' || $method eq 'star' || $method eq 'njst' || $method eq 'none'){
	die ("\n$usage\nERROR: $method is not an option for --method. Please specify 'astral', 'mpest', 'star', 'njst', or 'none'.\n\n");
}

if ($method eq 'astral'){
	$astral_jar or die ("\n$usage\nERROR: The --astral_jar option must be used to specify the path and filename of the ASTRAL jar file when --method=astral.\n\n");
} 

if ($method eq 'mpest'){
	$mpest_file or die ("\n$usage\nERROR: The --mpest_file option must be used to specify the path and filename of the MP-EST executable file when --method=mpest.\n\n");
	$mpest_ctl or die ("\n$usage\nERROR: The --mpest_ctl option must be used to specify the path and filename of the MP-EST control file template when --method=mpest.\n\n");
	looks_like_number($mpest_num) or die ("\n$usage\nERROR: The specified value for --mpest_num is not numeric ($mpest_num).\n\n");
} 

if ($method eq 'star' || $method eq 'njst'){
	which $Rscript or die ("\n$usage\nERROR: --method=$method was specified but $Rscript does not exist. Provide correct path to the Rscript application on your system with --Rscript.\n\n");
}

if ($method eq 'star'){
	-e $star_r_code or die ("\n$usage\nERROR: --method=star was specified but cannot find $star_r_code in local directory. Specify a full path to $star_r_code with --star_r_code.\n\n");
	$outgroup or die ("\n$usage\nERROR: --method=star was specified but no outgroup was provided. Specify the name of the outgroup with --outgroup.\n\n");
} 

if ($method eq 'njst'){
	-e $njst_r_code or die ("\n$usage\nERROR: --method=star was specified but cannot find $njst_r_code in local directory. Specify a full path to $njst_r_code with --njst_r_code.\n\n");
} 


$reps or die ("\n$usage\nERROR: Must specify a number of pseudoreplicates with --reps.\n\n");

looks_like_number($reps) or die ("\n$usage\nERROR: Must specify a number of pseudoreplicates with --reps. The specified value ($reps\) is not numeric.\n\n");


$sampling or die ("\n$usage\nERROR: Must specify a sampling approach with --sampling. Options are bs (for boostrap) or jk (for jackknife).\n\n");

unless ($sampling eq 'bs' || $sampling eq 'jk'){
	die ("\n$usage\nERROR: The specified value ($sampling\) is not a valid option for --sampling. Options are bs (for boostrap) or jk (for jackknife).\n\n");
}

if ($sampling eq 'jk'){
	$jkprob or die ("\n$usage\nERROR: Jackknife sampling was chosen (--sampling=jk) but no resampling probability was provided with --jkprob.\n\n");
	$jkprob eq 'Farris' and $jkprob = 1/exp(1);
	unless (looks_like_number($jkprob) && $jkprob > 0 && $jkprob < 1){
		die ("\n$usage\nERROR: Specified value for jackknife sampling probability ($jkprob\) is not allowable. Use --jkprob to specify a value greater than 0 and less than 1.\n\n");
	}
}

mkdir ($output_dir);

#Define random number that will be used in saving temp file names to avoid clashing if multiple runs are done in parallel.
my $randnum = int(rand(10000000));

#concatenate files from $gt_dir if user provided a directory of gene trees
if ($gt_dir){
	substr ($gt_dir, -1) eq '/' or $gt_dir .= '/';
	my @files = get_file_names ($gt_dir);
	my $FH_GT_OUT = open_output (".$randnum\_genetrees.txt");
	$gt_file = ".$randnum\_genetrees.txt";
	foreach (@files){
		my $FH_GT = open_file($gt_dir.$_);
		my $first_line = <$FH_GT>;
		chomp $first_line;
		print $FH_GT_OUT "$first_line\n";
	}
}


#extract trees
my @original_trees;

my $FH = open_file($gt_file);

while (<$FH>){
	chomp $_;
	$_ =~ /^\s*$/ and next;
	push (@original_trees, $_)
}
my $tree_num = scalar(@original_trees);
close $FH;

#if $gt_file was created as temp file from $gt_dir, then delete it (but do not delete user provided file).
$gt_dir and unlink($gt_file);


print STDERR "\n\n". (localtime) ."\nRunning $0\n\n";

if ($sampling eq "bs" && $jkprob){
	print STDERR "WARNING: --jkprob was specified but is being ignored because bootstrap sampling was requested (--sampling=bs). If you want to run a jackknife analysis, use --sampling=jk.\n\n";
}
my $gt_source;
if ($gt_file){$gt_source = $gt_file;}else{$gt_source = $gt_dir;}
print STDERR "Analyzing gene trees in $gt_source\. Found $tree_num trees.\n\nGenerating $reps pseudoreplicate datasets with $sampling sampling method. Output files being written to $output_dir\.\n\n";

#implement bootstrap sampling by selection trees based on a random integer between (0 to number of trees - 1)
if ($sampling eq "bs"){
	for (my $i = 1; $i <= $reps; ++$i){
		my $FH_BSO = open_output("$output_dir\/Rep_$i\.tre");
		for (my $j = 1; $j <= $tree_num; ++$j){
			my $tree = $original_trees[int(rand($tree_num))];
			print $FH_BSO "$tree\n";
		}
		close $FH_BSO;
	}
}

#implement jackknife by comparing random number between 0 and 1 to specified sampling probability
if ($sampling eq "jk"){
	for (my $i = 1; $i <= $reps; ++$i){
		my $FH_JKO = open_output("$output_dir\/Rep_$i\.tre");
		foreach (@original_trees){
			rand() > $jkprob && print $FH_JKO "$_\n";	
		}
		close $FH_JKO;
	}
}


#run multispecies coalescent analysis on each pseudoreplicate dataset
unless ($method eq 'none'){
	
	print STDERR "Performing multispecies coalescent analysis on each pseudoreplicate dataset using $method method. Output being written to $output_dir\/MSC_trees.tre\n\n";
	
	my $mpest_ctl_text;
	if ($method eq 'mpest'){
		$mpest_ctl_text = file_to_string($mpest_ctl);
	}	
	
	#store all MSC trees in a single file in the order of the pseudoreplicate datasets
	my $FH_TO = open_output("$output_dir\/MSC_trees.tre");
	for (my $i = 1; $i <= $reps; ++$i){
		
		#if astral was selected, simply run under default settings. Dump resulting trees (stdout) into the above output file and stderr into a log file
		if ($method eq 'astral'){
			system ("java -jar $astral_jar -i $output_dir\/Rep_$i\.tre >> $output_dir\/MSC_trees.tre 2>> $output_dir\/ASTRAL_log.txt");
		}
		
		#if mp-est was selected, construct a control file. The user provides the bottom part of the file. 
		#The scripts adds top 3 lines (input tree file name, 0 to turn off counting of triple distances [could be changed to 1 if desired], and random seed).
		#Run mpest with system call. Then open output (just has an extra appended .tre or _besttree.tre extension, depending on the version of MP-EST), convert to a simple newick string, and add to main MSC_tree output file
		#Finally, delete temp control file and the individual replicate output files.
		elsif ($method eq 'mpest'){
			my $FHCTL = open_output(".$randnum\_TEMP_CONTROL_FILE");
			print $FHCTL "$output_dir\/Rep_$i\.tre\n0\n", int(rand(1000000000)), "\n$mpest_num\n$mpest_ctl_text\n0";
			system ("$mpest_file .$randnum\_TEMP_CONTROL_FILE >> $output_dir\/MP-EST_log.txt");
			my $mpest_output_file;
			if (-e "$output_dir\/Rep_$i\.tre\.tre"){
				$mpest_output_file = "$output_dir\/Rep_$i\.tre\.tre";
			}elsif(-e "$output_dir\/Rep_$i\.tre\_besttree.tre"){
				$mpest_output_file = "$output_dir\/Rep_$i\.tre\_besttree.tre";
			}else{
				die ("\nERROR: Could not identify MP-EST file output file. Looked for Rep_$i\.tre\_besttree.tre or Rep_$i\.tre\.tre in $output_dir\n\n");
			}
			
			my @mpest_output = file_to_array($mpest_output_file);
			my $block_flag = 0;
			my %taxa_hash;
			my $best_tree;
			my $best_score;
			
			foreach my $mpest_line (@mpest_output){
				if ($mpest_line =~ /^\s+translate/){
					$block_flag = 1;
					next;
				}
				
				if ($block_flag){
					if ($mpest_line =~ /^\s+tree\smpest\s\[([\-\d\.]+)\]\s\=\s(\(.+\)\;)/){
						my $lnl_score = $1;
						my $newick_string = $2;
						
						unless ($best_tree){
							$best_score = $lnl_score;
							$best_tree = $newick_string;
							next;
						}
						
						if ($lnl_score > $best_score){
							$best_score = $lnl_score;
							$best_tree = $newick_string;							
						}
						
					}elsif($mpest_line =~ /^\s+(\d+\s[\w\s]+)[\,\;]/){
						my $transl_line = $1;
						my @sl = split (/\s+/, $transl_line);
						my $species_name;
						for (my $array_pos = 1; $array_pos < scalar(@sl); ++$array_pos){
							$species_name .= $sl[$array_pos];
						}
						$taxa_hash{$sl[0]} = $species_name;
					}
				}
			}
			$best_tree or die ("\nERROR: Could not find a tree line in $output_dir\/$mpest_output_file\n\n");

			foreach my $key (sort keys %taxa_hash){
				my $species = $taxa_hash{$key};
				$best_tree =~ s/\($key\:/\($species\:/g;
				$best_tree =~ s/\,$key\:/\,$species\:/g;
			}


			print $FH_TO "$best_tree\n";
			unlink ".$randnum\_TEMP_CONTROL_FILE";
			unlink $mpest_output_file;
			-e "$output_dir\/Rep_$i\.tre\_output.tre" and unlink "$output_dir\/Rep_$i\.tre\_output.tre";
		}
		
		##If STAR was selected, generate the species tree by calling the corresponding R code distributed with this script.
		elsif($method eq 'star'){
			system("$Rscript $star_r_code $output_dir\/Rep_$i\.tre $outgroup $output_dir\/MSC_trees.tre");
		}

		##If NJst was selected, generate the species tree by calling the corresponding R code distributed with this script.
		elsif($method eq 'njst'){
			system("$Rscript $njst_r_code $output_dir\/Rep_$i\.tre $output_dir\/MSC_trees.tre");
		}
		my $num = $i + 1;
		
		$num % 10 or print STDERR (localtime) . "\tCompleted $num replicates.\n";
	
	}
	close $FH_TO;
}

print STDERR "\n" . (localtime) . "\tRun Completed\n\n";


sub open_file {
	use strict;
	use warnings;

    my($filename) = @_;
    my $fh;

    unless(open($fh, $filename)) {
        print "Cannot open file $filename\n";
        exit;
    }
    return $fh;
}


sub open_output {
	use strict;
	use warnings;

    my($filename) = @_;
    my $fh_output;

    unless(open($fh_output, ">$filename")) {
        print "Cannot open file $filename\n";
        exit;
    }
    return $fh_output;
}


sub file_to_string {
	use strict;
	use warnings;

    my($filename) = @_;

    # Initialize variables
    my $filedata;

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }
    
    while (<GET_FILE_DATA>){
    	$filedata .= $_;
    }
    
    close GET_FILE_DATA;

    return $filedata;
}


sub file_to_array {
	use strict;
	use warnings;

    my($filename) = @_;

    # Initialize variables
    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}