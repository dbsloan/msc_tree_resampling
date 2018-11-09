# msc_tree_resampling
Tree-level bootstrap and jackknife resampling for gene-tree datasets in multispecies coalescent analysis

Scripts for automating gene-tree-level jackknife and bootstrap resampling for multispecies coalescent (MSC) methods.


## Overview: 

The main Perl script (msc_tree_resampling.pl) implements a resampling scheme (either bootstrapping or jackknifing) to generate pseudoreplicates of an original dataset consisting of multiple gene trees. Optionally, the script can then call one of four MSC programs (ASTRAL, MP-EST, NJst, or STAR) to infer a species tree for each of the resulting pseudoreplicate datasets. This tool is also distributed with two additional R scripts (star.r and njst.r), which are called by the main script when running STAR or NJst and should not be run directly by the user.


## Citation: 

Simmons MP, Sloan DB, Springer MS, Gatesy J. In Press. [Gene-wise resampling outperforms site-wise resampling in phylogenetic coalescence analyses](https://www.sciencedirect.com/science/article/pii/S1055790318301404). Molecular Phylogenetics and Evolution.


## Minimum Requirements: 

This automation is implemented with a Perl script that has been designed for a Unix environment (Mac OSX or Linux). It has been tested in Mac OSX 10.12 and Linux CentOS 6, but it should work in most Unix environments.

Perl - The provided Perl script should be called by users (msc_tree_resampling.pl). Perl is pre-installed in most Mac OSX and Linux distributions.


## Requirements for Specific MSC Programs:

Calling one of the four optional MSC analysis programs will have additional dependencies.


#### ASTRAL Requirements

ASTRAL - The Perl script calls [ASTRAL](https://github.com/smirarab/ASTRAL), which must be installed, and the user must provide the full path and file name for the ASTRAL jar file. The script has been tested with ASTRAL 4.10.5 and 4.11.2 but would likely work with other versions.

Java - ASTRAL is written in Java, so Java JDK should be installed and in your PATH.

NOTE - ASTRAL does not consistently root the inferred species trees among pseudoreplicates. So it is important to ensure that all of these trees are consistently rooted before you calculate a consensus tree.  You may do so using the -rr and -names commands from [Phyutility](https://github.com/blackrim/phyutility), for example.


#### MP-EST Requirements

MP-EST - The Perl script calls [MP-EST](http://faculty.franklin.uga.edu/lliu/content/mp-est), which must be installed, and the user must provide the full path and file name for the MP-EST executable. The script has been tested with MP-EST v1.5 and v1.6 but would likely work with other versions. The user must also provide the main content for the control file that is used to run MP-EST. See format in the sample_data: gene_trees_mpest.ctl.txt. The input file should be missing the first 4 lines and the last line of a typical MP-EST control file (see MP-EST documentation), which will be added in automatically by the script when calling MP-EST. This will include a random seed generated with the Perl "rand" function for each individual pseudoreplicate dataset. The user can specify the number of tree searches to perform for each MP-EST analysis using the --mpest_num option (default is to perform a single search). Only the best (highest likelihood) tree will be reported.


#### NJst and STAR Requirements

R - The NJst and STAR methods are part of the R package [Phybase](http://faculty.franklin.uga.edu/lliu/phybase). Therefore, R and the Phybase package must be installed. This also requires the installation of additional R packages: ape, Matrix, and methods. The Rscript application that allows running R code must be either in your PATH or you must provide a full path to it when calling msc_tree_resampling.pl (see below). The scripts have been tested on Phybase 1.4 and 1.5 but does **not** currently work on v2.0. Note that these methods are distance-based and, therefore, require that every possible pair of species co-occur in at least one gene tree. In cases of resampling datasets with a lot of missing data, it is possible that a given pseudoreplicate dataset will not meet this requirement. In such cases, STAR and NJst report a warning rather than a species tree, and this will appear in the final output file.


## Output Format. 

The script will produce pseudoreplicate tree files in the specified output directory. They will be named Rep_XXX.tre, where XXX ranges from 1 to the number of specified resampling replicates. Each file contains one newick-formatted gene tree on each line (the total number of gene trees will depend on the size of the original dataset and the resampling method). If one of the four MSC programs was specified with --method, the output directory will also contain a file called MSC_trees.tre. This file will contain one newick-formatted species tree on each line. The number of species trees should equal the number of pseudoreplicate datasets, and the ordering of the species trees within MSC_trees.tre follows the numerical order of the pseudoreplicate resampling files. You may use the MSC_trees.tre file to calculate a majority-rule consensus of the bootstrap or jackknife replicates using another program, such as [Phyutility](https://github.com/blackrim/phyutility), [Mesquite](http://mesquiteproject.org/), and [PAUP*](http://paup.phylosolutions.com/).  Finally, you may map the bootstrap or jackknife percentages onto the optimal tree for the original data set using the "Add support values..." function of [TreeGraph2](http://treegraph.bioinfweb.info/).


## Running msc_tree_resampling.pl:

The script can be called from the command line, and the user must specify a number of arguments, including the gene tree file(s) to input, the resampling parameters, and the MSC method (if any) to apply. A summary of required and optional arguments is provided below, followed by example calls.


Usage: perl msc_tree_resampling.pl [arguments]
   
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


   Specify jackknife deletion probability [required if --sampling=jk]
   
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
   
   	     --Rscript      - default: Rscript

   	     
   Outgroup species:
   
   	     --outgroup     - name of outgroup taxon (must be a single species, 
   	                      not a clade). [required if --method=star]


   Path and filename for the star.r code on local machine
     [required if --method=star and not in local directory]
   
   	     --star_r_code  - default: r_code/star.r

   	     
   Path and filename for the njst.r code on local machine
     [required if --method=njst and not in local directory]
   
   	     --njst_r_code  - default: r_code/njst.r

   	     
   Directory for output files (will be created and should not 
   exist already) [required]:

   	     --output_dir   - name of directory for output files

   
   EXAMPLES 

         perl msc_tree_resampling.pl 
            --gt_file=sample_data/gene_trees.tre  
            --output_dir=sample_output 
            --reps=100 
            --sampling=jk 
            --jkprob=0.5 
            --method=astral 
            --astral_jar=/PATH/TO/ASTRAL/astral.4.10.5.jar

        perl msc_tree_resampling.pl 
            --gt_file=sample_data/gene_trees.tre  
            --output_dir=sample_output 
            --reps=100 
            --sampling=jk 
            --jkprob=0.5 
            --method=mpest 
            --mpest_file=/PATH/TO/MPEST/mpest 
            --mpest_ctl=sample_data/gene_trees_mpest.ctl.txt
            --mpest_num=5

        perl msc_tree_resampling.pl 
            --gt_file=sample_data/gene_trees.tre  
            --output_dir=sample_output 
            --reps=100 
            --sampling=jk 
            --jkprob=0.5 
            --method=star 
            --outgroup=xantusia_vigilis

        perl msc_tree_resampling.pl 
            --gt_file=sample_data/gene_trees.tre  
            --output_dir=sample_output 
            --reps=100 
            --sampling=jk 
            --jkprob=0.5 
            --method=njst

