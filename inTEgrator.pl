#!/usr/bin/perl

use strict;
use warnings;
use Sort::Fields;
use Getopt::Long;
use Data::Dumper qw(Dumper);
use List::Util qw[min max shuffle];
use Array::Utils qw(:all);
use Scalar::Util qw(looks_like_number);
use Switch ;

my $header = qq{
        __________              __          
  (_)__/_  __/ __/__ ________ _/ /____  ____
 / / _ \\/ / / _// _ `/ __/ _ `/ __/ _ \\/ __/
/_/_//_/_/ /___/\\_, /_/  \\_,_/\\__/\\___/_/   
               /___/                      v0.1            

inTEgrator by A.Coudray (LVG-EPFL 2019)
} ;

my $usage_general = qq{

inTEgrator can either be launched using 'inTEgrator'
mode and therefore all at once

--mode inTEgrator   Launch the 3-step process to get
                    from a single bed file to the 
                    age estimate per TE locus. 

    ./inTEgrator --mode inTEgrator --help


Otherwise you can use separated steps used by inTE-
grator

--mode liftover_te  Launch first part of the scrip-
                    where TEs are liftOver towards
                    49 foreign species from the inp-
                    bed file (e.g. TE subfamily).

    ./inTEgrator --mode liftover_te --help

--mode pairwise_age Launch multi-sequence alignments
                    on input bed file and compute all
                    pairwise distance / ages using 
                    Kimura 2p model. Gives as an out-
                    put a matrix with pairwise age 
                    caculated from substitution dens-
                    ity.
    
    ./inTEgrator --mode pairwise_age --help


--mode merge_age    Merge liftover_age and pairwise
                    age made from substitutions as
                    does inTEgrator mode.

    ./inTEgrator --mode merge_age --help
};

my $usage_liftover_te = qq{
--mode liftover_te  Launch first part of the scrip-
                    where TEs are liftOver towards
                    49 foreign species from the inp-
                    bed file (e.g. TE subfamily).

Options :

--bed               bed file containing regions to 
                    translate in foreign genome co-
                    ordinates.

--liftover_dir      Directory containing liftOver
                    chains (by default :
                    path/to/inTEgrator/chains

--convert_to_hg38   [needs to be either 'T' or 'F']
                    if 'T' (=true), it will use bed
                    file as if it is in hg19 coord-
                    inates, then convert it to hg38
                    before translating to foreign
                    species. 

--te_anno_dir       Name of TE_anno_species directo-
                    needed to check whether TEs are
                    annotated as the same family in
                    foreign species (improve results
                    compare with liftover alone. 
                    WARNING : weights 5.1Gb !!! 
                   
--fam               To work we suppose that all elem-
                    ents comes from same TE family.
                    Should be precise by this argu-
                    ent.
};

my $usage_pairwise_age = qq{
--mode pairwise_age Launch multi-sequence alignments
                    on input bed file and compute all
                    pairwise distance / ages using 
                    Kimura 2p model. Gives as an out-
                    put a matrix with pairwise age 
                    caculated from substitution dens-
                    ity.

Options :

--ref_genome        reference genome use to convert
                    sequences in bed input in fasta
                    format for multi-seq alignment.
                    e.g. hg19.fa  


};

my $usage_merge_age = qq{
--dir               path to directory containing 
                    pairwise substitution matrices 

--suffix            suffix of matrices in --dir

--out_dir           path to output directory

--max_iter          maximum number of iteration to 
                    resolve the tree.
                   
--prefix            prefix of matrices in --dir

--max_n_align       maximum number of alignments 
                    used to resolve a node in the 
                    tree

--subfam_subset     give a file containing subfam 
                    names to use in the analysis 
                    (one per line)

--summary_dir       path to directory containing 
                    summary file of liftOver ana-
                    lysis. Should contains file
                    with .summary.txt suffix. e.g
                    summary_dir/SVA_A.summary.txt
};

my $help; 
my $pwd = `pwd` ; chomp $pwd ;
my $mode = 'empty' ;
my $bed = 'empty' ;
my $liftover_dir="empty" ;
my $convert_to_hg38 = 'T' ;
my $te_anno_dir = 'empty' ;
my $fam = 'empty' ;

# merge_age options ...
#my $dir = '' ;
#my $suffix = '' ;
#my $out_dir = '../results/default_resolve_age_2' ;
#my $max_iter = 1000 ;
#my $prefix = '' ;
#my $max_n_align = 50 ;
#my $subfam_subset = 'empty' ;
#my $summary_dir = '../data/all_summaries' ; 
#my $mode = 'empty' ;

GetOptions(
    "help" => \$help,
    "mode=s" => \$mode,
    "bed=s" => \$bed,
    "liftover_dir=s" => \$liftover_dir,
    "convert_to_hg38=s" => \$convert_to_hg38,
    "te_anno_dir=s" => \$te_anno_dir,
    "fam=s" => \$fam,
#    "dir=s" => \$dir,
#    "suffix=s" => \$suffix,
#    "out_dir=s" => \$out_dir,
#    "max_iter=s" => \$max_iter,
#    "prefix=s" => \$prefix,
#    "max_n_align=i" => \$max_n_align,
#    "subfam_subset=s" => \$subfam_subset,
#    "summary_dir=s" => \$summary_dir,
);

# Print Help and exit
if ($help) {
    if ( $mode eq 'empty') {
        print "$header\n$usage_general";
    }
    elsif ( $mode eq 'liftover_te' ){
        print "$header\n$usage_liftover_te" ;
    }
    elsif ( $mode eq 'merge_age' ){
        print "$header\n$usage_merge_age" ;
    }
    exit(0);
}

if ( $mode eq 'empty' ){
    print "$header\n$usage_general";
}
elsif ( $mode eq "liftover_te" || $mode eq "inTEgrator" ){

    if ( $liftover_dir eq 'empty' || $te_anno_dir eq 'empty' || $bed eq 'empty' || $fam eq 'emtpy' ){
        print "ERROR : liftover_dir , te_anno_dir, fam and bed should be given as arguments...\n" ;
        exit(0);
    }

    ### LIFTOVER MODE ### 
    # Get species list with/without hg38 
    my  @all_species = ("panTro5", "gorGor4", "ponAbe2", "nomLeu3", "rheMac8", "macFas5", "calJac3", "tarSyr2", "otoGar3", "micMur2", "tupBel1", "mm10", "speTri2", "cavPor3", "pteVam1", "vicPac2","turTru2", "susScr3", "felCat8", "dasNov3", "loxAfr3", "triMan1", "macEug2", "monDom5", "ornAna2", "galGal5", "anoCar2", "latCha1",  "danRer10") ;
    if ( $convert_to_hg38 eq 'True' ){ $convert_to_hg38 = 'T' }
    if ( $convert_to_hg38 eq 'T' ){
        #@all_species = ("hg38", "panTro5", "gorGor4", "ponAbe2", "nomLeu3", "rheMac8", "macFas5") ;
        @all_species = ("hg38", "panTro5", "gorGor4", "ponAbe2", "nomLeu3", "rheMac8", "macFas5", "calJac3", "tarSyr2", "otoGar3", "micMur2", "tupBel1", "mm10", "speTri2", "cavPor3", "pteVam1", "vicPac2","turTru2", "susScr3", "felCat8", "dasNov3", "loxAfr3", "triMan1", "macEug2", "monDom5", "ornAna2", "galGal5", "anoCar2", "latCha1",  "danRer10") ;
    }

    ## Create folders
    my $dir_out = "$pwd/liftover_age" ;
    `mkdir -p $dir_out` ;
    #`mkdir -p results/TE_liftOver` ;

    ### Iterate over all elements 
    my $name = $bed ; # get sample name - basename without suffix
    $name =~ s/(\.bed|.*\/|\/)//g ;
    print "$bed,$name\n";
     
    `mkdir -p $dir_out/$name` ;

    #my $name_hg19 = `cut -f1-4,6 $name_file > $dir_out/$name/hg19_${name}.bed` ; 
    my $file = "$dir_out/$name/hg19_${name}.bed" ;

    ### write original file in new format suitable for liftOver ### 
    print "rewrite original entry...\n" ;
    open my $original_tab, '<', $bed or die "can't open $bed: $!" ;
    open my $out_tab, '>', $file or die "can't open $file: $!" ;
    while(<$original_tab>){
        chomp ;
        my @line = split /\t/ ;
        if ( $line[1] eq 'start' ){ next }
        my $key = "$line[0]:$line[1]-$line[2]" ;
        print $out_tab "$line[0]\t$line[1]\t$line[2]\t$key\n" ;
    }
    close $original_tab ; close $out_tab ;

    my $file_hg19 = $file ;

    my %original; my %coord ; 
    my %name ; my %fam ;

    foreach my $species ( @all_species ){
        print "Doing liftOver to $species\n" ;

        my $out_dir = "$dir_out/$name/$species" ;
        `mkdir -p $out_dir` ;

        my $chain_file ;
        if ( $convert_to_hg38 eq 'T' && $species eq 'hg38' ){
            $chain_file="${liftover_dir}/hg19ToHg38.over.chain.gz" ;
        }
        else
        {
            #Species name pattern to grep
            my $pattern= substr $species, 1 ;
            $pattern.='.over.chain.gz' ;
            # get chain file for first liftover
            $chain_file=`ls ${liftover_dir}/hg38*|grep $pattern` ;
            chomp $chain_file ;
        }

        ### LiftOver process ### 
        my $minMatch = 0.5 ;
        my $name_out = "$out_dir/${name}_$minMatch.bed" ;
        my $name_unmap = "$out_dir/${name}_${minMatch}_unmap.bed" ;
        `liftOver -minMatch=$minMatch $file $chain_file $name_out $name_unmap 2> logs/liftOver.log` ;
        
        my $found_lines = `wc -l $name_out|sed -e 's/^[ \t]*//'|tr " " "\t"|cut -f1`  ;
        if ( ! looks_like_number($found_lines) ){ $found_lines = 0 }

        if ( $found_lines == 0 ){
            `cp $name_out $out_dir/all_${name}.bed` ;
            `rm $dir_out/$name/$species/*_unmap.bed` ;
            print "skipping $species for $name as nothing was found\n" ;
        }
        else
        {
            for( my $i = 0; $i < 5; $i+=2 ){

                my $minMatch = 0.9 - ( $i / 10 ) ;

                if ( $i == 0 ){
                    $minMatch = 0.9 ;
                    my $name_out = "$out_dir/${name}_$minMatch.bed" ;
                    my $name_unmap = "$out_dir/${name}_${minMatch}_unmap.bed" ;
                    `liftOver -minMatch=$minMatch $file $chain_file $name_out $name_unmap 2> logs/liftOver.log` ;
                }
                elsif ( $i < 10 )
                {
                    my $prev_minMatch ;
                    if ( $minMatch == 0.9 )
                    { 
                        $prev_minMatch = 0.99  ;
                    }
                    else
                    {
                        $prev_minMatch = 1 - ( ( $i - 1 ) / 10 ) ;
                    }
                    my $name_out = "$out_dir/${name}_$minMatch.bed" ;
                    my $name_unmap = "$out_dir/${name}_${minMatch}_unmap.bed" ;
                    my $name_prev_unmap = "$out_dir/${name}_${prev_minMatch}_unmap.bed" ; 
                    `liftOver -minMatch=$minMatch $name_prev_unmap $chain_file $name_out $name_unmap 2> logs/liftOver.log` ;
                }
            }
            `rm $dir_out/$name/$species/*_unmap.bed` ;

            print "concatenate liftOver files...\n" ;
            my @all_lo_files = `ls $out_dir/${name}*.bed` ;
            my $concat_file = "$out_dir/all_${name}.bed" ;
            if ( -e $concat_file ){ `rm $concat_file` } 
            open my $concat_out, '>>', $concat_file or die "can't open $concat_file: $!" ;
            foreach ( @all_lo_files ){
                chomp ;
                open my $this_tab, '<', $_ or die "can't open $_: $!" ;
                while(<$this_tab>){
                    chomp ;
                    print $concat_out "$_\n" ; 
                }
                close $this_tab ;
            }
            close $concat_out ;
        }

        ## Getting TEs intersecting with same fam in foreign species
        my $TEfam_annotation_file = "${te_anno_dir}/${species}/fam/${fam}.bed.gz" ;
        if ( -e $TEfam_annotation_file ){
            print "Intersect with $fam TEs in $species ...\n" ;
            my $file_temp_fam = "$out_dir/temp_$species.$fam.intersect.bed" ;
            `bedtools intersect -a $out_dir/all_${name}.bed -b $TEfam_annotation_file -f 0.5 -F 0.5 -e > $file_temp_fam` ;
            open my $tab_temp_fam, '<', $file_temp_fam or die "can't open $file_temp_fam: $!" ;
            while(<$tab_temp_fam>){
                chomp ;
                my @line = split /\t/ ;
                my $key = $line[3] ;
                if ( ! exists $fam{$key}->{$species} ){
                    $fam{$key}->{$species} = 1 ;
                }
            }
            close $tab_temp_fam ;
        }

        ### Make liftOver summary by iterating over files previously created 
        # 1st, put original file into hash for comparison ...
        open my $tab, '<', $file_hg19 or die "can't open $file_hg19: $!";
        while (<$tab>){
            chomp;
            my @line = split /\t/;
            my $key = "$line[0]:$line[1]-$line[2]";     # for bed file, 3 first fields
            $original{$key}->{$species} = 0 ;
            $coord{$key}->{$species} = "NA" ;
        }
        close $tab;

        my @files = `ls $out_dir/${name}*.bed` ;
        foreach my $lo_file ( @files ){
            chomp $lo_file ;
            my @two_fields = split(/_([^_]+)$/, $lo_file) ;
            my $alignment_range = "$two_fields[1]" ;
            $alignment_range =~ s/.bed//g ;

            open my $tab2, '<', $lo_file or die "cant open $lo_file: $!";
            while (<$tab2>){
                chomp;
                my @line = split /\t/;
                my $key = $line[$id_field] ; 
                my $foreign_key = "$line[0]:$line[1]-$line[2]" ;

                if ( exists $original{$key} ){ 
                    $original{$key}->{$species} = $alignment_range ;
                    $coord{$key}->{$species} = $foreign_key ;
                }
            }
            close $tab2;
        }

        # if we want to convert to hg38, take hg38 with 90% id as a starting file for next steps
        if ( $convert_to_hg38 eq 'T' && $species eq 'hg38' ){
            `cat $out_dir/$name* > $out_dir/hg38_ref_seq_1.0.bed` ;
            $file = "$out_dir/hg38_ref_seq_1.0.bed" ;
        }

        ### Printing results ###
        # 2 outputs : one normalize by consensus length, the other not
        `mkdir -p $dir_out/tables` ;
        my $output = "$dir_out/tables/${name}.txt" ;
        open my $out_reg, '>', $output or die "can't open $output: $!" ;
        print $out_reg "integrant\t", join("\t",@all_species), "\n" ;

        my $output_coord = "$dir_out/tables/${name}.coord.txt" ;
        open my $out_coord, '>', $output_coord or die "can't open $output_coord: $!" ;
        print $out_coord "integrant\t", join("\t",@all_species), "\n" ;

        my $output_lo_age = "$dir_out/tables/${name}.lo_age.txt" ;
        open my $out_lo_age, '>', $output_lo_age or die "can't open $output_lo_age: $!" ;
        print $out_lo_age "integrant\tlo_age\n" ; 

        foreach (sort keys %original){
            my $key = $_ ;
            print $out_reg "$key" ;
            print $out_coord "$key" ;

            my $elem = $key ;
            $elem =~ s/.*://g ;

            ## Measuring integrant length
            my @sep1 = split(":",$key) ;
            my @sep2 = split("-",$sep1[1]) ;
            my $elem_length = $sep2[1] - $sep2[0] ;
            my $max_age = 0 ; my $i = 1 ;

            foreach my $spec ( @all_species ){
                if ( exists $original{$_}->{$spec} ){
                    if ( exists $fam{$_}->{$spec} ){
                        my $final_reg = $original{$_}->{$spec} + 1 ; #  + $to_add ;
                        print $out_coord "\t$coord{$_}->{$spec}" ;
                        print $out_reg "\t$final_reg" ;
                        my $this_age = &get_lo_age($i) ;
                        if ( $this_age > $max_age ){
                            $max_age = $this_age ;
                        }
                    } else {
                        my $final_reg = $original{$_}->{$spec} ; 
                        print $out_coord "\t0" ;
                        print $out_reg "\t$final_reg" ;
                    }
                }
                else
                {
                    print $out_reg "\t0" ;
                    print $out_coord "\t0" ;
                }
                $i++ ;
            }
            print $out_reg "\n" ;
            print $out_coord "\n" ;
            print $out_lo_age "$key\t$max_age\n" ;
        }
        close $out_reg ; close $out_coord ;
    }
    `rm -rf $dir_out/$name` ;
}
if( $mode eq "inTEgrator" || $mode eq "merge_age" ){
    #    # Record subfam_subset
    #    my %sfam_todo ;
    #    open my $sfam_subset, '<', $subfam_subset or die "can't open $subfam_subset: $!" ;
    #    while(<$sfam_subset>){
    #        chomp;
    #        $sfam_todo{$_} = 1 ;
    #    }
    #    close $sfam_subset ;
    #
    #
    #    my @files = <$dir/$prefix*$suffix>;
    #    foreach my $file (@files) {
    #
    #        my $subfam_name = $file ; # get sample name - basename without suffix
    #        $subfam_name =~ s/($dir|$suffix|\/)//g ;
    #
    #        if ( exists $sfam_todo{$subfam_name} ){
    #            print "Launching analysis for $subfam_name\n";
    #        }
    #        else
    #        {
    #            print "skipping $subfam_name, not in subset...\n" ;
    #            next ;
    #        }
    #
    #        ### Put summary file with cluster and LO age information 
    #        print "getting summary files...\n" ;
    #        my $summary_file = "$summary_dir/$subfam_name.summary.txt" ;
    #        if ( ! -e $summary_file ){
    #            print "no summary file for $subfam_name, skipping...\n" ;
    #            next ;
    #        }
    #        my %clusters ; my %lo_age ; my %summary ;
    #        open my $summ, '<', $summary_file or die "can't open $summary_file: $!" ;
    #        while(<$summ>){
    #            chomp ;
    #            my @line = split /\t/ ;
    #            my $key = $line[0] ;
    #            $lo_age{$key} = $line[5] ;
    #            $summary{$key} = $_ ;
    #        }
    #        close $summ ;
    #
    #        ### STEP 1 : Get youngest age for each integrants and get pairs ###
    #        my %all_pairs ; my %groups ;
    #        my %youngest ; my @all_names ; my $first_line = 'Y' ; my $line_count = -1 ; 
    #        my %youngest_age ; my %idx_to_name ; my %name_to_idx ; my %status ; 
    #        my $mean_age = 0 ; my $total_count = 0 ;
    #        # First pass through matrix to get pairs
    #        open my $tab, '<', $file or die "cant open $file: $!";
    #        while (<$tab>){
    #            chomp;
    #            my @line = split / / ;
    #            my $minimum = 1000 ; my $field_count = -1 ;
    #            my $name = '' ; my $first = 'Y' ;
    #            foreach ( @line ){
    #                if ( $_ eq 'NA' ){ next } # ignore NAs
    #                if ( $first eq 'F' && $first_line eq 'F' ){ 
    #                    $mean_age += $_ ; $total_count++ ;
    #                    $all_pairs{$name}->{$all_names[$field_count]} = $_ ;
    #                }
    #                if ( $first_line eq 'Y' && $first eq 'F' ){ # record header
    #                    push @all_names, $_ ;
    #                    my @seq_name = ($_) ;
    #                    $groups{$_} = [ @seq_name ] ; # initilize groups, all seq has itself
    #                    $status{$_} = 'NA' ;
    #                }
    #                elsif ( $first eq 'Y' ){ ## avoid first field
    #                    $name = $_ ; $first = 'F' ;
    #                }
    #                elsif ( $_ < $minimum && $name ne $all_names[$field_count] ){ 
    #                    $minimum = $_ ;
    #                    $youngest{$name} = $all_names[$field_count] ; 
    #                    $youngest_age{$name} = $_ ;
    #                }
    #                $name_to_idx{$name} = $field_count ;
    #                $idx_to_name{$field_count} = $name ;
    #                $field_count++ ;
    #
    #            }
    #            $first_line = 'F' ; $line_count++ ;
    #        }
    #        close $tab;
    #
    #        $mean_age = $mean_age / $total_count ;
    #        my $NA_replacement = $mean_age ;
    #        my $size_subfam = $#all_names + 1 ;
    #
    #        ### STEP 2 : GET SEED WITH BEST MUTUAL MATCH ### 
    #        my %age ; my %progeny ; my %ancestry ; 
    #        my $resolved_child = 0 ; my $resolved_parent = 0 ;
    #        for (my $i = 0; $i <= $#all_names; $i ++) {
    #            my $qry_name = $all_names[$i] ; 
    #            if ( $qry_name eq 'names' ){ next }
    #            my $junior = $youngest{$qry_name} ;
    #
    #            if ( $youngest{$junior} eq $qry_name ){ 
    #                my $whos_the_child = 'NA' ;
    #                my $len_r = get_length($qry_name) ;
    #                my $len_junior = get_length($junior) ;
    #
    #                # Make test to determine who's the child/parent sequence
    #                if ( $len_r < 0.5*$len_junior ){
    #                    $whos_the_child = 'query' ;
    #                }
    #                elsif ( $len_junior < 0.5*$len_r ){
    #                    $whos_the_child = 'subject' ;
    #                }
    #                elsif ( $lo_age{$qry_name} ne 'NA' && $lo_age{$junior} ne 'NA' && $lo_age{$qry_name} < $lo_age{$junior} )
    #                {
    #                    $whos_the_child = 'query' ;
    #                }
    #                elsif ( $lo_age{$qry_name} ne 'NA' && $lo_age{$junior} ne 'NA' && $lo_age{$qry_name} > $lo_age{$junior} )
    #                {
    #                    $whos_the_child = 'subject' ;
    #                }
    #                elsif ( $lo_age{$qry_name} eq 'NA' || $lo_age{$junior} eq 'NA' || $lo_age{$qry_name} == $lo_age{$junior} )
    #                {
    #                    if ( $len_r > $len_junior ){
    #                        $whos_the_child = 'subject' ;
    #                    }
    #                    elsif( $len_r < $len_junior )
    #                    {
    #                        $whos_the_child = 'query' ;
    #                    }
    #                    elsif( $len_r == $len_junior ){
    #                        if ( index($qry_name,"_") != -1 && index($junior,"_") != -1 ){
    #                            $whos_the_child = 'subject' ;
    #                        }
    #                        elsif ( index($qry_name,"_") != -1 ){
    #                            $whos_the_child = 'query' ;
    #                        }
    #                        elsif ( index($junior,"_") != -1 ){
    #                            $whos_the_child = 'subject' ;
    #                        }
    #                        else
    #                        {
    #                            $whos_the_child = 'query' ;
    #                        }
    #                    }
    #                }
    #                if ( $whos_the_child eq 'query' ){
    #                    $age{$qry_name} = $youngest_age{$qry_name} ;
    #                    $status{$qry_name} = 'child' ; $status{$junior} = 'parent' ;
    #                    $progeny{$junior}->{$qry_name} = $youngest_age{$qry_name} ;
    #                    $ancestry{$qry_name} = $junior ;
    #                    my @union = ($qry_name, $junior) ;
    #                    $groups{$qry_name} = [ @union ] ; 
    #                    $groups{$junior} = [ @union ] ;
    #                }
    #                elsif ( $whos_the_child eq 'subject' ){
    #                    $age{$junior} = $youngest_age{$junior} ;
    #                    $status{$junior} = 'child' ; $status{$qry_name} = 'parent' ; 
    #                    $progeny{$qry_name}->{$junior} = $youngest_age{$junior} ;
    #                    $ancestry{$junior} = $qry_name ;
    #                    my @union = ($qry_name, $junior) ;
    #                    $groups{$qry_name} = [ @union ] ; 
    #                    $groups{$junior} = [ @union ] ;
    #                }
    #                else
    #                {
    #                    print "PROBLEM with line $i : who's the child is $whos_the_child\n" ;
    #                }
    #                $resolved_child++ ; $resolved_parent++ ;
    #            }
    #        }
    #
    #        print "Pass1, resolved child: $resolved_child, resolved parents: $resolved_parent\n"; 
    #
    #
    #
    #        ### STEP3 : ITERATE N TIMES OVER THE MATRIX TO RESOLVE CHILD/PARENTS RELATION ### 
    #
    #        my $prev_n_child = $resolved_child ; my $prev_n_parent = $resolved_parent ;
    #        my $round_above_limit = 0 ; my $magic_pass = 'F' ;
    #        # Iterate N times until convergence  
    #        for( my $x = 1 ; $x <= $max_iter; $x++ ){
    #
    #            # Re-establish youngest 
    #            $first_line = 'T' ;
    #            open my $tab2, '<', $file or die "cant open $file: $!";
    #            while (<$tab2>){
    #                chomp;
    #                my @line = split / / ;
    #                if ($first_line eq 'T' ){ $first_line = 'F' ; next }
    #                if ( $status{$line[0]} eq 'child' ){ next }
    #                my $minimum = 10000 ; my $field_count = -1 ;
    #                my $name = $line[0] ; ## name of the query TE 
    #                for ( my $i = 1; $i <= $#line; $i++ ) { 
    #                    if ( $line[$i] eq 'NA' ){ next }#$line[$i] = $NA_replacement } # ignore NAs
    #                    my $name_comp = $all_names[($i-1)] ; #name of the subject TE for comparison
    #
    #                    if ( $line[$i] < $minimum && $name ne $name_comp && $status{$name_comp} ne 'child' ){
    #                        $minimum = $line[$i] ;
    #                        $youngest{$name} = $name_comp ; 
    #                        $youngest_age{$name} = $line[$i] ;
    #                    }
    #                }
    #                $line_count++ ;
    #            }
    #            close $tab2;
    #
    #            for (my $i = 0; $i <= $#all_names; $i ++) {
    #                my $qry_name = $all_names[$i] ; 
    #                if ( $qry_name eq 'names' ){ next }
    #                if ( $status{$qry_name} eq 'child' ){ next }
    #                my $junior = $youngest{$qry_name} ;
    #                my $status_junior = $status{$junior} ;
    #                my $trial_to_resolve = 0 ;
    #                while( $status_junior eq 'child' ){ 
    #                    $trial_to_resolve++ ;
    #                    $junior = $ancestry{$junior} ;
    #                    $status_junior = $status{$junior} ;
    #                    if ( $status_junior eq 'child' ){
    #                        #print "junior had child status, his direct ancestry has $status_junior status!!\n" ;
    #                        #print "trial number $trial_to_resolve\n" ;
    #                    } 
    #                }
    #                ##print "$qry_name,$junior\n" ;
    #                if ( $youngest{$junior} eq $qry_name || $magic_pass eq 'T' ){ 
    #
    #                    my $whos_the_child = 'NA' ;
    #                    my $len_r = get_length($qry_name) ;
    #                    my $len_junior = get_length($junior) ;
    #
    #                    ### Make test to determine who's and child and who's the parent sequence
    #                    if ( $len_r < 0.5*$len_junior ){
    #                        $whos_the_child = 'query' ;
    #                    }
    #                    elsif ( $len_junior < 0.5*$len_r ){
    #                        $whos_the_child = 'subject' ;
    #                    }
    #                    elsif ( $lo_age{$qry_name} ne 'NA' && $lo_age{$junior} ne 'NA' && $lo_age{$qry_name} < $lo_age{$junior} )
    #                    {
    #                        $whos_the_child = 'query' ;
    #                    }
    #                    elsif ( $lo_age{$qry_name} ne 'NA' && $lo_age{$junior} ne 'NA' && $lo_age{$qry_name} > $lo_age{$junior} )
    #                    {
    #                        $whos_the_child = 'query' ;
    #                    }
    #                    elsif ( $lo_age{$qry_name} eq 'NA' || $lo_age{$junior} eq 'NA' || $lo_age{$qry_name} == $lo_age{$junior} )
    #                    {
    #                        if ( $len_r > $len_junior ){
    #                            $whos_the_child = 'subject' ;
    #                        }
    #                        elsif( $len_r < $len_junior )
    #                        {
    #                            $whos_the_child = 'query' ;
    #                        }
    #                        elsif( $len_r == $len_junior ){
    #                            if ( index($qry_name,"_") != -1 && index($junior,"_") != -1 ){
    #                                $whos_the_child = 'subject' ;
    #                            }
    #                            elsif ( index($qry_name,"_") != -1 ){
    #                                $whos_the_child = 'query' ;
    #                            }
    #                            elsif ( index($junior,"_") != -1 ){
    #                                $whos_the_child = 'subject' ;
    #                            }
    #                            else
    #                            {
    #                                $whos_the_child = 'query' ; #random ..
    #                            }
    #                        }
    #                    }
    #
    #                    # resolve age of the pair - use alignment of all group1 vs group2
    #                    my $age = 0 ; my $count = 0 ; my $count_qu = 0 ;
    #                    my @group_query = @{$groups{$qry_name}} ;
    #                    my @group_subject = @{$groups{$junior}} ;
    #                    #my $max_n_align = 50 ; #max(5, int(20 - 0.05*$size_subfam)) ;
    #                    #print "max alignment per shuffle : $max_n_align\n" ;
    #                    #print "comparing @group_query with @group_subject \n" ;
    #                    print "group_query_n : $#group_query , group_subject_n : $#group_subject \n" ; 
    #                    for(my $i = 0; $i < $max_n_align; $i++){
    #                        #foreach my $gr_qu ( shuffle(@group_query) ){
    #                        my $rand1 = int(rand($#group_query+1)) ;
    #                        my $gr_qu = $group_query[$rand1] ;
    #                        for(my $j = 0; $j < $max_n_align; $j++ ){
    #                            my $rand2 = int(rand($#group_subject+1)) ;
    #                            my $gr_su = $group_subject[$rand2] ;
    #                            #foreach my $gr_su ( shuffle(@group_subject) ){
    #                            if ( $all_pairs{$gr_qu}->{$gr_su} ne 'NA' ){
    #                                $age += $all_pairs{$gr_qu}->{$gr_su} ;
    #                                $count++ ;
    #                                #if ( $count > $max_n_align ){ last }
    #                            }
    #                            if ( $gr_qu eq $gr_su ){ print "\n\n\ndiscrepency : $gr_qu appear in 2 groups !!\n\n\n" } 
    #                        }
    #                        #$count_qu++ ;
    #                        #if ( $count_qu > $max_n_align ){ last }
    #                    }
    #                    $age = $age / $count ;
    #                    # resolve groups
    #                    my @union = unique(@group_query, @group_subject);
    #                    $groups{$qry_name} = [ @union ] ; # fuse groups
    #                    $groups{$junior} = [ @union ] ;
    #                    #print "age resolved\n" ;
    #
    #                    # save results after decision 
    #                    if ( $whos_the_child eq 'query' ){
    #                        $age{$qry_name} = $age ;
    #                        $status{$qry_name} = 'child' ; $status{$junior} = 'parent' ;
    #                        $progeny{$junior}->{$qry_name} = $youngest_age{$qry_name} ;
    #                        $ancestry{$qry_name} = $junior ; 
    #                    }
    #                    elsif ( $whos_the_child eq 'subject' ){
    #                        $age{$junior} = $age ;
    #                        $status{$junior} = 'child' ; $status{$qry_name} = 'parent' ; 
    #                        $progeny{$qry_name}->{$junior} = $youngest_age{$junior} ;
    #                        $ancestry{$junior} = $qry_name ;
    #                    }
    #                    else
    #                    {
    #                        print "PROBLEM with line $i : who's the child is $whos_the_child\n" ;
    #                        print "Query: $qry_name\t$status{$qry_name}\t$youngest{$qry_name}\t$youngest_age{$qry_name}\n" ;
    #                    }
    #                    $resolved_child++ ; 
    #                    if ( $status_junior eq 'parent' && $status{$qry_name} eq 'parent' ){ $resolved_parent-- }
    #                    elsif ( $status_junior ne 'parent' && $status{$qry_name} ne 'parent' ){ $resolved_parent++ }
    #                }
    #            }
    #            if ( $magic_pass eq 'T' ){ $magic_pass = 'F' }
    #
    #            my $pass = $x + 1 ;
    #            my $n_parent = 0 ;
    #            foreach ( keys %status ){
    #                if ( $status{$_} eq 'parent' ){ $n_parent++ }
    #            }
    #            print "Pass$pass, resolved child: $resolved_child, parents left: $n_parent\n";
    #            if ( $resolved_child eq $prev_n_child ){ 
    #                $round_above_limit++ ;
    #
    #                # force pairing even if not reciprocal 
    #                $magic_pass = 'T' ;
    #
    #                if ( $round_above_limit > $max_iter  ){
    #                    print "Converged at pass$pass\nExiting...\n" ;
    #                    last ;
    #                } 
    #            }
    #            else 
    #            {
    #                $prev_n_child = $resolved_child ;
    #                $prev_n_parent = $resolved_parent ;
    #            } 
    #            if ( $n_parent < 5 ){
    #                print "Converged at pass$pass\nExiting...\n" ;
    #                last ;
    #            }
    #        }
    #
    #        ### PRINTING OUT RESULTS ###
    #        my $out_file = "$out_dir/$subfam_name.KimIntAge.txt" ;
    #        open my $out, '>', $out_file or die "can't open $out_file: $!" ;
    #        foreach ( @all_names ){
    #            my $age = 'NA' ;
    #            if ( exists $age{$_} ){
    #                $age = $age{$_} ;
    #            }
    #            elsif ( exists $youngest_age{$_} )
    #            {
    #                $age = $youngest_age{$_} ;
    #            }
    #            print $out "$summary{$_}\t$age\t$status{$_}\n" ;
    #        }
    #        close $out ;
    #    }
}

# liftover_age subs
sub mean_length {
    my $list_file = $_[0] ;
    open my $tab_list, '<', $list_file or die "can't open $list_file: $!" ;
    my $total = 0 ;
    my $n = 0 ;
    while(<$tab_list>){
        chomp ;
        my @line = split /\t/ ;
        my $length = $line[2] - $line[1] ;
        $total += $length ;
        $n++ ;
    }
    close $tab_list ;
    my $mean = $total / $n ;
    return $mean ;
}

sub get_lo_age {
    my $age = 'NA';
    switch($_[0])
    {
        case [1] { $age = 0.0 }
        case [2] { $age = 6.7 }
        case [3] { $age = 9.1 }
        case [4] { $age = 15.8 }
        case [5] { $age = 20.2 }
        case [6,7] { $age = 29.4 }
        case [8] { $age = 43.2 }
        case [9] { $age = 67.1 }
        case [10,11] { $age = 74.0 }
        case [12] { $age = 82.0 }
        case [13,14,15] { $age = 90.0 }
        case [16,17,18,19,20] { $age = 96.0 }
        case [21,22,23] { $age = 105 }
        case [24,25] { $age = 159 }
        case [26] { $age = 177 }
        case [27,28] { $age = 312 }
        case [29] { $age = 413 }
        case [30] { $age = 435 }
    }
    $age;
}


# merge_age subs
#sub get_length {
#    my $guy = $_[0] ;
#    #print $guy ;
#    my @guys = split(/[:-]/,$guy) ;
#    #print @guys ;
#    my $len = $guys[2] - $guys[1] ;
#    return $len ;
#}
#
## shuffle a list in place using a pointer to it
#sub fisher_yates_shuffle {
#    my $deck = shift;  # $deck is a reference to an ar
#    my $i = @$deck;
#    while ($i--) {
#        my $j = int rand ($i+1);
#        @$deck[$i,$j] = @$deck[$j,$i];
#    }
#}
#
