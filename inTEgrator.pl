#!/usr/bin/perl

use strict;
use warnings;
use Sort::Fields;
use Getopt::Long;
use Data::Dumper qw(Dumper);
use List::Util qw[min max shuffle];
use Array::Utils qw(:all);

my $header = qq{
        __________              __          
  (_)__/_  __/ __/__ ________ _/ /____  ____
 / / _ \\/ / / _// _ `/ __/ _ `/ __/ _ \\/ __/
/_/_//_/_/ /___/\\_, /_/  \\_,_/\\__/\\___/_/   
               /___/                      v0.1            

inTEgrator by A.Coudray (LVG-EPFL 2019)
} ;

my $usage_general = qq{
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
my $dir = '' ;
my $suffix = '' ;
my $out_dir = '../results/default_resolve_age_2' ;
my $max_iter = 1000 ;
my $prefix = '' ;
my $max_n_align = 50 ;
my $subfam_subset = 'empty' ;
my $summary_dir = '../data/all_summaries' ; 

GetOptions(
    "help" => \$help,
    "dir=s" => \$dir,
    "suffix=s" => \$suffix,
    "out_dir=s" => \$out_dir,
    "max_iter=s" => \$max_iter,
    "prefix=s" => \$prefix,
    "max_n_align=i" => \$max_n_align,
    "subfam_subset=s" => \$subfam_subset,
    "summary_dir=s" => \$summary_dir,
);

# Print Help and exit
if ($help) {
    print "$header\n$usage_general";
    exit(0);
}

`mkdir -p $out_dir` ;

# Record subfam_subset
my %sfam_todo ;
open my $sfam_subset, '<', $subfam_subset or die "can't open $subfam_subset: $!" ;
while(<$sfam_subset>){
    chomp;
    $sfam_todo{$_} = 1 ;
}
close $sfam_subset ;


my @files = <$dir/$prefix*$suffix>;
foreach my $file (@files) {
    
    my $subfam_name = $file ; # get sample name - basename without suffix
    $subfam_name =~ s/($dir|$suffix|\/)//g ;

    if ( exists $sfam_todo{$subfam_name} ){
        print "Launching analysis for $subfam_name\n";
    }
    else
    {
        print "skipping $subfam_name, not in subset...\n" ;
        next ;
    }
    
    ### Put summary file with cluster and LO age information 
    print "getting summary files...\n" ;
    my $summary_file = "$summary_dir/$subfam_name.summary.txt" ;
    if ( ! -e $summary_file ){
        print "no summary file for $subfam_name, skipping...\n" ;
        next ;
    }
    my %clusters ; my %lo_age ; my %summary ;
    open my $summ, '<', $summary_file or die "can't open $summary_file: $!" ;
    while(<$summ>){
        chomp ;
        my @line = split /\t/ ;
        my $key = $line[0] ;
        $lo_age{$key} = $line[5] ;
        $summary{$key} = $_ ;
    }
    close $summ ;

    ### STEP 1 : Get youngest age for each integrants and get pairs ###
    my %all_pairs ; my %groups ;
    my %youngest ; my @all_names ; my $first_line = 'Y' ; my $line_count = -1 ; 
    my %youngest_age ; my %idx_to_name ; my %name_to_idx ; my %status ; 
    my $mean_age = 0 ; my $total_count = 0 ;
    # First pass through matrix to get pairs
    open my $tab, '<', $file or die "cant open $file: $!";
    while (<$tab>){
        chomp;
        my @line = split / / ;
        my $minimum = 1000 ; my $field_count = -1 ;
        my $name = '' ; my $first = 'Y' ;
        foreach ( @line ){
            if ( $_ eq 'NA' ){ next } # ignore NAs
            if ( $first eq 'F' && $first_line eq 'F' ){ 
                $mean_age += $_ ; $total_count++ ;
                $all_pairs{$name}->{$all_names[$field_count]} = $_ ;
            }
            if ( $first_line eq 'Y' && $first eq 'F' ){ # record header
                push @all_names, $_ ;
                my @seq_name = ($_) ;
                $groups{$_} = [ @seq_name ] ; # initilize groups, all seq has itself
                $status{$_} = 'NA' ;
            }
            elsif ( $first eq 'Y' ){ ## avoid first field
                $name = $_ ; $first = 'F' ;
            }
            elsif ( $_ < $minimum && $name ne $all_names[$field_count] ){ 
                $minimum = $_ ;
                $youngest{$name} = $all_names[$field_count] ; 
                $youngest_age{$name} = $_ ;
            }
            $name_to_idx{$name} = $field_count ;
            $idx_to_name{$field_count} = $name ;
            $field_count++ ;
            
        }
        $first_line = 'F' ; $line_count++ ;
    }
    close $tab;
    
    $mean_age = $mean_age / $total_count ;
    my $NA_replacement = $mean_age ;
    my $size_subfam = $#all_names + 1 ;

    ### STEP 2 : GET SEED WITH BEST MUTUAL MATCH ### 
    my %age ; my %progeny ; my %ancestry ; 
    my $resolved_child = 0 ; my $resolved_parent = 0 ;
    for (my $i = 0; $i <= $#all_names; $i ++) {
        my $qry_name = $all_names[$i] ; 
        if ( $qry_name eq 'names' ){ next }
        my $junior = $youngest{$qry_name} ;
        
        if ( $youngest{$junior} eq $qry_name ){ 
            my $whos_the_child = 'NA' ;
            my $len_r = get_length($qry_name) ;
            my $len_junior = get_length($junior) ;
            
            # Make test to determine who's the child/parent sequence
            if ( $len_r < 0.5*$len_junior ){
                $whos_the_child = 'query' ;
            }
            elsif ( $len_junior < 0.5*$len_r ){
                $whos_the_child = 'subject' ;
            }
            elsif ( $lo_age{$qry_name} ne 'NA' && $lo_age{$junior} ne 'NA' && $lo_age{$qry_name} < $lo_age{$junior} )
            {
                $whos_the_child = 'query' ;
            }
            elsif ( $lo_age{$qry_name} ne 'NA' && $lo_age{$junior} ne 'NA' && $lo_age{$qry_name} > $lo_age{$junior} )
            {
                $whos_the_child = 'subject' ;
            }
            elsif ( $lo_age{$qry_name} eq 'NA' || $lo_age{$junior} eq 'NA' || $lo_age{$qry_name} == $lo_age{$junior} )
            {
                if ( $len_r > $len_junior ){
                    $whos_the_child = 'subject' ;
                }
                elsif( $len_r < $len_junior )
                {
                    $whos_the_child = 'query' ;
                }
                elsif( $len_r == $len_junior ){
                    if ( index($qry_name,"_") != -1 && index($junior,"_") != -1 ){
                        $whos_the_child = 'subject' ;
                    }
                    elsif ( index($qry_name,"_") != -1 ){
                        $whos_the_child = 'query' ;
                    }
                    elsif ( index($junior,"_") != -1 ){
                        $whos_the_child = 'subject' ;
                    }
                    else
                    {
                        $whos_the_child = 'query' ;
                    }
                }
            }
            if ( $whos_the_child eq 'query' ){
                $age{$qry_name} = $youngest_age{$qry_name} ;
                $status{$qry_name} = 'child' ; $status{$junior} = 'parent' ;
                $progeny{$junior}->{$qry_name} = $youngest_age{$qry_name} ;
                $ancestry{$qry_name} = $junior ;
                my @union = ($qry_name, $junior) ;
                $groups{$qry_name} = [ @union ] ; 
                $groups{$junior} = [ @union ] ;
            }
            elsif ( $whos_the_child eq 'subject' ){
                $age{$junior} = $youngest_age{$junior} ;
                $status{$junior} = 'child' ; $status{$qry_name} = 'parent' ; 
                $progeny{$qry_name}->{$junior} = $youngest_age{$junior} ;
                $ancestry{$junior} = $qry_name ;
                my @union = ($qry_name, $junior) ;
                $groups{$qry_name} = [ @union ] ; 
                $groups{$junior} = [ @union ] ;
            }
            else
            {
                print "PROBLEM with line $i : who's the child is $whos_the_child\n" ;
            }
            $resolved_child++ ; $resolved_parent++ ;
        }
    }
    
    print "Pass1, resolved child: $resolved_child, resolved parents: $resolved_parent\n"; 
    

    
    ### STEP3 : ITERATE N TIMES OVER THE MATRIX TO RESOLVE CHILD/PARENTS RELATION ### 
    
    my $prev_n_child = $resolved_child ; my $prev_n_parent = $resolved_parent ;
    my $round_above_limit = 0 ; my $magic_pass = 'F' ;
    # Iterate N times until convergence  
    for( my $x = 1 ; $x <= $max_iter; $x++ ){

        # Re-establish youngest 
        $first_line = 'T' ;
        open my $tab2, '<', $file or die "cant open $file: $!";
        while (<$tab2>){
            chomp;
            my @line = split / / ;
            if ($first_line eq 'T' ){ $first_line = 'F' ; next }
            if ( $status{$line[0]} eq 'child' ){ next }
            my $minimum = 10000 ; my $field_count = -1 ;
            my $name = $line[0] ; ## name of the query TE 
            for ( my $i = 1; $i <= $#line; $i++ ) { 
                if ( $line[$i] eq 'NA' ){ next }#$line[$i] = $NA_replacement } # ignore NAs
                my $name_comp = $all_names[($i-1)] ; #name of the subject TE for comparison

                if ( $line[$i] < $minimum && $name ne $name_comp && $status{$name_comp} ne 'child' ){
                    $minimum = $line[$i] ;
                    $youngest{$name} = $name_comp ; 
                    $youngest_age{$name} = $line[$i] ;
                }
            }
            $line_count++ ;
        }
        close $tab2;

        for (my $i = 0; $i <= $#all_names; $i ++) {
            my $qry_name = $all_names[$i] ; 
            if ( $qry_name eq 'names' ){ next }
            if ( $status{$qry_name} eq 'child' ){ next }
            my $junior = $youngest{$qry_name} ;
            my $status_junior = $status{$junior} ;
            my $trial_to_resolve = 0 ;
            while( $status_junior eq 'child' ){ 
                $trial_to_resolve++ ;
                $junior = $ancestry{$junior} ;
                $status_junior = $status{$junior} ;
                if ( $status_junior eq 'child' ){
                    #print "junior had child status, his direct ancestry has $status_junior status!!\n" ;
                    #print "trial number $trial_to_resolve\n" ;
                } 
            }
            ##print "$qry_name,$junior\n" ;
            if ( $youngest{$junior} eq $qry_name || $magic_pass eq 'T' ){ 

                my $whos_the_child = 'NA' ;
                my $len_r = get_length($qry_name) ;
                my $len_junior = get_length($junior) ;

                ### Make test to determine who's and child and who's the parent sequence
                if ( $len_r < 0.5*$len_junior ){
                    $whos_the_child = 'query' ;
                }
                elsif ( $len_junior < 0.5*$len_r ){
                    $whos_the_child = 'subject' ;
                }
                elsif ( $lo_age{$qry_name} ne 'NA' && $lo_age{$junior} ne 'NA' && $lo_age{$qry_name} < $lo_age{$junior} )
                {
                    $whos_the_child = 'query' ;
                }
                elsif ( $lo_age{$qry_name} ne 'NA' && $lo_age{$junior} ne 'NA' && $lo_age{$qry_name} > $lo_age{$junior} )
                {
                    $whos_the_child = 'query' ;
                }
                elsif ( $lo_age{$qry_name} eq 'NA' || $lo_age{$junior} eq 'NA' || $lo_age{$qry_name} == $lo_age{$junior} )
                {
                    if ( $len_r > $len_junior ){
                        $whos_the_child = 'subject' ;
                    }
                    elsif( $len_r < $len_junior )
                    {
                        $whos_the_child = 'query' ;
                    }
                    elsif( $len_r == $len_junior ){
                        if ( index($qry_name,"_") != -1 && index($junior,"_") != -1 ){
                            $whos_the_child = 'subject' ;
                        }
                        elsif ( index($qry_name,"_") != -1 ){
                            $whos_the_child = 'query' ;
                        }
                        elsif ( index($junior,"_") != -1 ){
                            $whos_the_child = 'subject' ;
                        }
                        else
                        {
                            $whos_the_child = 'query' ; #random ..
                        }
                    }
                }
                
                # resolve age of the pair - use alignment of all group1 vs group2
                my $age = 0 ; my $count = 0 ; my $count_qu = 0 ;
                my @group_query = @{$groups{$qry_name}} ;
                my @group_subject = @{$groups{$junior}} ;
                #my $max_n_align = 50 ; #max(5, int(20 - 0.05*$size_subfam)) ;
                #print "max alignment per shuffle : $max_n_align\n" ;
                #print "comparing @group_query with @group_subject \n" ;
                print "group_query_n : $#group_query , group_subject_n : $#group_subject \n" ; 
                for(my $i = 0; $i < $max_n_align; $i++){
                #foreach my $gr_qu ( shuffle(@group_query) ){
                    my $rand1 = int(rand($#group_query+1)) ;
                    my $gr_qu = $group_query[$rand1] ;
                    for(my $j = 0; $j < $max_n_align; $j++ ){
                        my $rand2 = int(rand($#group_subject+1)) ;
                        my $gr_su = $group_subject[$rand2] ;
                        #foreach my $gr_su ( shuffle(@group_subject) ){
                        if ( $all_pairs{$gr_qu}->{$gr_su} ne 'NA' ){
                            $age += $all_pairs{$gr_qu}->{$gr_su} ;
                            $count++ ;
                            #if ( $count > $max_n_align ){ last }
                        }
                        if ( $gr_qu eq $gr_su ){ print "\n\n\ndiscrepency : $gr_qu appear in 2 groups !!\n\n\n" } 
                    }
                    #$count_qu++ ;
                    #if ( $count_qu > $max_n_align ){ last }
                }
                $age = $age / $count ;
                # resolve groups
                my @union = unique(@group_query, @group_subject);
                $groups{$qry_name} = [ @union ] ; # fuse groups
                $groups{$junior} = [ @union ] ;
                #print "age resolved\n" ;

                # save results after decision 
                if ( $whos_the_child eq 'query' ){
                    $age{$qry_name} = $age ;
                    $status{$qry_name} = 'child' ; $status{$junior} = 'parent' ;
                    $progeny{$junior}->{$qry_name} = $youngest_age{$qry_name} ;
                    $ancestry{$qry_name} = $junior ; 
                }
                elsif ( $whos_the_child eq 'subject' ){
                    $age{$junior} = $age ;
                    $status{$junior} = 'child' ; $status{$qry_name} = 'parent' ; 
                    $progeny{$qry_name}->{$junior} = $youngest_age{$junior} ;
                    $ancestry{$junior} = $qry_name ;
                }
                else
                {
                    print "PROBLEM with line $i : who's the child is $whos_the_child\n" ;
                    print "Query: $qry_name\t$status{$qry_name}\t$youngest{$qry_name}\t$youngest_age{$qry_name}\n" ;
                }
                $resolved_child++ ; 
                if ( $status_junior eq 'parent' && $status{$qry_name} eq 'parent' ){ $resolved_parent-- }
                elsif ( $status_junior ne 'parent' && $status{$qry_name} ne 'parent' ){ $resolved_parent++ }
            }
        }
        if ( $magic_pass eq 'T' ){ $magic_pass = 'F' }
        
        my $pass = $x + 1 ;
        my $n_parent = 0 ;
        foreach ( keys %status ){
            if ( $status{$_} eq 'parent' ){ $n_parent++ }
        }
        print "Pass$pass, resolved child: $resolved_child, parents left: $n_parent\n";
        if ( $resolved_child eq $prev_n_child ){ 
            $round_above_limit++ ;

            # force pairing even if not reciprocal 
            $magic_pass = 'T' ;
            
            if ( $round_above_limit > $max_iter  ){
                print "Converged at pass$pass\nExiting...\n" ;
                last ;
            } 
        }
        else 
        {
            $prev_n_child = $resolved_child ;
            $prev_n_parent = $resolved_parent ;
        } 
        if ( $n_parent < 5 ){
            print "Converged at pass$pass\nExiting...\n" ;
            last ;
        }
    }
    
    ### PRINTING OUT RESULTS ###
    my $out_file = "$out_dir/$subfam_name.KimIntAge.txt" ;
    open my $out, '>', $out_file or die "can't open $out_file: $!" ;
    foreach ( @all_names ){
        my $age = 'NA' ;
        if ( exists $age{$_} ){
            $age = $age{$_} ;
        }
        elsif ( exists $youngest_age{$_} )
        {
            $age = $youngest_age{$_} ;
        }
        print $out "$summary{$_}\t$age\t$status{$_}\n" ;
    }
    close $out ;
}

sub get_length {
    my $guy = $_[0] ;
    #print $guy ;
    my @guys = split(/[:-]/,$guy) ;
    #print @guys ;
    my $len = $guys[2] - $guys[1] ;
    return $len ;
}

# shuffle a list in place using a pointer to it
sub fisher_yates_shuffle {
    my $deck = shift;  # $deck is a reference to an ar
    my $i = @$deck;
    while ($i--) {
        my $j = int rand ($i+1);
        @$deck[$i,$j] = @$deck[$j,$i];
    }
}

