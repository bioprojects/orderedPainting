#! /usr/bin/perl

use Data::Dumper;
use Getopt::Std;
use FindBin;
#use Tie::IxHash;
use strict;
use warnings;

use vars qw ($opt_d $opt_l $opt_t $opt_s $opt_c $opt_n $opt_r $opt_m $opt_g $opt_v);
&getopts('d:l:t:s:c:n:rm:g:v');

my $usage = <<_EOH_;

#
# [-g file.hap]
#
#  -d prefix_rnd_1-10_both_dirs.list 
#    (list of ordered dirs where copyprobsperlocus.cat.sort.gz files are stored)
#
#  -l prefix_rnd_1-10_both_dirs_strainOrder.list
#
#  -s prefix_fine_forward.strainOrder 
#
#
# -n 1 or 2 or ... (to conduct only n-th ordering in the "-d" file, for parallel computing) 
#    or
# -r               (final step)
#
# [-m pos2missingInd.txt]
#

_EOH_
;

###########################################################################################################

#
# env
#
my $env_file = "$FindBin::Bin/env_func.bashrc";

#
# path of commands/scripts
#
my $sort_path = "$FindBin::Bin/sort";
my $sort_opt = " -m -n --batch-size=300 --parallel=8"; # "-m" makes the sorting much faster when all input files are sorted (ascending order)
                                                       # --batch-size=100 is enough for 10 orderings
                                                       # --batch-size=300 is required for more than 50 orderings

my $postprocess_path = "$FindBin::Bin/postprocess/pp";

#
# output to each ordering dir
#
#   $out_each_dir_averave_matrix
#
my $out_each_dir_averave_matrix  = "average.matrix.txt";
my $out_each_dir_site_distScore  = "site_distScore.txt";

my $sortfile_catOrderings_copyprob    = "copyprobsperlocus.cat.sort";
my $gz_sortfile_catOrderings_copyprob = "copyprobsperlocus.cat.sort.gz";


my $out_each_dir_site_minus_average_matrix_summary = "site_minus_average.matrix.txt"; # created by -r

#
# output to results dir
#
my $out_results             = "results_siteStats.txt"; # gzipped at the end 
my $out_results_summary_pos = "results_siteStats_summary.pos.txt";

my $out_sum_site_minus_average_summary        = "sum_site_minus_average.summary.txt"; # gzipped at the end 
my $out_sum_site_minus_average_summary_range  = "sum_site_minus_average.summary.range.txt";

#
# number of sites in sum_site_minus_average.summary.txt and for visualization
#
my $prop_top_sites_summary           = 0.01; # proportion of top sites to be written in sum_site_minus_average.summary.txt 
#my $num_top_sites_visualization      =  200;

my $num_other_sites_summary           = 10;  # num of middle/botton sites to be written in sum_site_minus_average.summary.txt 
my $num_other_sites_visualization     = 10;


#
# IN
#
my $dir_ordering_listFile = $opt_d or die $usage;
my $strainHapOrder_listFile = $opt_l or die $usage;

my $hap_file = "";
if ($opt_g) {
  $hap_file = $opt_g;
} else {
  $hap_file = $dir_ordering_listFile;
  $hap_file =~ s/_orderedS[0-9].*\.list$/.hap/g; 
  if (! -s $hap_file) {
    die "Error: $hap_file doesn't exist or empty";
  }
}

my $type_painting = 2; #  -t 2 (ordering) | 1 (all against everyone else, not supported in orderedPainting.sh) 
if ($opt_t) {
  $type_painting = $opt_t;
}

my $contrast_min = "";
my $contrast_max = "";
#if ($opt_c) { 
#  $contrast_min  = 0;
#  $contrast_max  = $opt_c;
#  if (!($contrast_max > 0)) {
#    die "Error: positive values must be specified by -c ";
#  }
#  $out_each_dir_site_distScore =~ s/\.txt/_c$contrast_max.txt/g; 
#  $out_results =~ s/\.txt/_c$contrast_max.txt/g;
#} else {
  $contrast_min  = 1;
  $contrast_max  = 1;
#}

if ($type_painting == 1) {
  $out_each_dir_site_distScore =~ s/\.txt$/_t1.txt/g; 
  $out_results =~ s/\.txt$/_t1.txt/g;
}

if (! -s $dir_ordering_listFile) {
  die "Error: $dir_ordering_listFile doesn't exist or empty";
}
if (! -s $strainHapOrder_listFile) {
  die "Error: $strainHapOrder_listFile doesn't exist or empty";
}

my @arr_ind_fineOrdering = ();
my %hash_ind_fineOrdering = ();
my $strainFineOrderFile = $opt_s or die $usage; # if this is not specified, site_contrast_distScore.txt will not be created
if (! -s $strainFineOrderFile) {
  die "Error: $strainFineOrderFile doesn't exist or empty";
}
open(FINESTRUCT_ORDER, $strainFineOrderFile);
while (my $line = <FINESTRUCT_ORDER>) {
  chomp($line);
  my @arr_line = split(/\t/, $line);
  
  my $strainName = "";
  if (scalar(@arr_line) == 1) {
    $strainName = $arr_line[0];
  } else {
    $strainName = $arr_line[1];
  }
  
  push(@arr_ind_fineOrdering, $strainName);
  $hash_ind_fineOrdering{$strainName} = 1;
}
close(FINESTRUCT_ORDER);

my $out_fine_header = "";
foreach my $header_name (@arr_ind_fineOrdering) {
  $out_fine_header .= "$header_name ";
}
$out_fine_header =~ s/ $//g;

my $out_dir_results = $dir_ordering_listFile;
   $out_dir_results =~ s/\.list$/_results/g;

#
# global variables
#
my $cmd = "";
my $cmd_common;
my $loop_part;
my $stamp = "";

my $num_site = `head -3 $hap_file | tail -1`;
chomp($num_site);

my $num_ind = `head -2 $hap_file | tail -1`;
chomp($num_ind);

my $num_dir_orderings = `wc -l $dir_ordering_listFile | awk '{print \$1}'`;
chomp($num_dir_orderings);

# for restore and estimation
my %hash_sum_site_minus_ave = ();
my %hash_sum_site_distScore = ();
my %hash_sum_site_bootstrapped_distScore = ();
my %hash_sum_site_infoContent = ();

#
# qsub
#
my $QSUB = "";
my $QSTAT = "qstat";

my $QUEUE_TYPE = `grep QUEUE_TYPE $env_file | grep -v '#'`;
chomp($QUEUE_TYPE);
$QUEUE_TYPE =~ s/"//g;
$QUEUE_TYPE =~ s/QUEUE_TYPE=//g;
$QUEUE_TYPE =~ s/ //g;

if ($QUEUE_TYPE eq "SGE") {
  $QSUB = "qsub -cwd -N ";
} elsif ($QUEUE_TYPE eq "LSF") {
  $QSUB = "bsub -J ";
} else {
  die "Error: uknown QUEUE_TYPE $QUEUE_TYPE";
}

#
# file about missing data
#
my $pos2missingInd_File = "";
if ($opt_m) {
  $pos2missingInd_File = $opt_m;
}

#
# file about constraint of donor and recipient
#
my $constraint_File = "";
if ($opt_c) {
  $constraint_File = $opt_c;
}


####################################################################
# main - 1
#   for each ordering
####################################################################
if (!$opt_n && !$opt_r) {
  die "Please specify -n integer (to conduct only n-th ordering in the -d file, for parallel computing) or -r (to restore)";
}

if (!$opt_r) {

  my $cnt_in_ordering_listFile = 1;
  open(DIR_ORDERING, $dir_ordering_listFile);
  while (my $dir_each_ordering = <DIR_ORDERING>) {
    chomp($dir_each_ordering);
    $dir_each_ordering  =~ s/\/$//g;
    if (! -d $dir_each_ordering) {
      die "Error: $dir_each_ordering doesn't exist";
    }

    #
    # remove *.hap files 
    #   (if thery still remain because the pipeline is re-executed and some hap files are not re-painted)
    #
    my @arr_dothap_files = glob("$dir_each_ordering/*.hap");
    foreach my $each_dothap_file (@arr_dothap_files) {
      unlink($each_dothap_file);
    }

    #
    # start
    #
    $stamp = `date +%Y%m%d_%T`;
    chomp($stamp);

    if ($opt_n) {
      if ($cnt_in_ordering_listFile != $opt_n) {
        print("$stamp $dir_each_ordering skipped (because it is not specified by -n) \n");
        $cnt_in_ordering_listFile++;
        next;
      }
    }

    print("$stamp $dir_each_ordering started\n");

    #
    # prepare strainHapOrderFile
    #
    my $strainHapOrderFile = `head -$cnt_in_ordering_listFile $strainHapOrder_listFile | tail -1`;
    chomp($strainHapOrderFile);
    if (! -s $strainHapOrderFile) {
      die "Error: $strainHapOrderFile doesn't exist or empty";
    }

    #
    # prepare cmd
    #
    $cmd_common  = "gzip -dc $dir_each_ordering/$gz_sortfile_catOrderings_copyprob |";
    $cmd_common .= " $postprocess_path ";
    $cmd_common .= " -d $dir_each_ordering ";
    $cmd_common .= " -l $strainHapOrderFile ";
    $cmd_common .= " -o $strainFineOrderFile ";
    if ($opt_m) {
      $cmd_common .= " -m $pos2missingInd_File ";
    }
    $cmd_common .= " -r $out_dir_results ";
    $cmd_common .= " -t $type_painting ";

    ################################################################
    # calculate average matrix (long loop)
    ################################################################
    $loop_part = 1;
    
    my $nrow_ave_matrix = 0;
    if (-s "$dir_each_ordering/$out_each_dir_averave_matrix") {
      $nrow_ave_matrix = `wc -l $dir_each_ordering/$out_each_dir_averave_matrix | awk '{print \$1}'`;
      chomp($nrow_ave_matrix);
    }

    if ($nrow_ave_matrix == scalar(@arr_ind_fineOrdering)+1) {
      print "$dir_each_ordering/$out_each_dir_averave_matrix already exists. Skipped.\n";
    } else {
      #
      # calculate an average matrix by processing each site in $gz_sortfile_catOrderings_copyprob
      #
      $stamp = `date +%Y%m%d_%T`;
      chomp($stamp);
      print "$stamp preparation of average matrix of this ordering started\n";

      $cmd = $cmd_common . " -p $loop_part";
      print("$cmd\n");
      if( system("$cmd") != 0) { die("Error: $cmd failed"); };

      $stamp = `date +%Y%m%d_%T`;
      chomp($stamp);
      print "$stamp preparation of average matrix of this ordering ended\n";
    }

    #############################################################################
    # calculate the distance statistic and its bootstrappd samples
    # for each site (long loop)
    #############################################################################
    $loop_part = 2;
    
    $stamp = `date +%Y%m%d_%T`;
    chomp($stamp);
    print "$stamp Calculation of distance to the average  ( $dir_each_ordering/$out_each_dir_site_distScore ) started \n";

    my $nrow_distScore = 0;
    if (-s "$dir_each_ordering/$out_each_dir_site_distScore") {
      $nrow_distScore = `wc -l $dir_each_ordering/$out_each_dir_site_distScore | awk '{print \$1}'`;
      chomp($nrow_distScore);
    }

    if ($nrow_distScore == $num_site) {
      print "$dir_each_ordering/$out_each_dir_site_distScore already exists. Skipped.\n";
    } else {
      #
      # calculate distance statistic and its bootstrappd samples 
      # by processing each site in $gz_sortfile_catOrderings_copyprob
      #   (always using the same seed internally)
      #
      $stamp = `date +%Y%m%d_%T`;
      chomp($stamp);
      print "$stamp Calculation of distance statistic of this ordering started\n";

      $cmd = $cmd_common . " -p $loop_part";
      if ($opt_c) {
        $cmd .= " -c $constraint_File";
      }
      print("$cmd\n");
      if( system("$cmd") != 0) { die("Error: $cmd failed"); };

      $stamp = `date +%Y%m%d_%T`;
      chomp($stamp);
      print "$stamp Calculation of distance statistic of this ordering ended\n";
    }

    $cnt_in_ordering_listFile++;

  } # ordering loop 
  close(DIR_ORDERING);

} else {

  #
  # pre-check for restore
  #
  my $err_msg  = "";
  open(DIR_ORDERING, $dir_ordering_listFile);
  while (my $dir_each_ordering = <DIR_ORDERING>) {
    chomp($dir_each_ordering);
    $dir_each_ordering  =~ s/\/$//g;
    if (! -d $dir_each_ordering) {
      die "Error: $dir_each_ordering doesn't exist";
    }

    if (! -s     "$dir_each_ordering/$out_each_dir_site_distScore") { 
      $err_msg .= "Error: $dir_each_ordering/$out_each_dir_site_distScore doesn't exist or empty\n";
    }
  }
  close(DIR_ORDERING);
  
  if ($err_msg ne "") {
    die "$err_msg";
  }

}


####################################################################
# main - 2
#   across orderings
####################################################################
if ($opt_n) {
  
  print "Because -n option was used to specify an ordering to be processesed, this program stops here.\n";
  print "Please proceed to the final step by using -r option\n";
  exit(0);

} else {

  if (!-d $out_dir_results) {
    mkdir($out_dir_results);
  }

  my $header = "";
  my $top_threshold = "";
  my $N_bootstrap = "";

  if ($opt_r) {
    #
    # preparation
    #
    $stamp = `date +%Y%m%d_%T`;
    chomp($stamp);
    print("$stamp preparing calculation of distScore across the orderings ... \n");


    ###################################################################################################################
    # $out_dir_results/$out_results
    # $out_dir_results/$out_results_summary_pos
    ###################################################################################################################

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # calculate summation of statistics across orderings
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    #
    # cat and sort $out_each_dir_site_distScore across orderings
    #
    my $cmd_sort1  = "$sort_path $sort_opt "; 
       $cmd_sort1 .= " -T $out_dir_results ";
    open(DIR_ORDERING, $dir_ordering_listFile);
    while (my $dir_each_ordering = <DIR_ORDERING>) {
      print("$dir_each_ordering");
      chomp($dir_each_ordering);
      $dir_each_ordering  =~ s/\/$//g;
      if (! -d $dir_each_ordering) {
        die "Error: $dir_each_ordering doesn't exist";
      }
      $cmd_sort1 .= " $dir_each_ordering/$out_each_dir_site_distScore ";
    }
    close(DIR_ORDERING);

    #
    # %hash_sum_site_distScore (*pos*=>contrast=>type=>value)
    #
    $stamp = `date +%Y%m%d_%T`;
    chomp($stamp);
    print("$stamp calculating summation of distScore across the orderings ... \n");

    print "$cmd_sort1\n";
    open(SITE_DIST_CAT_SORT, "$cmd_sort1 |");
    while (my $line = <SITE_DIST_CAT_SORT>) {
      #if ($line !~ /^[0-9]/) { # not required
      #  next;
      #}
      chomp $line;

      my @arr_line = split(/\t/, $line);
      my $pos = $arr_line[0];

      my $contrast = 1;
      my $distScore_per_mat = $arr_line[1];
      #my $donorInfoContent = $arr_line[2];

      if ($contrast_max < $contrast) {
        next;
      }

      if (!defined($distScore_per_mat)) {
        print "Warning: distScore_per_mat is empty and skipped at pos=$pos\n";
      } else {
        if (!defined($hash_sum_site_distScore{$pos})) { 
          $hash_sum_site_distScore{$pos}  = $distScore_per_mat;
          
          # bootstrapped distScore
          foreach (my $i=2; $i<scalar(@arr_line); $i++) {
            $hash_sum_site_bootstrapped_distScore{$pos}{$i-1} = $arr_line[$i];
          }
        } else {
          $hash_sum_site_distScore{$pos} += $distScore_per_mat;
          
          # bootstrapped distScore
          foreach (my $i=2; $i<scalar(@arr_line); $i++) {
            $hash_sum_site_bootstrapped_distScore{$pos}{$i-1} += $arr_line[$i];
          }
        }
      }

      #if (!defined($donorInfoContent)) {
      #  print "Warning: donotInfoContent is empty and skipped at pos=$pos\n";
      #} else {
      #  if (!defined($hash_sum_site_infoContent{$pos})) { 
      #    $hash_sum_site_infoContent{$pos}  = $donorInfoContent;
      #  } else {
      #    $hash_sum_site_infoContent{$pos} += $donorInfoContent;
      #  }
      #}

    }
    close(SITE_DIST_CAT_SORT);

    $stamp = `date +%Y%m%d_%T`;
    chomp($stamp);
    print("$stamp calculating summation of distScore across the orderings ... finished \n");


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # extract position of summary sites in terms of the distance statistic summed across the orderings
    #  and
    # output statistics in a descending order
    #  (always overwrite)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    $stamp = `date +%Y%m%d_%T`;
    chomp($stamp);
    print("$stamp output $out_dir_results/$out_results_summary_pos ... \n");

    #
    # extract position of summary sites and record the threshold of "top"
    #
    my $prefix_of_dirs_for_visualization = "$out_dir_results/visualize_";
    
    my %hash_visType_count = ();
    my %hash_pos_visType = ();
    my @arr_type = ('top','middle','bottom');
    foreach my $each_type (@arr_type) {
      my $dir = $prefix_of_dirs_for_visualization . $each_type; 
      if (! -d $dir) {
        mkdir("$dir");
        print("$dir was created\n")
      }
      $hash_visType_count{$each_type} = 0;
    }

    my %hash_summaryPos2Type = ();
    my %hash_summaryPos_rank = ();

    my $i_site = 1; 
    my $num_top_sites_summary = sprintf("%0.f",$num_site*$prop_top_sites_summary);

    # for each position sorted by distScore (descending)
    foreach my $pos (sort { $hash_sum_site_distScore{$b} <=> $hash_sum_site_distScore{$a} } keys %hash_sum_site_distScore){
     
      my $type = "";
      if ( $i_site <  $num_top_sites_summary ) {
        $type = "top";
        $hash_summaryPos2Type{$pos} = "$type";
        $hash_summaryPos_rank{$pos} = "$i_site";
        $hash_visType_count{$type}++;
        
        #if ($hash_visType_count{$type} <= $num_top_sites_visualization) {
          $hash_pos_visType{$pos} = $type;
        #}
        
        #
        # record the threshold of "top"
        #
        if ($i_site == $num_top_sites_summary-1) {
          $top_threshold = $hash_sum_site_distScore{$pos};
        }
        
      } elsif (
           ($i_site > $num_site/2 - $num_other_sites_summary/2) && 
           ($i_site < $num_site/2 + $num_other_sites_summary/2) 
              ) {
        $type = "middle";
        $hash_summaryPos2Type{$pos} = "$type";
        $hash_summaryPos_rank{$pos} = "$i_site";
        $hash_visType_count{$type}++;
        if ($hash_visType_count{$type} <= $num_other_sites_visualization) {
          $hash_pos_visType{$pos} = $type;
        }
      } elsif ($i_site >  $num_site - $num_other_sites_summary ) {
        $type = "bottom";
        $hash_summaryPos2Type{$pos} = "$type";
        $hash_summaryPos_rank{$pos} = "$i_site";
        $hash_visType_count{$type}++;
        if ($hash_visType_count{$type} <= $num_other_sites_visualization) {
          $hash_pos_visType{$pos} = $type;
        }
      }

      $i_site++;
    }


    #
    # output
    #
    $stamp = `date +%Y%m%d_%T`;
    chomp($stamp);
    print("$stamp output $out_dir_results/$out_results and $out_dir_results/$out_results_summary_pos ... \n");
    
    open(OUT_RESULTS, "> $out_dir_results/$out_results");
    open(OUT_SUMMAY_POS, "> $out_dir_results/$out_results_summary_pos");

    print OUT_RESULTS "pos";
    print OUT_RESULTS "\t" . "D_i";
    #print OUT_RESULTS "\t" . "results_of_c";
    #print OUT_RESULTS "\t" . "sum_donorInfoContent";
    print OUT_RESULTS "\t" . "bootstrap";
    print OUT_RESULTS "\n";

    print OUT_SUMMAY_POS "pos";
    print OUT_SUMMAY_POS "\t" . "rank";
    print OUT_SUMMAY_POS "\t" . "type";
    print OUT_SUMMAY_POS "\t" . "D_i";
    print OUT_SUMMAY_POS "\t" . "bootstrap";
    print OUT_SUMMAY_POS "\n";

    # for each position sorted by distScore (descending)
    $i_site = 1;
    foreach my $pos (sort { $hash_sum_site_distScore{$b} <=> $hash_sum_site_distScore{$a} } keys %hash_sum_site_distScore){

      #my $contrast = 1;
      print OUT_RESULTS $pos;
      print OUT_RESULTS "\t".$hash_sum_site_distScore{$pos};
      #print OUT_RESULTS "\t".$contrast;
      #print OUT_RESULTS "\t".$hash_sum_site_infoContent{$pos};  # to be commented out
      
      #
      # bootstrap support by using the threshold of the top percentile 
      # of the distribution of Di for all sites without bootstrapping
      #
      my $bootstrap_support = "NA";
      #if ($hash_summaryPos2Type{$pos} eq "top") { # to output only for atypical sites
        my $cnt_bootstrapped_in_top = 0;
        
        my @arr_i_boot = keys %{$hash_sum_site_bootstrapped_distScore{$pos}};
        foreach my $i_boot (@arr_i_boot) {
          if ($top_threshold <= $hash_sum_site_bootstrapped_distScore{$pos}{$i_boot}) {
            $cnt_bootstrapped_in_top++;
          }
        }
        $bootstrap_support = ( $cnt_bootstrapped_in_top / scalar(@arr_i_boot) ) * 100;
      #}
      print OUT_RESULTS "\t" . $bootstrap_support;
      
      print OUT_RESULTS "\n";

      if (defined($hash_summaryPos2Type{$pos})) {
        print OUT_SUMMAY_POS  $pos;
        #print OUT_SUMMAY_POS "\t" . "dist" . $i_site;
        print OUT_SUMMAY_POS "\t" . $i_site;
        print OUT_SUMMAY_POS "\t" . $hash_summaryPos2Type{$pos};
        print OUT_SUMMAY_POS "\t" . $hash_sum_site_distScore{$pos};
        print OUT_SUMMAY_POS "\t" . $bootstrap_support;
        print OUT_SUMMAY_POS "\n";
      }

      $i_site++;
    }
    close(OUT_RESULTS);
    close(OUT_SUMMAY_POS);

    print("$stamp output $out_dir_results/$out_results and $out_dir_results/$out_results_summary_pos ... finished \n");


    $cmd = "gzip -f $out_dir_results/$out_results";
    if( system("$cmd") != 0) { die("Error: $cmd failed"); };

    if (-s "$out_dir_results/$out_results.gz") {
      $stamp = `date +%Y%m%d_%T`;
      chomp($stamp);
      print("$stamp gzipping to $out_dir_results/$out_results.gz finished\n");
    } else {
      print("$stamp gzipping to $out_dir_results/$out_results.gz failed\n");
    }



    ###################################################################################################################
    # output sum of Sij - Mj (long loop)
    #   only for the top and middle/bottom positions extracted above
    ###################################################################################################################
    $loop_part = 3;

    my $p3_job_name = `date +%d%H%S`;
    chomp($p3_job_name);
    $p3_job_name = "p3_" . $p3_job_name;

    $stamp = `date +%Y%m%d_%T`;
    chomp($stamp);
    print("$stamp output $out_each_dir_site_minus_average_matrix_summary for each ordering  ... \n");

    my $cnt_in_ordering_listFile = 0;
    open(DIR_ORDERING, $dir_ordering_listFile);
    while (my $dir_each_ordering = <DIR_ORDERING>) {
      print("$dir_each_ordering");
      chomp($dir_each_ordering);
      $dir_each_ordering  =~ s/\/$//g;
      if (! -d $dir_each_ordering) {
        die "Error: $dir_each_ordering doesn't exist";
      }

      $cnt_in_ordering_listFile++;
      my $strainHapOrderFile = `head -$cnt_in_ordering_listFile $strainHapOrder_listFile | tail -1`;
      chomp($strainHapOrderFile);
      if (! -s $strainHapOrderFile) {
        die "Error: $strainHapOrderFile doesn't exist or empty";
      }

      #
      # prepare cmd for this ordering
      #
      $cmd_common  = "gzip -dc $dir_each_ordering/$gz_sortfile_catOrderings_copyprob |";
      $cmd_common .= " $postprocess_path ";
      $cmd_common .= " -d $dir_each_ordering ";
      $cmd_common .= " -l $strainHapOrderFile ";
      $cmd_common .= " -o $strainFineOrderFile ";
      if ($opt_m) {
        $cmd_common .= " -m $pos2missingInd_File ";
      }
      $cmd_common .= " -r $out_dir_results ";
      $cmd_common .= " -t $type_painting ";

      $cmd = "$QSUB $p3_job_name <<< '$cmd_common -p $loop_part'";
      print("$cmd\n");
      system("$cmd");
    } # dir_each_ordering
    close(DIR_ORDERING);

    #
    # execute in parallel for the orderings
    # 
    while () {
      my $check = `$QSTAT | grep $p3_job_name | wc -l`;
      #print "$check";
      chomp($check);
      
      if ($check == 0) {
        last;
      }
      sleep 10;
    }
    
    my @arr_p3_job_logs = glob("$p3_job_name.*");
    foreach my $each_p3_job_log (@arr_p3_job_logs) {
      unlink($each_p3_job_log);
    }

    $stamp = `date +%Y%m%d_%T`;
    chomp($stamp);
    print("$stamp output $out_each_dir_site_minus_average_matrix_summary for each ordering  ... finished \n");


    # 
    # combine and output (overwrite)
    #   %hash_sum_site_minus_ave (pos=>recipient_name=>donor_name)
    #                                  * fine ordering * 
    #
    my $cmd_sort2  = "$sort_path $sort_opt "; 
       $cmd_sort2 .= " -T $out_dir_results ";
    open(DIR_ORDERING, $dir_ordering_listFile);
    while (my $dir_each_ordering = <DIR_ORDERING>) {
      print("$dir_each_ordering");
      chomp($dir_each_ordering);
      $dir_each_ordering  =~ s/\/$//g;
      if (! -d $dir_each_ordering) {
        die "Error: $dir_each_ordering doesn't exist";
      }
      $cmd_sort2 .= " $dir_each_ordering/$out_each_dir_site_minus_average_matrix_summary ";
    }
    close(DIR_ORDERING);


    my $i_ordering_of_this_pos_recipient = 0;
    open(OUT_SUM_DEV, "> $out_dir_results/$out_sum_site_minus_average_summary");
    print OUT_SUM_DEV "type distRankDesc pos $out_fine_header\n";

    $stamp = `date +%Y%m%d_%T`;
    chomp($stamp);
    print("$stamp output $out_dir_results/$out_sum_site_minus_average_summary ... \n");

    my $max = 0;
    my $min = 0;

    # note that for each pos, the number of lines = scalar(@arr_ind_fineOrdering) * $num_dir_orderings 
    print "$cmd_sort2\n";
    open(DIV_CAT_SORT, "$cmd_sort2 |");
    while (my $line = <DIV_CAT_SORT>) {
      if ($line !~ /^[0-9]/) {
        next;
      }
      chomp $line;

      $i_ordering_of_this_pos_recipient++;

      my @arr_line = split(/ /, $line);

      my $pos = $arr_line[0];
      my $recipient_name = $arr_line[1]; # 2nd column is required to distinguish rows with the same recipient

      for (my $i_donor=0; $i_donor<scalar(@arr_ind_fineOrdering); $i_donor++) { 
        my $donor_name     = $arr_ind_fineOrdering[$i_donor];
        
        if (!defined($hash_sum_site_minus_ave{$recipient_name}{$donor_name})) {
          $hash_sum_site_minus_ave{$recipient_name}{$donor_name}  = $arr_line[$i_donor+2]; # note: pos recipient_name values ... (num_dir_orderings lines)
        } else {
          $hash_sum_site_minus_ave{$recipient_name}{$donor_name} += $arr_line[$i_donor+2]; # note: pos recipient_name values ... (num_dir_orderings lines)
        }
      }

      #
      # after $hash_sum_site_minus_ave of this site was prepared above, 
      # output
      #
      
      #print "$pos\t$recipient_name\t$i_ordering_of_this_pos_recipient\t$num_dir_orderings\n";
      if ($i_ordering_of_this_pos_recipient != scalar(@arr_ind_fineOrdering)*$num_dir_orderings) {
        next;
      } else {
        my $type = $hash_summaryPos2Type{$pos};
        
        if (defined($hash_pos_visType{$pos})) {
          open(OUT_SUM_DEV_SITE, "> $prefix_of_dirs_for_visualization".$hash_summaryPos2Type{$pos}."/$pos.txt");
          print OUT_SUM_DEV_SITE "type distRankDesc pos $out_fine_header\n";
        }

        # output matrix of a site
        foreach my $recipient_name (@arr_ind_fineOrdering) {
          if (!defined($hash_summaryPos2Type{$pos})) {
            die "Error: $pos is not defined in hash_summaryPos2Type";
          }
          my $out_line_pos_recipient = $type . " " . $hash_summaryPos_rank{$pos} . " $pos ";
          foreach my $donor_name (@arr_ind_fineOrdering) {
            if (!defined($hash_sum_site_minus_ave{$recipient_name}{$donor_name})) {
              print Dumper(\%hash_sum_site_minus_ave);
              die "Error: hash_sum_site_minus_ave of pos=$pos, recipient=$recipient_name, donor=$donor_name is undefined";
            }
            $out_line_pos_recipient .= $hash_sum_site_minus_ave{$recipient_name}{$donor_name} . " ";
            
            if ($hash_sum_site_minus_ave{$recipient_name}{$donor_name} > $max) {
              $max = $hash_sum_site_minus_ave{$recipient_name}{$donor_name};
            }

            if ($hash_sum_site_minus_ave{$recipient_name}{$donor_name} < $min) {
              $min = $hash_sum_site_minus_ave{$recipient_name}{$donor_name};
            }
          } # donor
          
          $out_line_pos_recipient =~ s/ $//g;
          print OUT_SUM_DEV $out_line_pos_recipient . "\n";
          
          if (defined($hash_pos_visType{$pos})) {
            print OUT_SUM_DEV_SITE $out_line_pos_recipient . "\n";
          }
        } # recipient

        if (defined($hash_pos_visType{$pos})) {
          close(OUT_SUM_DEV_SITE);
        }

        #
        # for the next pos
        #
        $i_ordering_of_this_pos_recipient = 0;
        %hash_sum_site_minus_ave = ();
      }
    }
    close(OUT_SUM_DEV);

    $stamp = `date +%Y%m%d_%T`;
    chomp($stamp);
    print("$stamp output $out_dir_results/$out_sum_site_minus_average_summary ... finished \n");


    $cmd = "gzip -f $out_dir_results/$out_sum_site_minus_average_summary";
    if( system("$cmd") != 0) { die("Error: $cmd failed"); };
    
    $stamp = `date +%Y%m%d_%T`;
    chomp($stamp);
    if (-s "$out_dir_results/$out_sum_site_minus_average_summary.gz") {
      print("$stamp gzipping to $out_dir_results/$out_sum_site_minus_average_summary.gz finished \n");
    } else {
      print("Error: gzipping to $out_dir_results/$out_sum_site_minus_average_summary.gz failed");
    }

    $cmd = "echo '$min $max' > $out_dir_results/$out_sum_site_minus_average_summary_range";
    print("$cmd\n");
    if( system("$cmd") != 0) { die("Error: $cmd failed"); };

    $stamp = `date +%Y%m%d_%T`;
    chomp($stamp);
    print("$stamp done. \n");

  } # opt_r
}

