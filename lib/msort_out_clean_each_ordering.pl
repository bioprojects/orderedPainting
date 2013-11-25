#! /usr/bin/perl

use Data::Dumper;
use Getopt::Std;
use FindBin;
#use Tie::IxHash;
use strict;
use warnings;

use vars qw ($opt_d $opt_t $opt_g $opt_v);
&getopts('d:t:g:uv');

my $usage = <<_EOH_;

#
# -d dir_each_ordering
# -g phasefile.hap
#
# [-t 1]
#

_EOH_
;


#
# IN
#

my $type_painting = 2; #  -t 2 (ordering) | 1 (all against everyone else, not supported in orderedPainting.sh) 
if ($opt_t) {
  $type_painting = $opt_t;
}

my $dir_each_ordering   = $opt_d or die $usage;
my $phasefile           = $opt_g or die $usage;

#
# env
#
my $sort_path = "$FindBin::Bin/sort";
my $sort_opt = " -m --batch-size=100 --parallel=8"; # "-m" makes the sorting much faster when all input files are sorted (ascending order)

my $sortfile_catOrderings_copyprob    = "copyprobsperlocus.cat.sort";
my $gz_sortfile_catOrderings_copyprob = "copyprobsperlocus.cat.sort.gz";

my $arrayJobID = "";
if ($ENV{LSB_JOBINDEX} ne "") {
  $arrayJobID = $ENV{LSB_JOBINDEX};
} elsif ($ENV{LSB_JOBINDEX} ne "") {
  $arrayJobID = $ENV{LSB_JOBINDEX};
}

#
# vars
#
my $cmd = "";
my $stamp = "";

##############################################################################
# main
#   execute sort -m
#   for alreadly split *.copyprobsperlocus.out_?? files in $dir_each_ordering
##############################################################################

if ($type_painting == 1) {
  #
  # cannot be called from orderedPainting.sh
  #
  my @arr_split_copyprobsperlocus = glob("$dir_each_ordering/*.copyprobsperlocus.out_??");
  if (scalar(@arr_split_copyprobsperlocus) == 0) {
    die "Error: there is no .copyprobsperlocus.out file in $dir_each_ordering ";
  } else {
    foreach my $each_split_copyprobsperlocus (@arr_split_copyprobsperlocus) {
      # get recipient indices from file names, and add them to the 1st column 
      # (otherwise there is no way to know recipient names in the all-against-everyone else painting condition)
      my $i_recipient = $each_split_copyprobsperlocus;
      $i_recipient =~ s/\.copyprobsperlocus\.out_[a-z]+$//g;
      $i_recipient =~ s/^.*_//g;
      $i_recipient =~ s/^0+//g;

      $cmd = " perl -i -pe 's/^([0-9]+) /\$1\t$i_recipient /g' $each_split_copyprobsperlocus ";
      print("$cmd\n");
      if( system("$cmd") != 0) { die("Error: $cmd failed"); };
    }
  }
}

#
# remove unnecessary .out files in each ordering dir
# 
my @arr_dotout_files = glob("$dir_each_ordering/*.out");
foreach my $each_dotout_file (@arr_dotout_files) {
  if ($each_dotout_file !~ /copyprobsperlocus/) {
    unlink($each_dotout_file);
    print "$dir_each_ordering/$each_dotout_file was removed\n";
  }
}

#
# msort (arrayjob from 1 to 9)
#
if ($arrayJobID > 0) {

  my $each_suffix = sprintf("%02d", $arrayJobID);

  $stamp = `date +%Y%m%d_%T`;
  chomp($stamp);
  print "$stamp $each_suffix\n";
  
  $cmd  = "$sort_path $sort_opt ";
  $cmd .= " -T $dir_each_ordering $dir_each_ordering/*.copyprobsperlocus.out_$each_suffix";
  $cmd .= " | gzip > $dir_each_ordering/$gz_sortfile_catOrderings_copyprob" . "." . "$each_suffix"; # split gz files
  print "$cmd\n";
  system("$cmd");

#
# msort (non-arrayjob)
#
} else {

  my $nrecip = `head -2 $phasefile | tail -1`;
  chomp($nrecip);
  if ($type_painting == 2) {
    $nrecip--;
  }
  
  my $nsnp = `head -3 $phasefile | tail -1`;
  chomp($nsnp);

  my $correct_nrow_catCopyprob = $nsnp*$nrecip+($nrecip*2); # ($nrecip*2) is the number of headers

  #
  # get suffix of split files
  #
  my $first_recipient_file = `ls $dir_each_ordering/*.copyprobsperlocus.out_?? | head -1`;
  chomp($first_recipient_file);

  my $first_recipient_prefix = $first_recipient_file;
     $first_recipient_prefix =~ s/\.copyprobsperlocus\.out_..//g;

  my $str_split_suffixes = `ls $first_recipient_prefix*.copyprobsperlocus.out_?? | perl -pe 's/^.*\.copyprobsperlocus\.out_//g'`;
  my @arr_split_suffixes = split(/\n/, $str_split_suffixes);

  while () {

    #
    # if there is an incomplete $sortfile_catOrderings_copyprob files created previously, remove it
    #
    if (  -f "$dir_each_ordering/$gz_sortfile_catOrderings_copyprob") {
      unlink("$dir_each_ordering/$gz_sortfile_catOrderings_copyprob");
    }

    #
    # msort across recipients
    #
    # save them as a concatenated gzip file
    #   (gzipping does not give difference of computational time)
    #
    foreach my $each_suffix (@arr_split_suffixes) {

      $stamp = `date +%Y%m%d_%T`;
      chomp($stamp);
      print "$stamp $each_suffix\n";
      
      $cmd  = "$sort_path $sort_opt ";
      $cmd .= " -T $dir_each_ordering $dir_each_ordering/*.copyprobsperlocus.out_$each_suffix";
      $cmd .= " | gzip -c >> $dir_each_ordering/$gz_sortfile_catOrderings_copyprob";
      print "$cmd\n";
      system("$cmd");

    }

    #
    # check
    #
    $stamp = `date +%Y%m%d_%T`;
    print "$stamp";
   
    print("checking $dir_each_ordering/$gz_sortfile_catOrderings_copyprob ... \n");
    my $nrow = `gzip -dc $dir_each_ordering/$gz_sortfile_catOrderings_copyprob | wc -l`;
    chomp($nrow);
    
    if ($nrow == $correct_nrow_catCopyprob) {
      print("OK, there are $nrow rows in $dir_each_ordering/$gz_sortfile_catOrderings_copyprob\n");
      last;
    } else {
      print "gzip -dc $dir_each_ordering/$gz_sortfile_catOrderings_copyprob | wc -l: $nrow, but must be $correct_nrow_catCopyprob.  Do the msort again.\n";
    }

    print "\n";

  }

}
