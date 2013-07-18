#! /usr/bin/perl

use Data::Dumper;
use Getopt::Std;
#use Tie::IxHash;
use strict;
use warnings;

use vars qw ($opt_f $opt_v);
&getopts('f:v');

my $usage = <<_EOH_;

# -f file.hap.txt

# [-v]

_EOH_
;

#
# IN
#

my $inFile   = $opt_f or die $usage;
open(IN, $inFile);
while (my $line = <IN>) {
  chomp $line;
  print length($line) . "\n";
}
close(IN);
