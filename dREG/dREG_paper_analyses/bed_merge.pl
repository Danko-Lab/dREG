#!/usr/bin/perl
#	Program outputs the column specificed from a TSV file (0 = first column ...)

# Perl trim function to remove whitespace from the start and end of the string

my $p_chrom;
my $p_chromStart;
my $p_chromEnd;

$join= int($ARGV[0]);
$doprint=0;

while(<STDIN>) {
  @SPL = split(/[\s\t]/);
  $chrom= $SPL[0];
  $chromStart=$SPL[1];
  $chromEnd=$SPL[2];

  if(($p_chrom eq $chrom) && (($p_chromEnd+$join)>=$chromStart)) {
    $p_chromEnd=$chromEnd;
  } else {
    if($doprint == 1) {
      print $p_chrom."\t".$p_chromStart."\t".$p_chromEnd."\n";
    }
    $p_chrom= $chrom;
    $p_chromStart= $chromStart;
    $p_chromEnd= $chromEnd;
  }
  $doprint = 1;
}
## Print the last entry.
print $p_chrom."\t".$p_chromStart."\t".$p_chromEnd."\n";

