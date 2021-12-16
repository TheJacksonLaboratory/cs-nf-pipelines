#!/usr/bin/perl
use strict;
use warnings;

open(FILEIN1, $ARGV[0]) || die "cannot open the file"; ####filter *_stat
open(FILEIN2, $ARGV[1]) || die "cannot open the file"; ####rsem aln stat
open(FILEIN3, $ARGV[2]) || die "cannot open the file"; ####picard stats

my ($value1, $value2, $value3);


while(my $readFile1 = <FILEIN1>)
{
  if($readFile1 =~ /^\s*Percentage of HQ reads\s+(.*?)\s+(.*)/)
   {
           $value1 = $1;
   }
  elsif($readFile1 =~ /^\s*Total number of reads\s+(.*?)\s+(.*).*$/)
   {
           $value2 = $1;

   }
  elsif($readFile1 =~ /^\s*Total number of HQ filtered reads\s+(.*?)\s+(.*).*$/)
   {
           $value3 = $1;
   }

}

print "Total number of Reads\t$value2\n";
print "Total number of HQ filtered reads\t$value3\n";
print "Percentage of HQ reads\t$value1\n";

my $flag = 0;



while(my $readFile2 = <FILEIN2>)
{
   if(($readFile2 =~ /^\s*(\d+)\s+.*$/) && ($flag == 0))
    {
        print "Total number of input reads for RSEM transcriptome Alignment\t$1\n";
        $flag = 1;
    }
   elsif(($readFile2 =~ /^\s*(.*?)\s+\((.*)\).*$/) && ($flag == 1))
    {
        print "Total number of un-paired reads for RSEM transcriptome Alignment\t$1\t$2\n";
        $flag = 2;
    }
   elsif(($readFile2 =~ /^\s*(.*?)\s+\((.*)\).*$/) && ($flag == 2))
    {
        print "Total number of reads aligned 0 times from RSEM transcriptome Alignment\t$1\t$2\n";
        $flag = 3;
    }
   elsif(($readFile2 =~ /^\s*(.*?)\s+\((.*)\).*$/) && ($flag == 3))
    {
        print "Total number of reads aligned exactly 1 time from RSEM transcriptome Alignment\t$1\t$2\n";
        $flag = 4;
    }
   elsif(($readFile2 =~ /^\s*(.*?)\s+\((.*)\).*$/) && ($flag == 4))
    {
        print "Total number of reads aligned  >1 time from RSEM transcriptome Alignment\t$1\t$2\n";
        $flag = 5;
    }
   elsif(($readFile2 =~ /^\s*(.*?)\%\s+.*$/) && ($flag == 5))
    {
        print "Overall alignment rate\t$1\%\n";
    }


}

$flag = 0;

my @splitHeader = ();
my @splitValue  = ();

while(my $readFile3 = <FILEIN3>)
{

   if(($readFile3 =~ /^\s*##\s+METRICS\s+CLASS\s+picard.analysis.RnaSeqMetrics.*$/) && ($flag == 0))
    {
           $flag = 1;
           next;
    }
   elsif(($readFile3 =~ /^\s*PF_BASES.*$/) && ($flag == 1))
    {
           $flag = 2;
           chomp $readFile3;
           @splitHeader = split("\t", $readFile3);
           next;
    }
   elsif($flag == 2)
    {
           $flag = 3;
           chomp $readFile3;
           @splitValue = split("\t", $readFile3);
           next;
    }

}

for (my $i =0; $i<=$#splitHeader-3; $i++)
{
   print "$splitHeader[$i]\t$splitValue[$i]\n";

}
