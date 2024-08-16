#!/usr/bin/perl -w
use POSIX;
use File::Basename;

# This script annotates ensembl genes with copy number and breakpoints
# perl ensemblegenes_cnv_break.pl *.segments_raw.extend.txt mart_export_gene_chr1-Y.hg19ensembl75-85.08232016.txt

if ($#ARGV != 1) {
	print "This scripts requires: <file_cn> <file_gene> \n";
	exit(-1);
}

$file_cn = $ARGV[0];		
$file_gene = $ARGV[1];

$file_output = basename($file_cn,".txt").".ensgene_cnvbreak.txt";
open(OUTFILE, ">$file_output");

open(GENEFILE, "$file_gene") or die "can't open $file_gene: $!";
$gene = <GENEFILE>;
chomp($gene);

open(CNFILE, "$file_cn") or die "can't open $file_cn: $!";
@data = <CNFILE>;
close(CNFILE);
chomp(@data);

#print OUTFILE "$tmp\tstartext\tendext\tstartext_desc\tendext_desc\tCN_raw\tLOH\tparm_fraction\tqarm_fraction\tploidy\tcopydiff_2\tcopydiff_ploidy\tlogratio_2\tlogratio_ploidy\n";
print OUTFILE "$gene\tnum_cnv_seg\tseg_desc\tploidy\tnMajor\tnMinor\tnAraw\tnBraw\tCN_raw\tLOH\tcopydiff_2\tcopydiff_ploidy\tlogratio_2\tlogratio_ploidy\tnMajor_max\tnMinor_max\tnAraw_max\tnBraw_max\tCN_raw_max\tLOH_max\tcopydiff_2_max\tcopydiff_ploidy_max\tlogratio_2_max\tlogratio_ploidy_max\n";

while ($gene = <GENEFILE>) {
    
    chomp($gene);
    @line = split(/\t/, $gene);
    $chr = $line[2];
    $start = $line[3];
    $end = $line[4];
    
    #$cnraw1=999;
    $numseg=0;
    $region="";
    %segline = ();
    @n = ();
    
    for ($j=1; $j<=$#data; $j++) {
        @segment = split(/\t/, $data[$j]);
        
        $chr_cn = $segment[1];
        $pos1 = $segment[2];
        $pos2 = $segment[3];
        $pos1ext = $segment[9];
        $pos2ext = $segment[10];
        $left = $segment[11];
        $right = $segment[12];
        $cnraw = $segment[13];
        
        if (($chr_cn eq $chr) && ($start <= $pos2ext) && ($end >= $pos1ext)) { #overlap
            #$numseg++;
            push(@n, $cnraw);
            $segline{$cnraw} = [ @segment ];
            
            #check if overlap with regions with no call
            if (($start <= $pos1) && ($end >= $pos1ext)) {
                $region = $region.$left.";";
            }
            if (($start <= $pos2ext) && ($end >= $pos2)) {
                $region = $region.$right.";";
            }
            
            #if ($cnraw < $cnraw1) {
            #    $cnraw1 = $cnraw;
            #    $count = $j;
            #}
        }
    }
    
    if ($region eq "") {
        $region = "NA";
    }
    
    if ($#n >= 0) {
        
        $numseg = $#n +1;
        @sortn = sort{ $a <=> $b } @n;
    
        $nA = $segline{$sortn[0]}[4];
        $nB = $segline{$sortn[0]}[5];
        $rawA = $segline{$sortn[0]}[6];
        $rawB = $segline{$sortn[0]}[7];
        $cnraw = $segline{$sortn[0]}[13];
        $loh = $segline{$sortn[0]}[14];
        $ploidy= $segline{$sortn[0]}[17];
        $copydiff1 = $segline{$sortn[0]}[18];
        $copydiff2 = $segline{$sortn[0]}[19];
        $logratio1 = $segline{$sortn[0]}[20];
        $logratio2 = $segline{$sortn[0]}[21];
    
        $outline = "$gene\t$numseg\t$region\t$ploidy\t$nA\t$nB\t$rawA\t$rawB\t$cnraw\t$loh\t$copydiff1\t$copydiff2\t$logratio1\t$logratio2\t";
    
        if ($numseg > 1 ) {
            $nA = $segline{$sortn[$#sortn]}[4];
            $nB = $segline{$sortn[$#sortn]}[5];
            $rawA = $segline{$sortn[$#sortn]}[6];
            $rawB = $segline{$sortn[$#sortn]}[7];
            $cnraw = $segline{$sortn[$#sortn]}[13];
            $loh = $segline{$sortn[$#sortn]}[14];
            $copydiff1 = $segline{$sortn[$#sortn]}[18];
            $copydiff2 = $segline{$sortn[$#sortn]}[19];
            $logratio1 = $segline{$sortn[$#sortn]}[20];
            $logratio2 = $segline{$sortn[$#sortn]}[21];
        }
        else {
            $nA = "NA";
            $nB = "NA";
            $rawA = "NA";
            $rawB = "NA";
            $cnraw = "NA";
            $loh = "NA";
            $copydiff1 = "NA";
            $copydiff2 = "NA";
            $logratio1 = "NA";
            $logratio2 = "NA";

        }
    
        $outline = $outline."$nA\t$nB\t$rawA\t$rawB\t$cnraw\t$loh\t$copydiff1\t$copydiff2\t$logratio1\t$logratio2";
        print OUTFILE "$outline\n";
    }
}

close (GENEFILE);
close (OUTFILE);
