#!/usr/bin/perl -w
use POSIX;
use File::Basename;

# This script adds to segment file the arm fraction, LOH and CN diff and log ratio relative to 2 and ploidy
# The segments are extended

# perl segment_raw_annotate.pl *segments_raw.txt *ploidy.txt hg38_chromosome_arm.txt [male, female, unknown]

if ($#ARGV != 3) {
	print "This scripts requires: <file_cn> <file_ploidy> <file_arm> <gender> \n";
	exit(-1);
}

$file_cn = $ARGV[0];		
$file_ploidy = $ARGV[1];
$file_arm = $ARGV[2];
$gender = $ARGV[3];

$file_output = basename($file_cn,".txt").".extend.txt";

$ploidy = `cat $file_ploidy`;
chomp($ploidy);

# $gender = `cat $file_gender`;
# chomp($gender);

if (($gender eq "female") || ($gender eq "unknown")) {
    $cn_factor = 1;
}
elsif ($gender eq "male") {
    $cn_factor= 0.5;
}

$tmp = `cat $file_arm | awk 'NR>1'`;
@arm = split(/\n/,$tmp);
chomp(@arm);

open(CN, "$file_cn") or die "can't open $file_cn: $!";
$tmp = <CN>;
chomp($tmp);

open(OUTFILE, ">$file_output");
print OUTFILE "$tmp\tstartext\tendext\tstartext_desc\tendext_desc\tCN_raw\tLOH\tparm_fraction\tqarm_fraction\tploidy\tcopydiff_2\tcopydiff_ploidy\tlogratio_2\tlogratio_ploidy\n";

open(TMPFILE, ">tmp.txt");

#merge segments
$tmp = <CN>;
chomp($tmp);
@line = split(/\t/,$tmp);
print "@line\n";
$sample = $line[0];
$chromo = $line[1];
$n1 = $line[4];
$n2 = $line[5];
$cn1 = $line[6];
$cn2 = $line[7];
$start = $line[2];
$end = $line[3];
$num = $line[8];

print "$num\n";

while ($tmp = <CN>) {
    chomp($tmp);
    @line = split(/\t/,$tmp);
    
    if (($chromo eq $line[1]) && ($cn1 == $line[6]) && ($cn2 == $line[7])) {
        $end = $line[3];
        $num = $num + $line[8];
    }
    else {
        print TMPFILE "$sample\t$chromo\t$start\t$end\t$n1\t$n2\t$cn1\t$cn2\t$num\n";
        $sample = $line[0];
        $chromo = $line[1];
        $n1 = $line[4];
        $n2 = $line[5];
        $cn1 = $line[6];
        $cn2 = $line[7];
        $start = $line[2];
        $end = $line[3];
        $num = $line[8];
    }
}
#lastline
print TMPFILE "$sample\t$chromo\t$start\t$end\t$n1\t$n2\t$cn1\t$cn2\t$num\n";

close (CN);
close (TMPFILE);

open(CN, "tmp.txt") or die "can't open tmp.txt: $!";
@seg = <CN>;
chomp(@seg);
close (CN);
$n = 0;

for ($j=0; $j<$#seg; $j++) {
    
    @array1 = split(/\t/,$seg[$j]);
    @array2 = split(/\t/,$seg[$j+1]);
    #$x1 = $array1[2];
    $x2 = $array1[3];
    $y1 = $array2[2];
    #$y2 = $array2[3];
    
    if ($array1[1] ne $n) { #first line for chr
        
        $n = $array1[1];
        $left = 0;
        $left1 = "telomere";
        
        for ($i=1; $i<=$#arm; $i+=2) {
            @line = split(/\t/,$arm[$i]);
            if ($n eq substr($line[0],3)) {
                $a = $line[1];
                $b = $line[2];
            }
        }
        
        if ($array2[1] ne $n) { #last line for chr
            $right = $b;
            $right1 = "telomere";
        }
        elsif (($x2 < $a) && ($y1 > $a)) {
            $right = $a;
            $right1 = "centromere";
        }
        else {
            $right = floor(($x2 + $y1)/2);
            $right1 = "no_probe";
        }
    }
    else {
        
        $left = $right + 1;
        $left1 = $right1;
        
        if ($array2[1] ne $n) { #last line for chr
            
            $right = $b;
            $right1 = "telomere";
        }
        elsif (($x2 < $a) && ($y1 > $a)) {
            $right = $a;
            $right1 = "centromere";
        }
        else {
            $right = floor(($x2 + $y1)/2);
            $right1 = "no_probe";
        }
    }
    
    $copy = $array1[6] + $array1[7];
    if ($array1[6] >= 0.5 && $array1[7] <= 0.1) {
        $loh=1;
    }
    else {
        $loh=0;
    }
    
    for ($i=0; $i<=$#arm; $i+=2) {
        @line = split(/\t/,$arm[$i]);
        if ($n eq substr($line[0],3)) {
            if (($right>=$line[1]) && ($left<=$line[2])) {
                @tmp = ($left,$right,$line[1],$line[2]);
                @sorttmp = sort{ $a <=> $b } @tmp;
                $overlap1=($sorttmp[2]-$sorttmp[1])/($line[2]-$line[1]);
            }
            else {
                $overlap1=0;
            }
        }
    }
    
    for ($i=1; $i<=$#arm; $i+=2) {
        @line = split(/\t/,$arm[$i]);
        if ($n eq substr($line[0],3)) {
            if (($right>=$line[1]) && ($left<=$line[2])) {
                @tmp = ($left,$right,$line[1],$line[2]);
                @sorttmp = sort{ $a <=> $b } @tmp;
                $overlap2=($sorttmp[2]-$sorttmp[1])/($line[2]-$line[1]);
            }
            else {
                $overlap2=0;
            }
        }
    }

    if (($n eq "X") || ($n eq "Y")) {
        $diff1=$copy - ($cn_factor * 2);
        $diff2=$copy- ($cn_factor * $ploidy);
        $logratio1 = log(($copy+0.01)/($cn_factor * 2))/log(2);
        $logratio2 = log(($copy+0.01)/($cn_factor * $ploidy))/log(2);
    }
    else {
        $diff1=$copy-2;
        $diff2=$copy-$ploidy;
        $logratio1 = log(($copy+0.01)/2)/log(2);
        $logratio2 = log(($copy+0.01)/$ploidy)/log(2);
    }
    
    print OUTFILE "$seg[$j]\t$left\t$right\t$left1\t$right1\t$copy\t$loh\t$overlap1\t$overlap2\t$ploidy\t$diff1\t$diff2\t$logratio1\t$logratio2\n";
}

@array1 = split(/\t/,$seg[$#seg]);

if ($array1[1] ne $n) { #first line for chr
    
    $n = $array1[1];
    $left = 0;
    $left1 = "telomere";
    
    for ($i=1; $i<=$#arm; $i+=2) {
        @line = split(/\t/,$arm[$i]);
        if ($n eq substr($line[0],3)) {
            $a = $line[1];
            $b = $line[2];
        }
    }
    
    $right = $b;
    $right1 = "telomere";
    
}
else {
    
    $left = $right + 1;
    $left1 = $right1;
    
    $right = $b;
    $right1 = "telomere";
    
}

$copy = $array1[6] + $array1[7];
if ($array1[6] >= 0.5 && $array1[7] <= 0.1) {
    $loh=1;
}
else {
    $loh=0;
}

for ($i=0; $i<=$#arm; $i+=2) {
    @line = split(/\t/,$arm[$i]);
    if ($n eq substr($line[0],3)) {
        if (($right>=$line[1]) && ($left<=$line[2])) {
            @tmp = ($left,$right,$line[1],$line[2]);
            @sorttmp = sort{ $a <=> $b } @tmp;
            $overlap1=($sorttmp[2]-$sorttmp[1])/($line[2]-$line[1]);
        }
        else {
            $overlap1=0;
        }
    }
}

for ($i=1; $i<=$#arm; $i+=2) {
    @line = split(/\t/,$arm[$i]);
    if ($n eq substr($line[0],3)) {
        if (($right>=$line[1]) && ($left<=$line[2])) {
            @tmp = ($left,$right,$line[1],$line[2]);
            @sorttmp = sort{ $a <=> $b } @tmp;
            $overlap2=($sorttmp[2]-$sorttmp[1])/($line[2]-$line[1]);
        }
        else {
            $overlap2=0;
        }
    }
}

if (($n eq "X") || ($n eq "Y")) {
    $diff1=$copy - ($cn_factor * 2);
    $diff2=$copy- ($cn_factor * $ploidy);
    $logratio1 = log(($copy+0.01)/($cn_factor * 2))/log(2);
    $logratio2 = log(($copy+0.01)/($cn_factor * $ploidy))/log(2);
}
else {
    $diff1=$copy-2;
    $diff2=$copy-$ploidy;
    $logratio1 = log(($copy+0.01)/2)/log(2);
    $logratio2 = log(($copy+0.01)/$ploidy)/log(2);
}

print OUTFILE "$seg[$j]\t$left\t$right\t$left1\t$right1\t$copy\t$loh\t$overlap1\t$overlap2\t$ploidy\t$diff1\t$diff2\t$logratio1\t$logratio2\n";

close(CN);
close (OUTFILE);
