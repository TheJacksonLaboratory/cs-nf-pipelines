#!/usr/bin/perl -w
use POSIX;
use File::Basename;

# This script annotates ensembl transcripts with copy number and breakpoints
# perl ../ensembltranscript_cnv_break_sequenza.pl T100_1.fastq.gz_segments.txt ../ensembl93_transcript_081618.txt 

if ($#ARGV != 1) {
	print "This scripts requires: <file_cn> <file_gene> \n";
	exit(-1);
}

## The script was adjusted to print each overlapping CNV segment to a new line. Originally written to assume 0, 1, or 2 overlaps, the overlap printout was generalized to 'N' cases. 

$file_cn = $ARGV[0];		
$file_gene = $ARGV[1];

$file_output = basename($file_cn,".txt").".enstranscript_cnvbreak.txt";
open(OUTFILE, ">$file_output");

$file_output1 = basename($file_cn,".txt").".ensgene_cnvbreak.txt";
open(OUTFILE1, ">$file_output1");

open(GENEFILE, "$file_gene") or die "can't open $file_gene: $!";
$gene = <GENEFILE>;
chomp($gene);

open(CNFILE, "$file_cn") or die "can't open $file_cn: $!";
@data = <CNFILE>;
close(CNFILE);
chomp(@data);

print OUTFILE "cnv_seg_number\t$gene\tBf\tdepth.ratio\tCNt\tA\tB\n";

print OUTFILE1 "cnv_seg_number\t$gene\tBf\tdepth.ratio\tCNt\tA\tB\n";

$geneid="";

while ($gene = <GENEFILE>) {
    
    chomp($gene);
    @line = split(/\t/, $gene);
    
    $chr = $line[3];
    $start = $line[6];
    $end = $line[7];
    
    $numseg=0;
    %segline = ();
    @n = ();
    
    for ($j=1; $j<=$#data; $j++) {
        
        $data[$j] =~ s/"//g;
        @segment = split(/\t/, $data[$j]);
        
        $chr_cn = substr($segment[0],3);
        $pos1 = $segment[1];
        $pos2 = $segment[2];
        $DR = $segment[6];

        if (($chr_cn eq $chr) && ($start <= $pos2) && ($end >= $pos1)) { #overlap
            $seg_count++;
            push(@n, $DR);
            $segline{$DR} = [ @segment ];

        }
    }

    if ($#n >= 0) {
        $numseg = $#n +1;
        
        @sortn = sort{ $a <=> $b } @n;
    
        for($loop_index = 0; $loop_index <= $#sortn; $loop_index++) {
            $Bf = $segline{$sortn[$loop_index]}[3];
            $DR = $segline{$sortn[$loop_index]}[6];
            $CNt = $segline{$sortn[$loop_index]}[9];
            $A = $segline{$sortn[$loop_index]}[10];
            $B = $segline{$sortn[$loop_index]}[11];

            $cnv_count = $loop_index + 1;

            $outline = "$cnv_count\t$gene\t$Bf\t$DR\t$CNt\t$A\t$B";
            print OUTFILE "$outline\n";
        }
    }
    
    else {    
        print OUTFILE "0\t$gene\tNA\tNA\tNA\tNA\tNA\n";
    }
    
    $seg_count = 0;
    
    if ($line[0] ne $geneid) {
        
        $geneid=$line[0];
        $start = $line[4];
        $end = $line[5];
        
        #$cnraw1=999;
        $numseg=0;
        #$region="";
        %segline = ();
        @n = ();
        
        for ($j=1; $j<=$#data; $j++) {
            
            $data[$j] =~ s/"//g;
            @segment = split(/\t/, $data[$j]);
            
            $chr_cn = substr($segment[0],3);
            $pos1 = $segment[1];
            $pos2 = $segment[2];
            $DR = $segment[6];
            
            if (($chr_cn eq $chr) && ($start <= $pos2) && ($end >= $pos1)) { #overlap
                push(@n, $DR);
                $segline{$DR} = [ @segment ];
                
            }
        }
        
        if ($#n >= 0) {
            
            $numseg = $#n +1;
            @sortn = sort{ $a <=> $b } @n;
                
            for($loop_index = 0; $loop_index <= $#sortn; $loop_index++) {
                $Bf = $segline{$sortn[$loop_index]}[3];
                $DR = $segline{$sortn[$loop_index]}[6];
                $CNt = $segline{$sortn[$loop_index]}[9];
                $A = $segline{$sortn[$loop_index]}[10];
                $B = $segline{$sortn[$loop_index]}[11];

                $cnv_count = $loop_index + 1;

                $outline = "$cnv_count\t$gene\t$Bf\t$DR\t$CNt\t$A\t$B";
                print OUTFILE1 "$outline\n";
            }
        }
        
        else {
            print OUTFILE1 "0\t$gene\tNA\tNA\tNA\tNA\tNA\n";
        }

        
    }

    
}

close (GENEFILE);
close (OUTFILE);
close (OUTFILE1);

`cat $file_output1 | cut -f1,2,4-7,10,11,13,14,17,18,20,22,24-> tmp.txt`;
`mv tmp.txt $file_output1`;
## remove transcript information from gene file. 
