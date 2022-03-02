#!/usr/bin/perl
use strict;
use Getopt::Long;

################Written By Anuj Srivastava#######################################

my $usage = <<'USAGE';

####################Cosmic_Annotation.pl#############################

        usage: Cosmic_Annotation.pl [options]

                -i1= infile_1.fastq  [Cosmis Coding and Non-Coding Annotation file]
                -i2= infile_2.fastq  [Snpeff Annotated SNP and Indel files]

##################################################################################


USAGE

my ( $infile1 , $infile2);

my $result = GetOptions(
    "i1=s"  => \$infile1,
    "i2=s"  => \$infile2,
)
;


die $usage unless ($infile1); ##### Mandatory arguments
die $usage unless ($infile2); ##### Mandatory arguments


open(FILE_COSMIC,  $infile1); 
open(FILE_Variant, $infile2);


my $readInfile; 
my $readCosmicfile;
my $finalkey;  my $NameChr; my $Cor; my $CID;  my %Chash = ();

my $genomeMatch = 0;

Label1:while(my $readCosmicfile = <FILE_COSMIC>)
{
  if($readCosmicfile =~ /GRCh38/)
   {
       $genomeMatch = 1;
       last Label1;
 
   }

}

if($genomeMatch !=  1)
{
  die"genome does not match hg38\n";

}

close(FILE_COSMIC);


open(FILE_COSMIC,  $infile1);

while(my $readCosmicfile = <FILE_COSMIC>)
{
 if($readCosmicfile !~ /^#.*$/)
  {
      if($readCosmicfile =~ /^\s*(.*?)\s+(.*?)\s+(.*?)\s+.*$/)
      {
          $NameChr  = $1;
          $Cor      = $2;
          $CID      = $3;

          if($NameChr =~ /MT/)
           {
               $NameChr = 'M';
           }
          $finalkey = "chr"."$NameChr"."_"."$Cor";


          if(!exists($Chash{$finalkey}))
          {
            $Chash{$finalkey} = $CID;
            next;
          }
          else
          {
            $Chash{$finalkey} = "$Chash{$finalkey}".";"."$CID";

          }

          

      }


  }
} 

my $InputKey;
my $InputID;
my $chr; my $pos; my $etc; 

while(my $readInfile = <FILE_Variant>)
{
 if($readInfile =~ /^#.*$/)
  {
     print "$readInfile";
 
  }

 elsif($readInfile =~ /^\s*(.*?)\s+(.*?)\s+(.*?)\s+(.*)$/)
  {
     $InputKey =  "$1"."_"."$2";
     $InputID = $3;
     $chr = $1;
     $pos = $2; 
     $etc = $4;
     
     if(exists($Chash{$InputKey}))
      {
           if($InputID =~ /\./)  
           { 
             print "$chr\t$pos\t$Chash{$InputKey}\t$etc\n";
           }
           else
           {
             print "$chr\t$pos\t$InputID;$Chash{$InputKey}\t$etc\n";
              
           }


      }
     else
      {
        print "$readInfile";

      }

          
  }
}
