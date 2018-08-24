#!/usr/bin/perl -w
# perl TSD.extract.pl /home/zyb/Database/Fasta/Gmax.v1.single.fasta MULE.cor.txt
use strict;
use Cwd;
use Array::Utils qw(:all);

my($count,$num_threads,$config,$line,$dir,$DNAref,$Protref,$fasta_file,$blast,$DNA_transposon_detector,$einverted);
my($gFasta,$blast_out,$flanking_length,$temp1,$temp2,$TSD1,$TSD2,$job);
my(@config,@fasta,@inverted,@genome,@blast,@lmatch,@rmatch,@mule,@mariner,@LTR_repeat,@helitron,@cacta,@hat,@pong,@pif);
my(@temp1,@temp2,@temp3);
my(%gSeq,%temph,%tempn);

my $file=shift @ARGV;
my $file2=shift @ARGV;
open(F,"$file") or die $!;
open(F2,"$file2") or die $!;
@genome=<F>;
my @list=<F2>;
for(my $count=0;$count<$#genome;$count++)
{
        chomp $genome[$count];
        if($genome[$count]=~/>(.+)/)
        {
                my $key=$1;
                chomp $genome[$count+1];
                $gSeq{$key}=$genome[$count+1];
        }
}

for(my $i=0;$i<$#list;$i++)
{
        if($list[$i]=~/^(.+)\:(\d+)\.\.(\d+)/)
        {
                my $id=$1;
                my $star_line=$2;
                my $length_diff=$3-$2;
                my $end_line=$3;
                print $id,"\t".$star_line,"\t",$end_line,"\n";
                my $ltsd=substr($gSeq{$id},$star_line-12,12);
                my $rtsd=substr($gSeq{$id},$end_line-1,12);
                open(OUT3,">TSD.fasta") or die $!;
                print OUT3 ">left_TSD\n",$ltsd,"\n";
                print OUT3 ">right_TSD\n",$rtsd,"\n";
                print $ltsd,"\t",$rtsd,"\n";
                my @align=`mafft --quiet --clustalout --auto TSD.fasta`;
                print @align[$#align-3..$#align];

                my $a=substr($gSeq{$id},$star_line,120);
                my $b=substr($gSeq{$id},$end_line-120,120);
                $b=reverse $b;
                $b=~tr/ATCG/TAGC/;
                open(OUT4,">TSD.fasta") or die $!;
                print OUT4 ">left_TSD\n",$a,"\n";
                print OUT4 ">right_TSD\n",$b,"\n";
                @align=`mafft --quiet --clustalout --maxiterate 300 --localpair TSD.fasta`;
                print @align,"\n\n\n";

                open(OUT5,">>Identified.TIR.fasta") or die $!;
                print OUT5 ">",$id,"_",$star_line,"_",$end_line,"_F","\n",$a,"\n";
                print OUT5 ">",$id,"_",$star_line,"_",$end_line,"_R","\n",$b,"\n";

        }



}