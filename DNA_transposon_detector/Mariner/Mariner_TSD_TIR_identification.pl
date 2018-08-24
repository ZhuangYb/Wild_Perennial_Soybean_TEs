#!/usr/bin/perl -w
use strict;
use Cwd;
use Array::Utils qw(:all);
use List::MoreUtils qw(uniq);

my $genomeseq=shift @ARGV;

my($count,$num_threads,$config,$line,$dir,$DNAref,$Protref,$fasta_file,$blast,$DNA_transposon_detector,$einverted,$checker,$weight);
my($gFasta,$key,$blast_out,$flanking_length,$temp1,$temp2,$TSD1,$TSD2,$job,$tag1,$tag2,$tag3,$tag4,$tag_check);
my(@config,@fasta,@inverted,@genome,@blast,@lmatch,@rmatch,@mule,@mariner,@LTR_repeat,@helitron,@cacta,@hat,@pong,@pif);
my(@tir_F,@tir_R);
my(%gSeq,%temph,%tempn,%ref);

open(GENOME,"$genomeseq") or die $!;
@genome=<GENOME>;

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

LOOP:foreach $key(keys %gSeq)
{
    @tir_F=();
    @tir_R=();
    my $tag_count=0;
    $key=~/_(\d+)$/;
    my $match_l=$1;
    
    while ($gSeq{$key}=~/[ATCG]TCCCTCC/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/C[ATCG]CCCTCC/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CT[ATCG]CCTCC/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CTC[ATCG]CTCC/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CTCC[ATCG]TCC/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CTCCC[ATCG]CC/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CTCCCT[ATCG]C/g)
    {
        push @tir_F,pos($gSeq{$key});
    } 
    while ($gSeq{$key}=~/CTCCCTC[ATCG]/g)
    {
        push @tir_F,pos($gSeq{$key});
    }   


    while ($gSeq{$key}=~/[ATCG]GAGGGAG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/G[ATCG]AGGGAG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/GG[ATCG]GGGAG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/GGA[ATCG]GGAG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/GGAG[ATCG]GAG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/GGAGG[ATCG]AG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/GGAGGG[ATCG]G/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/GGAGGGA[ATCG]/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
   

    open (ALIGN,">>tir_tsd.alignment.txt") or die $!;
    open (SEQ,">>Mariner.identifed.fasta") or die $!;
    @tir_F=sort {$a<=>$b} @tir_F;
    @tir_R=sort {$a<=>$b} @tir_R;
    @tir_F=uniq(@tir_F);
    @tir_R=uniq(@tir_R);
    #print "Found\n";
    print join("\t",@tir_F),"\n";
    print join("\t",@tir_R),"\n";
    foreach my $pos1(@tir_F)
    {
        foreach my $pos2(@tir_R)
        {
            if($pos1<10000 &&  $pos2>10000+$match_l)
            {
                my $tsd_F=substr($gSeq{$key},$pos1-13,5);
                my $tsd_R=substr($gSeq{$key},$pos2,5);
                my $tir_F=substr($gSeq{$key},$pos1-8,20);
                my $tir_R=substr($gSeq{$key},$pos2-18,20);
                $tir_R=reverse $tir_R;
                $tir_R=~tr/ATCG/TAGC/;
                open(TSD,">TSD.candidate.txt") or die $!;
                print TSD ">tsd_F\n",$tsd_F,"\n";
                print TSD ">tsd_R\n",$tsd_R,"\n";
                close TSD;
                
                open(TIR,">TIR.candidate.txt") or die $!;
                print TIR ">tir_F\n",$tir_F,"\n";
                print TIR ">tir_R\n",$tir_R,"\n";
                my @align1=`mafft --quiet --maxiterate 1000 --clustalout --localpair --thread 1 TIR.candidate.txt`;
                my @match=$align1[5]=~/\*/g;
                if($tsd_F=~/TA/ && $tsd_R=~/TA/ && $#match>=12)
                {
                    my $seq=substr($gSeq{$key},$pos1-8,$pos2-$pos1+9);
                    print SEQ ">$key","_",$pos1-8,"_",$pos2,"\n";
                    print SEQ $seq,"\n";
                    print ALIGN "$key","\t",$pos1,"\t",$pos2,"\n";
                    print ALIGN "$tir_F\t","$tir_R\n";
                    print ALIGN @align1;
                }
                elsif($tsd_F=~/AT/ && $tsd_R=~/AT/ && $#match>=12)
                {
                    my $seq=substr($gSeq{$key},$pos1-8,$pos2-$pos1+9);
                    print SEQ ">$key","_",$pos1-8,"_",$pos2,"\n";
                    print SEQ $seq,"\n";
                    print ALIGN "$key","\t",$pos1,"\t",$pos2,"\n";
                    print ALIGN "$tir_F\t","$tir_R\n";
                    print ALIGN @align1;
                }
            }
        }
    }

}

