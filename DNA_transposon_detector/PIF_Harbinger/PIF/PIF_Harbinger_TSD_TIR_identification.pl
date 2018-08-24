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
    print "Processing sequence: ",$key,"\n";
    @tir_F=();
    @tir_R=();
    my $tag_count=0;
    $key=~/_(\d+)$/;
    my $match_l=$1;

    while ($gSeq{$key}=~/[ATCG]{5}[ATCG]GTTTGG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/[ATCG]{5}T[ATCG]TTTGG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/[ATCG]{5}TG[ATCG]TTGG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/[ATCG]{5}TGT[ATCG]TGG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/[ATCG]{5}TGTT[ATCG]GG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/[ATCG]{5}TGTTT[ATCG]G/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/[ATCG]{5}TGTTTG[ATCG]/g)
    {
        push @tir_F,pos($gSeq{$key});
    }



    while ($gSeq{$key}=~/[ATCG]CAAACA[ATCG]{5}/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/C[ATCG]AAACA[ATCG]{5}/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CC[ATCG]AACA[ATCG]{5}/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CCA[ATCG]ACA[ATCG]{5}/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CCAA[ATCG]CA[ATCG]{5}/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CCAAA[ATCG]A[ATCG]{5}/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CCAAAC[ATCG][ATCG]{5}/g)
    {
        push @tir_R,pos($gSeq{$key});
    }

    
    open (ALIGN,">>tir_tsd.alignment.txt") or die $!;
    open (SEQ,">>PIF.identifed.fasta") or die $!;
    @tir_F=sort {$a<=>$b} @tir_F;
    @tir_R=sort {$a<=>$b} @tir_R;
    @tir_F=uniq(@tir_F);
    @tir_R=uniq(@tir_R);
    #print "Found\n";
    #print join("\t",@tir_F),"\n";
    #print join("\t",@tir_R),"\n";
    foreach my $pos1(@tir_F)
    {
        foreach my $pos2(@tir_R)
        {
            my $pattern=0;
            if($pos1<10000 &&  $pos2>10000+$match_l)
            {
                my $tsd_F=substr($gSeq{$key},$pos1-15,3);
                my $tsd_R=substr($gSeq{$key},$pos2,3);
                my $tir_F=substr($gSeq{$key},$pos1-12,20);
                my $tir_R=substr($gSeq{$key},$pos2-20,20);
                $tir_R=reverse $tir_R;
                $tir_R=~tr/ATCG/TAGC/;
                
                $pattern=1 if $tsd_F eq "TAA" && $tsd_R eq "TAA";
                $pattern=1 if $tsd_F eq "TTA" && $tsd_R eq "TTA";
                open(TIR,">TIR.candidate.txt") or die $!;
                print TIR ">tir_F\n",$tir_F,"\n";
                print TIR ">tir_R\n",$tir_R,"\n";
                my @align1=`mafft --quiet --maxiterate 1000 --clustalout --localpair --thread 1 TIR.candidate.txt`;
                my @match=$align1[5]=~/\*/g;
                    
                if($#match>=12 && $pattern==1)
                {
                    my $seq=substr($gSeq{$key},$pos1-12,$pos2-$pos1+13);
                    print SEQ ">$key","_",$pos1-8,"_",$pos2,"\n";
                    print SEQ $seq,"\n";
                    print ALIGN "TSD: ","\t",$tsd_F,"\t",$tsd_R,"\n";
                    print ALIGN "$key","\t",$pos1,"\t",$pos2,"\n";
                    print ALIGN "$tir_F\t","$tir_R\n";
                    print ALIGN @align1;
                }
            }
        }
    }    
}

