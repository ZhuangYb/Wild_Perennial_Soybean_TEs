#!/usr/bin/perl -w
use strict;
use Cwd;
use Array::Utils qw(:all);
use List::MoreUtils qw(uniq);

my $genomeseq=shift @ARGV;

my($count,$num_threads,$config,$line,$dir,$DNAref,$Protref,$fasta_file,$blast,$DNA_transposon_detector,$einverted,$checker,$weight);
my($gFasta,$key,$blast_out,$flanking_length,$temp1,$temp2,$TSD1,$TSD2,$job,$tag1,$tag2,$tag3,$tag4,$tag_check);
my(@config,@fasta,@inverted,@genome,@blast,@lmatch,@rmatch,@mule,@mariner,@LTR_repeat,@helitron,@cacta,@hat,@pong,@pif);
my(@temp1,@temp2,@temp3);
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
    my $tag_count=0;
    $key=~/_(\d+)$/;
    my $match_l=$1;
    my $t1;
    my $t2;
    my $t; #internal seq
    my $lmatch;
    my $rmatch;
    my $flength=length($gSeq{$key});
    open(OUT1,">temp1.fasta") or die $!;
    open(OUT2,">temp2.fasta") or die $!;

    $t1=substr($gSeq{$key},0,10000);
    $t=substr($gSeq{$key},10000,$match_l);
    $t2=substr($gSeq{$key},10000+$match_l,$flength-10000-$match_l);
    print OUT1 ">Left\n",$t1,"\n";
    print OUT2 ">right\n",$t2,"\n";

    open(TIRALIGN,">>TIR.identified.alignment.txt") or die $!;
    `jellyfish count -m 9 -s 100000 temp1.fasta  -o left.jf`;
    `jellyfish count -m 9 -s 100000 temp2.fasta  -o right.jf`;
    `jellyfish dump -c left.jf >left.kmer`;
    `jellyfish dump -c right.jf >right.kmer`;
        
    open(F3,"left.kmer") or die $!;
    open(F4,"right.kmer") or die $!;
    my @left=<F3>;
    my @right=<F4>;
    my @left_final=();
    my @right_final=();

    foreach my $kleft(@left)
    {
        $kleft=~/(.+?)\s+/;
        push @left_final, $1;
    }
    foreach my $kright(@right)
    {
        $kright=~/(.+?)\s+/;
        push @right_final, $1;
    }

    my @shared=intersect(@left_final,@right_final);
    unless ($#shared==-1)
    {    
        foreach my $shared_kmer(@shared)
        {
            @lmatch=();
            @rmatch=();
            while($t1=~/$shared_kmer/g)
            {
                push @lmatch,pos($t1);
            }
            
            while($t2=~/$shared_kmer/g)
            {
                push @rmatch,pos($t2);
            }
            
            foreach my $lindex(@lmatch)
            {
                my $lseq=substr($t1,$lindex,120);
                my $lflank=substr($t1,$lindex,2);
                my $tsd_l=substr($t1,$lindex-9,9);
                my $lflank2=substr($t1,$lindex-11,2);
                
                foreach my $rindex(@rmatch)
                {   
                    my @head_star=();
                    open(OUT3,">TSD1.fasta") or die $!;
                    my $rflank=substr($t2,$rindex,2);
                    my $tsd_r=substr($t2,$rindex-9,9);
                    my $rflank2=substr($t2,$rindex-11,2);
                    my $rseq=substr($t2,$rindex-129,120);
                    my $rseq_ori=$rseq;
                    unless($lflank eq $rflank || $lflank2 eq $rflank2)
                    {
                        my $tsd_2l=substr($t1,$lindex-10,11);
                        my $tsd_2r=substr($t2,$rindex-10,11);
                        open(OUTW,">test.fasta");
                        print OUTW ">left_TSD\n",$tsd_2l,"\n";
                        print OUTW ">right_TSD\n",$tsd_2r,"\n";
                        my @alignw=`mafft --quiet --clustalout --auto --thread 1 test.fasta`;

                        $rseq=reverse $rseq;
                        $rseq=~tr/ATCG/TAGC/;
                        print OUT3 ">left_TSD\n",$lseq,"\n";
                        print OUT3 ">right_TSD\n",$rseq,"\n";
                                                        
                        my @align=`mafft --quiet --clustalout --auto --thread 1 TSD1.fasta`;
                        my $total_star=0;
                        @head_star=substr($align[5],16,30)=~/\*/g;
                        foreach my $star_line(@align)
                        {
                            my @stare=$star_line=~/\*/g if $star_line=~/\*/;
                            $total_star=$total_star+$#stare+1 if $star_line=~/\*/;
                        }
                        if($total_star>=70)
                        {   

                            $tag3=substr($lseq,0,60);
                            $tag4=substr($rseq,0,60);
       
                            my $timer=0;
                            my @timer=();
                            
                            @timer=$tag3=~/A/g;
                            $timer++ if $#timer/60 >=0.4;
                            @timer=$tag3=~/T/g;
                            $timer++ if $#timer/60 >=0.4; 
                            @timer=$tag3=~/C/g;
                            $timer++ if $#timer/60 >=0.4;
                            @timer=$tag3=~/G/g;
                            $timer++ if $#timer/60 >=0.4;

                            my @seq1=();
                            my @seq2=();
                            
                            
                            if($timer<=1 && $#head_star>=21)
                            {
                                open(TAGOUT,">>120bp.consesus.tir.fasta") or die $!;
                                print TAGOUT ">",$key,"_",$tag_count,"_F","\n",$lseq,"\n";
                                print TAGOUT ">",$key,"_",$tag_count,"_R","\n",$rseq_ori,"\n";
                                print TIRALIGN $tsd_l,"\n",$tsd_r,"\n";
                                print TIRALIGN @alignw,@align;
                                print TIRALIGN substr($align[5],16,30),"\n";
                                open (MULEH,">>MULE.Hconfident.candidates.txt") or die $!;
                                print MULEH $key,"\t",$lindex,"\t",$rindex-9,"\n";
                                close MULEH;
                                open (MUTOUTH,">>MULE.Hconfident.candidates.seq.fasta") or die $!;
                                print MUTOUTH ">",$key,"_",$lindex,"_",$rindex-9,"\n";
                                my $lout=substr($t1,$lindex);
                                my $rout=substr($t2,0,$rindex-9);
                                print MUTOUTH $lout,$t,$rout,"\n";
                                $tag_count++;
                            }
                            elsif($total_star>=100)
                            {
                                open (MULEW,">>MULE.noconserved.tag.candidates.txt") or die $!;
                                print MULEW $key,"\t",$lindex,"\t",$rindex-9,"\n";
                                close MULEW;
                                open (MUTOUTW,">>MULE.noconserved.tag.candidates.seq.fasta") or die $!;
                                print MUTOUTW ">",$key,"_",$lindex,"_",$rindex-9,"\n";
                                my $lout=substr($t1,$lindex);
                                my $rout=substr($t2,0,$rindex-9);
                                print MUTOUTW $lout,$t,$rout,"\n";
                            }
                        }
                    }
                }
                @lmatch=();
                @rmatch=();
            }
        }
        @left_final=();
        @right_final=();
    }
    close TIRALIGN;
}
unlink "TSD1.fasta";
unlink "TIR.fasta";
unlink "left.jf";
unlink "right.jf";
unlink "temp1.fasta";
unlink "temp2.fasta";
unlink "tag.blast.out";
