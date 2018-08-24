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

foreach $key(keys %gSeq)
{
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
    
    foreach my $shared_kmer1(@left_final)
    {
        foreach my $shared_kmer2(@right_final)
        {
            @lmatch=();
            @rmatch=();
            my @kmatch=();
            open(KMER,">kmer.fasta") or die $!;
            print KMER ">left\n",$shared_kmer1,"\n",">right\n",$shared_kmer2,"\n";
            my @kalign=`mafft --quiet --clustalout --maxiterate 300 --localpair --thread 1 kmer.fasta`; 
            foreach my $kline(@kalign)
            {
                @kmatch=$kline=~/\*/g if $kline=~/\*/;
            }
            if($#kmatch>=7)
            {
                #print @kalign;
                while($t1=~/$shared_kmer1/g)
                {
                    push @lmatch,pos($t1);
                }
            
                while($t2=~/$shared_kmer2/g)
                {
                    push @rmatch,pos($t2);
                }
            
                foreach my $lindex(@lmatch)
                {
                    my $lseq=substr($t1,$lindex,100);
                    my $lflank=substr($t1,$lindex,1);
                    #my $tsd_l=substr($t1,$lindex-9,9);
                    my $lflank2=substr($t1,$lindex-10,1);
                
                    foreach my $rindex(@rmatch)
                    {
                        $checker=0;
                        open(OUT3,">TSD1.fasta") or die $!;
                        my $rflank=substr($t2,$rindex,1);
                        #my $tsd_r=substr($t2,$rindex-9,9);
                        my $rflank2=substr($t2,$rindex-10,1);
                        my $rseq=substr($t2,$rindex-110,100);
                        unless($lflank eq $rflank && $lflank2 eq $rflank2)
                        {
                            $rseq=reverse $rseq;
                            $rseq=~tr/ATCG/TAGC/;
                            print OUT3 ">left_TSD\n",$lseq,"\n";
                            print OUT3 ">right_TSD\n",$rseq,"\n";
                                                        
                            my @align=`mafft --quiet --clustalout --maxiterate 300 --localpair --thread 1 TSD1.fasta`;
                            my $total_star=0;
                        
                            foreach my $star_line(@align)
                            {
                                my @stare=$star_line=~/\*/g if $star_line=~/\*/;
                                $total_star=$total_star+$#stare+1 if $star_line=~/\*/;
                            }
                            if($total_star>=70)
                            {   
                                $tag1=substr($lseq,0,10);
                                $tag2=substr($rseq,0,10);
                                open(OUT12,">10bp.fasta") or die $!;
                                print OUT12 ">left_TSD\n",$tag1,"\n";
                                print OUT12 ">right_TSD\n",$tag2,"\n";
                                my @align12=`mafft --quiet --clustalout --maxiterate 300 --localpair --thread 1 10bp.fasta`; 
                                my @head10=();
                                foreach my $line12(@align12)
                                {
                                    @head10=$line12=~/\*/g if $line12=~/\*/;
                                }

                                $tag3=substr($lseq,0,60);
                                $tag4=substr($rseq,0,60);
        
                                my $timer=0;
                                my @timer=();
                            
                                @timer=$tag3=~/A/g;
                                $timer=1 if $#timer/60 >=0.4;
                                @timer=$tag3=~/T/g;
                                $timer=1 if $#timer/60 >=0.4; 
                                @timer=$tag3=~/C/g;
                                $timer=1 if $#timer/60 >=0.4;
                                @timer=$tag3=~/G/g;
                                $timer=1 if $#timer/60 >=0.4;

                                #print @align12,@align,$#head10,"\n";
                                if($#head10>=6 && $timer==0)
                                {   
                                    print TIRALIGN @align;
                                    open (MULEH,">>MULE.Hconfident.candidates.txt") or die $!;
                                    print MULEH $key,"\t",$lindex,"\t",$rindex-9,"\n";
                                    close MULEH;
                                    open (MUTOUTH,">>MULE.Hconfident.candidates.seq.fasta") or die $!;
                                    print MUTOUTH ">",$key,"_",$lindex,"_",$rindex-9,"\n";
                                    my $lout=substr($t1,$lindex);
                                    my $rout=substr($t2,0,$rindex-9);
                                    print MUTOUTH $lout,$t,$rout,"\n";
                                }
                                elsif($total_star>=80)
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
                }
            }
        }
    }
    @left_final=();
    @right_final=();
    close TIRALIGN;
}
unlink "TSD1.fasta";
unlink "TIR.fasta";
unlink "left.jf";
unlink "right.jf";
unlink "temp1.fasta";
unlink "temp2.fasta";
unlink "tag.blast.out";
