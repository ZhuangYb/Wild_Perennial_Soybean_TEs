#!/usr/bin/perl -w
use strict;
use Cwd;
use Array::Utils qw(:all);
use List::MoreUtils qw(uniq);

my($genomeseq,$key);
my(@lmatch,@rmatch,@genome);
my(%gSeq);


$genomeseq=shift @ARGV;
open(GENO,"$genomeseq") or die $!;
@genome=<GENO>;

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
    my $t2_r;
    my $t; #internal seq
    my $lmatch;
    my $rmatch;
    my $flength=length($gSeq{$key});
    open(OUT1,">temp1.fasta") or die $!;
    open(OUT2,">temp2.fasta") or die $!;
    $t1=substr($gSeq{$key},0,10000);
    $t=substr($gSeq{$key},10000,$match_l);
    $t2=substr($gSeq{$key},10000+$match_l,$flength-10000-$match_l);
    $t2_r=reverse $t2;
    $t2_r=~tr/ATCG/TAGC/;
    print OUT1 ">Left\n",$t1,"\n";
    print OUT2 ">right\n",$t2_r,"\n";
    `jellyfish count -m 9 -s 100000 temp1.fasta  -o left.jf`;
    `jellyfish count -m 9 -s 100000 temp2.fasta  -o right.jf`;
    `jellyfish dump -c left.jf >left.kmer`;
    `jellyfish dump -c right.jf >right.kmer`;
    unlink "left.jf";
    unlink "right.jf";
    unlink "temp1.fasta";
    unlink "temp2.fasta";
    close OUT1;
    close OUT2;

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
            my $shared_kmer_r=reverse $shared_kmer;
            $shared_kmer_r=~tr/ATCG/TAGC/;

        	while($t1=~/$shared_kmer/g)
    		{
       			push @lmatch,pos($t1);
    		}
            
    		while($t2=~/$shared_kmer_r/g)
    		{
        		push @rmatch,pos($t2);
    		}

    		foreach my $lindex(@lmatch)
    		{
        		my($lseq,$rseq);
        		$lseq=substr($t1,$lindex-159,300);
        		foreach my $rindex(@rmatch)
        		{
                    my $timer=0;
        			$rseq=substr($t2,$rindex-159,300);

        			open(TIR1,">TIR.candidate.fasta") or die $!;
        			print TIR1 ">Combined_seq\n",$lseq,$rseq,"\n";
        			close TIR1;

                    open(TIR2,">TIR2.candidate.fasta") or die $!;
                    my $rseq_r=reverse $rseq;
                    $rseq_r=~tr/ATCG/TAGC/;
                    print TIR2 ">left_seq\n",$lseq,"\n",">right_seq\n",$rseq_r,"\n";
                    close TIR2;
                    my @aligntir=`mafft --quiet --clustalout --maxiterate 300 --localpair --thread 1 TIR2.candidate.fasta`;
                    print @aligntir;

        			my @quiet=`einverted  -maxrepeat 20000 -gap 12 -threshold  50 -match 3 -mismatch -4 -outfile Temp.invert -outseq Temp.TIR -sequence TIR.candidate.fasta `;
        			my $combined=$lseq.$rseq;
                    open(INPUT2,"Temp.invert");
                    my @invert=<INPUT2>;
        			for(my $count2=0;$count2<=$#invert;$count2++)
        			{
                        print $invert[$count2];
        				if($invert[$count2]=~/Score/)
        				{
                            #open(INPUT1,"TIR.candidate.fasta");
                            #foreach my $input1(<INPUT1>)
                            #{
                                #print $input1;
                            #}
        					my @start_F=split(" ",$invert[$count2+1]);
        					my @start_R=split(" ",$invert[$count2+3]);
                            my $tir_length=$start_F[2]-$start_F[0]+1;
                            my @timer=();
                            @timer=$start_F[1]=~/A/g;
                            $timer=1 if $#timer/$tir_length >=0.4;
                            @timer=$start_F[1]=~/T/g;
                            $timer=1 if $#timer/$tir_length >=0.4; 
                            @timer=$start_F[1]=~/C/g;
                            $timer=1 if $#timer/$tir_length >=0.4;
                            @timer=$start_F[1]=~/G/g;
                            $timer=1 if $#timer/$tir_length >=0.4;

        					if($timer==0 && $tir_length>=30 && $start_F[2]<300 && $start_R[2]>00)
                            {
                                open(TSD,">temp.tsd.fasta") or die $!;
                                print $start_F[0],"\t",$start_R[0],"\n";
                                my $tsd_l=substr($combined,$start_F[0]-12,11);
                                my $tsd_r=substr($combined,$start_R[0],11);
                                print TSD ">Left\n",$tsd_l,"\n",">right\n",$tsd_r,"\n";
                                close TSD;
                                my @align=`mafft --quiet --clustalout --maxiterate 300 --localpair --thread 1 temp.tsd.fasta`; 
                                open(ALIGN1,">>TSD.align");
                                print ALIGN1 @align;
                                close ALIGN1;
                            }
        				}
        			}
        		}
        	}
    	}
    }
}












































