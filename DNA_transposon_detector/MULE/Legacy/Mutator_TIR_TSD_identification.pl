#!/usr/bin/perl -w
use strict;
use Cwd;
use Array::Utils qw(:all);
use List::MoreUtils qw(uniq);

my($genomeseq,$key);
my(@lmatch,@rmatch,@genome);
my(%gSeq,%TIR_tag);


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

my $tag_count=0;
foreach $key(keys %gSeq)
{
    my @position1;
    my @position2;
    %TIR_tag=();
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
            my $shared_kmer_r=reverse $shared_kmer;
            $shared_kmer_r=~tr/ATCG/TAGC/;

        	while($t1=~/$shared_kmer/g)
    		{
       			push @position1,pos($t1);
    		}
            
    		while($t2=~/$shared_kmer_r/g)
    		{
        		push @position2,pos($t2);
    		}
        }
    }
    @position1=uniq(@position1);
    @position2=uniq(@position2);
    @position1=sort{$a<=>$b}@position1;
    @position2=sort{$a<=>$b}@position2;

    my @final_position1=();
    my @final_position2=();
    my $min=0;
    my @holder=();
    for (my $i=1;$i<=$#position1;$i++)
    {
        if($position1[$i]-$min<50)
        {
            push @holder,$position1[$i];
        }
        else
        {
            $min=$position1[$i];
        }
        
    }
    foreach my $ele(@position1)
    {
        my $finder=0;
        foreach my $ele2(@holder)
        {
            $finder=1 if $ele==$ele2;
        }
        push @final_position1,$ele unless $finder==1;
    }

    $min=0;
    @holder=();
    for (my $i=1;$i<$#position2;$i++)
    {
        if($position2[$i]-$min<50)
        {
            push @holder,$position2[$i];
        }
        else
        {
            $min=$position2[$i];
        }
    }
    foreach my $ele(@position2)
    {
        my $finder=0;
        foreach my $ele2(@holder)
        {
            $finder=1 if $ele==$ele2;
        }
        push @final_position2,$ele unless $finder==1;
    }
    foreach my $lindex(@final_position1)
    {
        my($lseq,$rseq);
        $lseq=substr($t1,$lindex-69,60);
        foreach my $rindex(@final_position2)
        {
            my $rseq=substr($t2,$rindex,60);
            open(OUT11,">temp11.fasta") or die $!;
            open(OUT12,">temp12.fasta") or die $!;
            print OUT11 ">Left\n",$lseq,"\n";
            print OUT12 ">right\n",$rseq,"\n";
            `jellyfish count -m 4 -s 100000 temp11.fasta  -o left11.jf`;
            `jellyfish count -m 4 -s 100000 temp12.fasta  -o right12.jf`;
            `jellyfish dump -c left11.jf >left11.kmer`;
            `jellyfish dump -c right12.jf >right12.kmer`;
            open(F11,"left11.kmer") or die $!;
            open(F12,"right12.kmer") or die $!;
            my @left11=<F11>;
            my @right12=<F12>;
            my @left_final11=();
            my @right_final12=();

            foreach my $kleft11(@left11)
            {
                $kleft11=~/(.+?)\s+/;
                push @left_final11, $1;
            }
            foreach my $kright12(@right12)
            {
                $kright12=~/(.+?)\s+/;
                push @right_final12, $1;
            }
            my @shared_sub=intersect(@left_final11,@right_final12);
            unless ($#shared_sub==-1)
            {    
                foreach my $shared_kmer11(@shared_sub)
                {
                    my @lmatch11=();
                    my @rmatch12=();

                    while($lseq=~/$shared_kmer11/g)
                    {
                        push @lmatch11,pos($lseq);
                    }
            
                    while($rseq=~/$shared_kmer11/g)
                    {
                        push @rmatch12,pos($rseq);
                    }
                    foreach my $lindex11(@lmatch11)
                    {
                        my($lseq_subL,$lseq_subR);
                        my($rseq_subL,$rseq_subR);
                        $lseq_subL=substr($lseq,$lindex11-9,9);
                        $lseq_subR=substr($lseq,$lindex11-4,9);
                        foreach my $rindex12(@rmatch12)
                        {   
                            my $timer=0;            
                            my $rseq_subL=substr($rseq,$rindex12-9,9);
                            my $rseq_subR=substr($rseq,$rindex12-4,9);
                            open(SUB1,">sub1.fasta") or die $!;
                            print SUB1 ">left_TSD\n",$lseq_subL,"\n";
                            print SUB1 ">right_TSD\n",$rseq_subL,"\n";
                            close SUB1;
                            my @align11=`mafft --quiet --clustalout --maxiterate 300 --localpair --thread 1 sub1.fasta`; 

                            open(SUB2,">sub2.fasta") or die $!;
                            print SUB2 ">left_TSD\n",$lseq_subR,"\n";
                            print SUB2 ">right_TSD\n",$rseq_subR,"\n";
                            close SUB2;
                            my @align12=`mafft --quiet --clustalout --maxiterate 300 --localpair --thread 1 sub2.fasta`;

                            my (@kmatch11,@kmatch12);
                            @kmatch11=$align11[5]=~/\*/g;
                            @kmatch12=$align12[5]=~/\*/g;
                                    
                            if($#kmatch11 >=7)
                            {
                                #print @align11;
                                my $subseq11=substr($t1,$lindex-69+$lindex11,120);
                                my $subseq12=substr($t2,$rindex+$rindex12-129,120);
                                my $left_final_tsd=substr($t1,$lindex-78+$lindex11,9);
                                my $right_final_tsd=substr($t2,$rindex+$rindex12-9,9);
                                #print $left_final_tsd,"\n",$right_final_tsd,"\n",$subseq11,"\n",$subseq12,"\n";
                                $subseq12=reverse $subseq12;
                                $subseq12=~tr/ATCG/TAGC/;
                                open(FINAL1,">final_tir.fasta") or die $!;
                                print FINAL1 ">left_TSD\n",$subseq11,"\n";
                                print FINAL1 ">right_TSD\n",$subseq12,"\n";
                                close FINAL1;
                                my @tir11=`mafft --quiet --clustalout --maxiterate 300 --localpair --thread 1 final_tir.fasta`;
                                my $final_match=0;
                                my @head_star=();
                                my @tail_star=();
                                print @tir11;
                                @head_star=substr($tir11[5],16,30)=~/\*/g;
                                #@tail_star=substr($tir11[5],47,30)=~/\*/g;
                                #print $#head_star,"\n",$#tail_star,"\n";
                                my @timer=();
                            
                                @timer=$subseq11=~/A/g;
                                $timer=1 if $#timer/60 >=0.4;
                                @timer=$subseq11=~/T/g;
                                $timer=1 if $#timer/60 >=0.4; 
                                @timer=$subseq11=~/C/g;
                                $timer=1 if $#timer/60 >=0.4;
                                @timer=$subseq11=~/G/g;
                                $timer=1 if $#timer/60 >=0.4;
                                if ($#head_star>=22 && $timer==0)
                                {
                                    unless($TIR_tag{$subseq11} && $TIR_tag{$subseq12})
                                    {
                                        open(FINALOUT,">>Candidate.Mutators.align") or die $!;
                                        print FINALOUT "$key\n",@tir11;
                                        close FINALOUT;
                                        open(TAGOUT,">>Candidate.TIRs.fasta") or die $!;
                                        print TAGOUT ">left.$key",$tag_count,"\n",$subseq11,"\n",">right.$key",$tag_count,"\n",$subseq12,"\n";
                                        $tag_count++;
                                        $TIR_tag{$subseq11}=1;
                                        $TIR_tag{$subseq12}=1
                                    }
                                }
                            }
                            elsif($#kmatch12>=7 )
                            {
                                #print @align12;
                                my $subseq11=substr($t1,$lindex-64+$lindex11,120);
                                my $subseq12=substr($t2,$rindex+$rindex12-124,120);
                                my $left_final_tsd=substr($t1,$lindex-73+$lindex11,9);
                                my $right_final_tsd=substr($t2,$rindex+$rindex12-4,9);
                                #print $left_final_tsd,"\n",$right_final_tsd,"\n",$subseq11,"\n",$subseq12,"\n";
                                $subseq12=reverse $subseq12;
                                $subseq12=~tr/ATCG/TAGC/;
                                open(FINAL2,">final_tir.fasta") or die $!;
                                print FINAL2 ">left_TSD\n",$subseq11,"\n";
                                print FINAL2 ">right_TSD\n",$subseq12,"\n";
                                close FINAL2;
                                my @tir12=`mafft --quiet --clustalout --maxiterate 300 --localpair --thread 1 final_tir.fasta`;
                                my $final_match=0;
                                my @head_star=();
                                my @tail_star=();
                                #print @tir12;
                                @head_star=substr($tir12[5],16,30)=~/\*/g;
                                #@tail_star=substr($tir12[5],47,30)=~/\*/g;
                                #print $#head_star,"\n",$#tail_star,"\n";

                                my @timer=();
                            
                                @timer=$subseq11=~/A/g;
                                $timer=1 if $#timer/60 >=0.4;
                                @timer=$subseq11=~/T/g;
                                $timer=1 if $#timer/60 >=0.4; 
                                @timer=$subseq11=~/C/g;
                                $timer=1 if $#timer/60 >=0.4;
                                @timer=$subseq11=~/G/g;
                                $timer=1 if $#timer/60 >=0.4;
                                if ($#head_star>=22 && $timer==0)
                                {
                                    unless($TIR_tag{$subseq11} && $TIR_tag{$subseq12})
                                    {
                                        open(FINALOUT,">>Candidate.Mutators.align") or die $!;
                                        print FINALOUT "$key\n",@tir12;
                                        close FINALOUT;
                                        open(TAGOUT,">>Candidate.TIRs.fasta") or die $!;
                                        print TAGOUT ">left.$key",$tag_count,"\n",$subseq11,"\n",">right.$key",$tag_count,"\n",$subseq12,"\n";
                                        $tag_count++;
                                        $TIR_tag{$subseq11}=1;
                                        $TIR_tag{$subseq12}=1;
                                    }
                                }
                            }
                        }
                    }    
                }
            }
        }
    }
            
}




