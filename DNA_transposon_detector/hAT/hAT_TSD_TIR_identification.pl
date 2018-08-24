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
    #conserved motif CAGTGGCG
    while ($gSeq{$key}=~/[ATCG][ATCG]GTGGCG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }while ($gSeq{$key}=~/[ATCG]A[ATCG]TGGCG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }while ($gSeq{$key}=~/[ATCG]AG[ATCG]GGCG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }while ($gSeq{$key}=~/[ATCG]AGT[ATCG]GCG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }while ($gSeq{$key}=~/[ATCG]AGTG[ATCG]CG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }while ($gSeq{$key}=~/[ATCG]AGTGG[ATCG]G/g)
    {
        push @tir_F,pos($gSeq{$key});
    }while ($gSeq{$key}=~/[ATCG]AGTGGC[ATCG]/g)
    {
        push @tir_F,pos($gSeq{$key});
    }



    while ($gSeq{$key}=~/C[ATCG][ATCG]TGGCG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/C[ATCG]G[ATCG]GGCG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/C[ATCG]GT[ATCG]GCG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/C[ATCG]GTG[ATCG]CG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/C[ATCG]GTGG[ATCG]G/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/C[ATCG]GTGGC[ATCG]/g)
    {
        push @tir_F,pos($gSeq{$key});
    }




    while ($gSeq{$key}=~/CA[ATCG][ATCG]GGCG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CA[ATCG]T[ATCG]GCG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CA[ATCG]TG[ATCG]CG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CA[ATCG]TGG[ATCG]G/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CA[ATCG]TGGC[ATCG]/g)
    {
        push @tir_F,pos($gSeq{$key});
    }




    while ($gSeq{$key}=~/CAG[ATCG][ATCG]GCG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CAG[ATCG]G[ATCG]CG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CAG[ATCG]GG[ATCG]G/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CAG[ATCG]GGC[ATCG]/g)
    {
        push @tir_F,pos($gSeq{$key});
    }




    while ($gSeq{$key}=~/CAGT[ATCG][ATCG]CG/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CAGT[ATCG]G[ATCG]G/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CAGT[ATCG]GC[ATCG]/g)
    {
        push @tir_F,pos($gSeq{$key});
    }


    while ($gSeq{$key}=~/CAGTG[ATCG][ATCG]G/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CAGTG[ATCG]C[ATCG]/g)
    {
        push @tir_F,pos($gSeq{$key});
    }

    while ($gSeq{$key}=~/CAGTGG[ATCG][ATCG]/g)
    {
        push @tir_F,pos($gSeq{$key});
    }
   

    
    while ($gSeq{$key}=~/[ATCG][ATCG]CCACTG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/[ATCG]G[ATCG]CACTG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/[ATCG]GC[ATCG]ACTG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/[ATCG]GCC[ATCG]CTG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/[ATCG]GCCA[ATCG]TG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/[ATCG]GCCAC[ATCG]G/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/[ATCG]GCCACT[ATCG]/g)
    {
        push @tir_R,pos($gSeq{$key});
    }




    while ($gSeq{$key}=~/C[ATCG][ATCG]CACTG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/C[ATCG]C[ATCG]ACTG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/C[ATCG]CC[ATCG]CTG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/C[ATCG]CCA[ATCG]TG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/C[ATCG]CCAC[ATCG]G/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/C[ATCG]CCACT[ATCG]/g)
    {
        push @tir_R,pos($gSeq{$key});
    }


    while ($gSeq{$key}=~/CG[ATCG][ATCG]ACTG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CG[ATCG]C[ATCG]CTG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CG[ATCG]CA[ATCG]TG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CG[ATCG]CAC[ATCG]G/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CG[ATCG]CACT[ATCG]/g)
    {
        push @tir_R,pos($gSeq{$key});
    }




    while ($gSeq{$key}=~/CGC[ATCG][ATCG]CTG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CGC[ATCG]A[ATCG]TG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CGC[ATCG]AC[ATCG]G/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CGC[ATCG]ACT[ATCG]/g)
    {
        push @tir_R,pos($gSeq{$key});
    }




    while ($gSeq{$key}=~/CGCC[ATCG][ATCG]TG/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CGCC[ATCG]C[ATCG]G/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CGCC[ATCG]CT[ATCG]/g)
    {
        push @tir_R,pos($gSeq{$key});
    }




    while ($gSeq{$key}=~/CGCCA[ATCG][ATCG]G/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
    while ($gSeq{$key}=~/CGCCA[ATCG]T[ATCG]/g)
    {
        push @tir_R,pos($gSeq{$key});
    }



    while ($gSeq{$key}=~/CGCCAC[ATCG][ATCG]/g)
    {
        push @tir_R,pos($gSeq{$key});
    }
  
    open (ALIGN,">>tir_tsd.alignment.txt") or die $!;
    open (SEQ,">>hAT.identifed.fasta") or die $!;
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
            if($pos1<10000 &&  $pos2>10000+$match_l)
            {
                my $tsd_F=substr($gSeq{$key},$pos1-16,8);
                my $tsd_R=substr($gSeq{$key},$pos2,8);
                my $tir_F=substr($gSeq{$key},$pos1-8,8);
                my $tir_R=substr($gSeq{$key},$pos2-8,8);
                $tir_R=reverse $tir_R;
                $tir_R=~tr/ATCG/TAGC/;
                open(TSD,">TSD.candidate.txt") or die $!;
                print TSD ">tsd_F\n",$tsd_F,"\n";
                print TSD ">tsd_R\n",$tsd_R,"\n";
                close TSD;
                my @align=`mafft --quiet --maxiterate 1000 --clustalout --localpair --thread 1 TSD.candidate.txt`;
                my @match=$align[5]=~/\*/g;
                
                open(TIR,">TIR.candidate.txt") or die $!;
                print TIR ">tir_F\n",$tir_F,"\n";
                print TIR ">tir_R\n",$tir_R,"\n";
                my @align1=`mafft --quiet --maxiterate 1000 --clustalout --localpair --thread 1 TIR.candidate.txt`;
                my @match1=$align1[5]=~/\*/g;

                if($#match>=6 && $#match1>=4)
                {
                    my $seq=substr($gSeq{$key},$pos1-8,$pos2-$pos1+9);
                    print SEQ ">$key","_",$pos1-8,"_",$pos2,"\n";
                    print SEQ $seq,"\n";
                    print ALIGN "$key","\t",$pos1,"\t",$pos2,"\n";
                    print ALIGN "$tir_F\t","$tir_R\n";
                    print ALIGN @align1;
                    print ALIGN @align,"\n\n\n";
                }
            }
        }
    }
    
}