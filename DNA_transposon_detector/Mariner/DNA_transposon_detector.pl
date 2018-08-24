#!/usr/bin/perl -w
use strict;
use Cwd;
use Array::Utils qw(:all);
my($count,$num_threads,$config,$line,$dir,$DNAref,$Protref,$fasta_file,$blast,$DNA_transposon_detector,$einverted);
my($gFasta,$blast_out,$flanking_length,$temp1,$temp2,$TSD1,$TSD2,$job,$MULE_tag);
my(@config,@fasta,@inverted,@genome,@blast,@lmatch,@rmatch,@mule,@mariner,@LTR_repeat,@helitron,@cacta,@hat,@pong,@pif);
my(@temp1,@temp2,@temp3);
my(%gSeq,%temph,%tempn);

$config=shift @ARGV;
unless(defined $config)
{
        print "You must provide config file containing configuration info for LTR identification!\n";
}

# readin config file
open(CONFIG,"$config") or die $!;
@config=<CONFIG>;
foreach $line(@config)
{
        $line=~/(.+)\=(.+)/;
        if($1 eq "Protref")
        {
                $Protref=$2;
        }
        elsif($1 eq "fasta_file")
        {
                $fasta_file=$2;
        }
        elsif($1 eq "einverted")
        {
                $einverted=$2;
        }
        elsif($1 eq "blast")
        {
                $blast=$2;
        }
        elsif($1 eq "DNA_transposon_detector")
        {
                $DNA_transposon_detector=$2;
        }
        elsif($1 eq "einverted")
        {
                $einverted=$2;
        }
        elsif($1 eq "num_threads")
        {
                $num_threads=$2;
        }
        elsif($1 eq "flanking_length")
        {
                $flanking_length=$2;
        }
        elsif($1 eq "MULE_tag")
        {
                $MULE_tag=$2;
        }

}


# Build temp folder and rename fasta file header
$dir=getcwd;
if( -d "Temp" )
{
        system("rm -rf Temp");
        mkdir "Temp";
}
else
{
        mkdir "Temp";
}

if( -d "Results" )
{
        system("rm -rf Results");
        mkdir "Results";
}
else
{
        mkdir "Results";
}

unless(-d "BlastDB")
{
        mkdir "BlastDB";
}
open(FASTA,"$fasta_file") or die $!;
open(SINGLE,">Temp/formated.fasta") or die $!;
@fasta=<FASTA>;

for(my $i=0;$i<=$#fasta;$i++)
{
        chomp $fasta[$i];
        if($fasta[$i]=~/>/)
        {
                if($i==0)
                {
                        print SINGLE $fasta[$i],"\n";
                }
                else
                {
                        print SINGLE "\n",$fasta[$i],"\n";
                }
        }
        else
        {
                print SINGLE $fasta[$i];
        }
}
close FASTA;
close SINGLE;
@fasta=();


`$blast/makeblastdb  -in "$dir/Temp/formated.fasta" -input_type fasta -dbtype nucl -out $dir/BlastDB/Target_Genome`;
`$blast/tblastn -query "$Protref"  -db "$dir/BlastDB/Target_Genome"   -word_size 3 -evalue 1e-10 -num_threads $num_threads -outfmt 6 >./Temp/Gmax_mule_domain_tblastn.out`;

chdir "$dir/Temp";
open(GENOME,"formated.fasta") or die $!;
open(BLAST,"Gmax_mule_domain_tblastn.out") or die $!;
open(OUT,">DNA_transposon_hit_flanking.fasta") or die $!;
@genome=<GENOME>;
@blast=<BLAST>;

#read in all the fasta files in hash;
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
$count=1;
foreach my $line(@blast)
{
        my @pos=split(" ",$line);
        my $header=$pos[1];
        my $start;
        my $ends;
        if($pos[8]<$pos[9])
        {
                $start=$pos[8]-1;
                $ends=$pos[9];
        }
        elsif($pos[8]>$pos[9])
        {
                $start=$pos[9]-1;
                $ends=$pos[8];
        }
        my $length=$ends-$start+1;
        my $elength=length($gSeq{$header})-$ends;
        if($start>$flanking_length && $elength>$flanking_length)
        {
                my $seq=substr($gSeq{$header},$start-$flanking_length-1,$length+$flanking_length+$flanking_length);
                print OUT ">",$header,"_",$start,"_",$ends,"_",$ends-$start,"\n";
                print OUT $seq,"\n";
        }
        elsif($start>$flanking_length && $elength<$flanking_length)
        {
                my $seq=substr($gSeq{$header},$start-$flanking_length-1,$length+$flanking_length+$elength);
                print OUT ">",$header,"_",$start,"_",$ends,"_",$ends-$start,"\n";
                print OUT $seq,"\n";
        }
        elsif($start<$flanking_length && $elength>$flanking_length)
        {
                my $seq=substr($gSeq{$header},0,$length+$flanking_length+$start);
                print OUT ">",$header,"_",$start,"_",$ends,"_",$ends-$start,"\n";
                print OUT $seq,"\n";
        }
        elsif($start<$flanking_length && $elength<$flanking_length)
        {
                my $seq=substr($gSeq{$header},0,$length+$start+$elength);
                print OUT ">",$header,"_",$start,"_",$ends,"_",$ends-$start,"\n";
                print OUT $seq,"\n";
        }
}
@blast=();
@genome=();
%gSeq=();

#Thing need to be mannually done
#after running Mutator_TSD_TIR, get high confi, get first 100bp tag and muanlly check the best boarders, and copy and fasta.then run homolog based search.pl
#`makeblastdb  -in formated.fasta -input_type fasta -dbtype nucl -out Genome_seq`

