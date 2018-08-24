#!/usr/bin/perl -w
use strict;
my ($genome_file,$line,$key,$tir_l,$tir_r,$max1,$max2);
my (@tags,@genome,@flank,@pos,@left,@right,@checked_tags);
my (%tags,%gSeq,%pos,%chr,%ref_tag,%hconfi,%gflank,%tir_tag,%uniq);

$genome_file=shift @ARGV;
open(GEN,"$genome_file") or die $!;

@genome=<GEN>;
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

open(CHECK,"MULE.Hconfident.candidates.seq.manully.checked.fasta") or die $!;
@checked_tags=<CHECK>;

for(my $i=0;$i<$#checked_tags;$i+=2)
{
	chomp $checked_tags[$i+1];
	$checked_tags[$i]=~/(.+)_[FR]/;
	if($tir_tag{$1})
	{
		$tir_tag{$1}=$tir_tag{$1}.":".$checked_tags[$i+1];
	}
	else
	{
		$tir_tag{$1}=$checked_tags[$i+1];
	}

}

foreach my $ref_tag(keys %tir_tag)
{
	my @tag_ref=split(":",$tir_tag{$ref_tag});
	open(REFTAG,">ref.tag") or die $!;
	print REFTAG ">Tir_F\n",$tag_ref[0],"\n",">Tir_R\n",$tag_ref[1],"\n";
	close REFTAG;
	`makeblastdb  -in ref.tag -input_type fasta -dbtype nucl -out Tag`;
	`blastn -query $genome_file  -db Tag  -num_threads 40 -outfmt 6 -word_size 7   >g2ref.blastn.out`;


	open(GBLAST,"g2ref.blastn.out") or die $!;
	my @gblast=<GBLAST>;

	foreach $line(@gblast)
	{
	my @temp=split(" ",$line);
	if ($chr{$temp[0]})
	{
		if($temp[6]<$temp[7])
		{
			my $position=$temp[6]."_".$temp[7];
			$chr{$temp[0]}=$chr{$temp[0]}.":".$position;
		}
		else
		{
			my $position=$temp[7]."_".$temp[6];
			$chr{$temp[0]}=$chr{$temp[0]}.":".$position;
		}
	}
	else
	{
		if($temp[6]<$temp[7])
		{
			my $position=$temp[6]."_".$temp[7];
			$chr{$temp[0]}=$position;
		}
		else
		{
			my $position=$temp[7]."_".$temp[6];
			$chr{$temp[0]}=$position;
		}
	}
	}


	foreach my $chromosome (keys %chr)
	{
	my %genome_coordinate=();
	my @coordinate=split(":",$chr{$chromosome});
	foreach my $gcor(@coordinate)
	{
		my @temp0=split("_",$gcor);
		$genome_coordinate{$temp0[0]}=$temp0[1];
	}
	my @start=sort{$a<=>$b} keys %genome_coordinate;
	print join("\t",@start),"\n";
	for(my $t=0;$t<$#start;$t++)
	{
		for(my $s=$t+1;$s<=$#start;$s++)
		{
			if($start[$s]-$start[$t] <20000 && $start[$s]-$start[$t] >100)
			{
				my $TSD_F=substr($gSeq{$chromosome},$start[$t]-30,30);
				my $TSD_R=substr($gSeq{$chromosome},$genome_coordinate{$start[$s]},30);
				open(DB,">DB.fa") or die $!;
				print DB ">DB\n",$TSD_F,"\n";
				open(QUERY,">Query.fa") or die $!;
				print QUERY ">Query\n",$TSD_R,"\n";
				close DB;
				close QUERY;

				`makeblastdb  -in DB.fa -input_type fasta -dbtype nucl -out DB`;
				`blastn -query Query.fa  -db DB  -num_threads 44 -word_size 4 -gapextend 2 -outfmt 6 >TSD.blast.out`;

				open(TSDCHECK,"TSD.blast.out") or die $!;
				my @tsd_blast=<TSDCHECK>;
				foreach my $tsd_line(@tsd_blast)
				{
					my @temp1=split(" ",$tsd_line);
					if ($temp1[3]-$temp1[4]>=8 && $temp1[3]<=11 && $temp1[8]<$temp1[9])
					{
						my $TIR_F=substr($gSeq{$chromosome},$start[$t]+$temp1[9]-30,30);
						my $TIR_R=substr($gSeq{$chromosome},$genome_coordinate{$start[$s]}+$temp1[6]-30,30);
						$TIR_R=reverse $TIR_R;
						$TIR_R=~tr/ATCG/TAGC/;
						open(TIRSEQ,">tirseq.fa") or die $!;
						print TIRSEQ ">tir_F\n",$TIR_F,"\n",">tir_R\n",$TIR_R,"\n";
						close TIRSEQ;
						my @align= `mafft --maxiterate 1000  --quiet --clustalout --localpair tirseq.fa`;
						my @tir_match=$align[5]=~/\*/g;
						my $tsd_subF=substr($TSD_F,$temp1[8]-1,$temp1[9]-$temp1[8]+2);
						my $tsd_subR=substr($TSD_R,$temp1[6]-1,$temp1[7]-$temp1[6]+2);
						open(TSDD,">TSD_mactched.fasta") or die $!;
						print TSDD ">tsd_subF\n",$tsd_subF,"\n",">tsd_subR\n",$tsd_subR,"\n";
						close TSDD;
						my @align0=`mafft --maxiterate 1000  --quiet --clustalout --localpair TSD_mactched.fasta`;
						my @tsd_match=$align0[5]=~/\*/g;
						if($#tir_match>=34 && $#tsd_match>=9)
						{
							print $tsd_line,"\n";
							print $TSD_F,"\n";
							print $TSD_R,"\n";
								
							print @align0;
							open(MUTATOR,">>Mutator.cadidates.coordinate.txt") or die $!;
							print MUTATOR $chromosome,"\t",$start[$t]+$temp1[9]-31,"\t",$genome_coordinate{$start[$s]}+$temp1[6],"\n";
							print $chromosome,"\t",$start[$t]+$temp1[9]-31,"\t",$genome_coordinate{$start[$s]}+$temp1[6],"\n";
							print @align;
						}
					}
				}
			}
			else
			{
				last;
			}		
		}
	}
	}
}


















