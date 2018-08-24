#!/usr/bin/perl -w
use strict;
my ($genome_file,$line,$tag_file,$key,$tir_l,$tir_r,$max1,$max2);
my (@tags,@genome,@flank,@pos,@left,@right,@checked_tags);
my (%tags,%gSeq,%pos,%chr,%hconfi,%gflank,%tir_tag,%flank);

$genome_file=shift @ARGV;
$tag_file=shift @ARGV;
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

open(CHECK,"$tag_file") or die $!;
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

foreach $key(keys %tir_tag)
{
	my @temp=();
	%pos=();
	%chr=();
	@temp=split(":",$tir_tag{$key});
	my $tir_F=$temp[0];
	my $tir_R=$temp[1];
	#print $tir_F,"\n",$tir_R,"\n";
	open(TARGET,">Paired.TIR.fasta") or die $!;
	print TARGET ">left_tir","\n",$tir_F,"\n";
	print TARGET ">right_tir","\n",$tir_R,"\n";
	close TARGET;
	##############################################################################################################THis part need to be modified#############################################################################################################
	`blastn -query Paired.TIR.fasta  -db /home/zyb/Wild_Soybean_Pangenome/Transposon_DB/MULE/Gmax/Temp/Genome_seq -num_threads 40  -outfmt 6 -word_size 8 >tag.blast.out`;
	##########################################################################################################################################################################################################################
	open(TAG,"tag.blast.out") or die $!;
	@tags=<TAG>;

	foreach $line(@tags)
	{
		my @temp=split(" ",$line);
		if ($chr{$temp[1]})
		{
			if($temp[8]<$temp[9])
			{
				my $position=$temp[8]."_".$temp[9];
				$chr{$temp[1]}=$chr{$temp[1]}.":".$position;
			}
			#else
			#{
				#my $position=$temp[9]."_".$temp[8];
				#$chr{$temp[1]}=$chr{$temp[1]}.":".$position;
			#}

		}
		else
		{
			if($temp[8]<$temp[9])
			{
				my $position=$temp[8]."_".$temp[9];
				$chr{$temp[1]}=$position;
			}
			#else
			#{
				#my $position=$temp[9]."_".$temp[8];
				#$chr{$temp[1]}=$position;
			#}
		}
	}

	foreach my $chromosome (keys %chr)
	{
		my $index=0;
		my %genome_coordinate=();
		my @coordinate=split(":",$chr{$chromosome});
		foreach my $gcor(@coordinate)
		{
			my @temp0=split("_",$gcor);
			$genome_coordinate{$temp0[0]}=$temp0[1];
		}
		my @start=sort{$a<=>$b} keys %genome_coordinate;
		#print $chromosome ,"\n";
		#print join("\t",@start),"\n";
		for(my $t=0;$t<$#start;$t++)
		{
			for(my $s=$t+1;$s<=$#start;$s++)
			{
				if($start[$s]-$start[$t] <20000 && $start[$s]-$start[$t]>100)
				{
					my $TSD_F=substr($gSeq{$chromosome},$start[$t]-25,30);
					my $TSD_R=substr($gSeq{$chromosome},$genome_coordinate{$start[$s]}-5,30);
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
						if ($temp1[3]-$temp1[4]>=8 && $temp1[3]<=11 && $temp1[8]<$temp1[9] && $temp1[8] <$temp1[9])
						{
							#print $tsd_line;
							#open(TSDD,">TSD_mactched.fasta") or die $!;
							#print TSDD ">tsd_subF\n",$TSD_F,"\n",">tsd_subR\n",$TSD_R,"\n";
							#close TSDD;
							#my @align0=`mafft --maxiterate 1000  --quiet --clustalout --localpair TSD_mactched.fasta`;
							#print @align0;
							#print $TSD_F,"\n",$TSD_R,"\n";

							my $TIR_F=substr($gSeq{$chromosome},$start[$t]-26+$temp1[8]+$temp1[3],60);
							my $TIR_R=substr($gSeq{$chromosome},$genome_coordinate{$start[$s]}-6+$temp1[6]-60,60);
							#print $TIR_F,"\n",$TIR_R,"\n";

							$TIR_R=reverse $TIR_R;
							$TIR_R=~tr/ATCG/TAGC/;
							open(TIRSEQ,">tirseq.fa") or die $!;
							print TIRSEQ ">tir_F\n",$TIR_F,"\n",">tir_R\n",$TIR_R,"\n";
							close TIRSEQ;
							my @align= `mafft --maxiterate 1000  --quiet --clustalout --localpair tirseq.fa`;
							my @tir_match=$align[5]=~/\*/g;
							my @head_star=substr($align[5],16,10)=~/\*/g;
							my $tsd_subF=substr($gSeq{$chromosome},$start[$t]-26+$temp1[8],$temp1[3]);
							my $tsd_subR=substr($gSeq{$chromosome},$genome_coordinate{$start[$s]}-6+$temp1[6],$temp1[3]);
							my $flank_F=substr($gSeq{$chromosome},$start[$t]-26+$temp1[8],$temp1[3]);
							my $flank_R=substr($gSeq{$chromosome},$genome_coordinate{$start[$s]}-6+$temp1[6],$temp1[3]);

							open(TSDD,">TSD_mactched.fasta") or die $!;
							print TSDD ">tsd_subF\n",$tsd_subF,"\n",">tsd_subR\n",$tsd_subR,"\n";
							close TSDD;
							my @align0=`mafft --maxiterate 1000  --quiet --clustalout --localpair TSD_mactched.fasta`;
							#print @align0;
							my @tsd_match=$align0[5]=~/\*/g;
							unless($flank{$flank_F} || $flank{$flank_R})
							{
								if($#tir_match>=34 && $#tsd_match>=7)
								{
									#print $tsd_line,"\n";
									#print $TSD_F,"\n";
									#print $TSD_R,"\n";
								
									#print @align0;
									open(MUTATOR,">>Mutator.cadidates.coordinate.txt") or die $!;
									print MUTATOR $chromosome,"_",$start[$t]-26+$temp1[8]+$temp1[3],"_",$genome_coordinate{$start[$s]}-6+$temp1[6],"\n";
									open(FALIGN,">>Mutator.candidates.align") or die $!;
									print FALIGN $chromosome,"_",$start[$t]-26+$temp1[8]+$temp1[3],"_",$genome_coordinate{$start[$s]}-6+$temp1[6],"\t",$#tsd_match+1,"\n",@align0,@align;
									#print $chromosome,"\t",$start[$t]+$temp1[9]-31,"\t",$genome_coordinate{$start[$s]}+$temp1[6],"\n";
									#print @align;
									open(FINAL,">>Final.identified.tir.fasta") or die $!;
									print FINAL ">",$chromosome,"_",$start[$t]-26+$temp1[8]+$temp1[3],"_",$genome_coordinate{$start[$s]}-6+$temp1[6],"_",$index,"_F\n",$TIR_F,"\n";
									print FINAL ">",$chromosome,"_",$start[$t]-26+$temp1[8]+$temp1[3],"_",$genome_coordinate{$start[$s]}-6+$temp1[6],"_",$index,"_R\n",$TIR_R,"\n";
									$flank{$flank_F}=1;
									$flank{$flank_R}=1;
									$index++;
								}
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


















