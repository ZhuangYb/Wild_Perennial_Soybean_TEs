#!/usr/bin/perl -w
use strict;
my ($genome_file,$line,$key,$tir_l,$tir_r,$max1,$max2);
my (@tags,@genome,@flank,@pos,@left,@right,@checked_tags);
my (%tags,%gSeq,%pos,%cor,%hconfi,%gflank,%tir_tag);

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

foreach $key(keys %tir_tag)
{
	my @temp=();
	@temp=split(":",$tir_tag{$key});
	my $tir_F=substr($temp[0],0,30);
	my $tir_R=substr($temp[1],length($temp[1])-30,30);
		
	open(TARGET,">Paired.TIR.fasta") or die $!;
	print TARGET ">left_tir","\n",$tir_F,"\n";
	print TARGET ">right_tir","\n",$tir_R,"\n";
	close TARGET;
	`blastn -query Paired.TIR.fasta  -db Genome_seq -num_threads 40  -outfmt 6 -word_size 8 >tag.blast.out`;
	open(TAG,"tag.blast.out") or die $!;
	@tags=<TAG>;

	foreach $line(@tags)
	{
		my @temp=split(" ",$line);
		if(defined $tags{$temp[1]})
		{
			$tags{$temp[1]}=$tags{$temp[1]}.":".$temp[8]."_".$temp[9] if $temp[6]<=3 ;
		}
		else
		{
			$tags{$temp[1]}=$temp[8]."_".$temp[9] if $temp[6]<=3 ;
		}
		$cor{$temp[8]."_".$temp[9]}=$temp[0] if $temp[6]<=3 ;
	}
	my $count=1;
	foreach $key(keys %tags)
	{
		%pos=();
		my @order=();
		my @pos=split(":",$tags{$key});
		unless($#pos==0)
		{
			foreach $line(@pos)
			{
				my @temp=split("_",$line);
				$pos{$temp[0]}=$temp[1];
				push @order,$temp[0];
			}
			@order=sort {$a <=> $b} @order;
			for(my $i=0;$i<$#order;$i++)
			{	
				for(my $j=$i+1;$j<=$#order;$j++)
				{
					my $check_rt1=1;
					my $check_rt2=2;
					open(TEMP,">TSD.check.fasta") or die $!;
					if($order[$i]<$pos{$order[$i]})
					{
						$check_rt1=1;
						my $tsd_l=substr($gSeq{$key},$order[$i]-12,11);
						print TEMP ">left_TSD\n",$tsd_l,"\n";
						$max1=$pos{$order[$i]};
					}
					else
					{
						$check_rt1=2;
						my $tsd_l=substr($gSeq{$key},$pos{$order[$i]}-12,11);
						print TEMP ">left_TSD\n",$tsd_l,"\n";
						$max1=$order[$i];
					}

					if($order[$j]<$pos{$order[$j]})
					{
						$check_rt2=1;
						my $tsd_l=substr($gSeq{$key},$pos{$order[$j]}-1,11);
						print TEMP ">right_TSD\n",$tsd_l,"\n";
						$max2=$pos{$order[$j]};
					}
					else
					{
						$check_rt2=2;
						my $tsd_l=substr($gSeq{$key},$order[$j]-1,11);
						print TEMP ">right_TSD\n",$tsd_l,"\n";
						$max2=$order[$j];
					}
				
					if($cor{$order[$i]."_".$pos{$order[$i]}} ne $cor{$order[$j]."_".$pos{$order[$j]}}  && abs($max1-$max2)<20000 && $check_rt1==$check_rt2)
					{
						my @align=`mafft --quiet --clustalout --maxiterate 300 --localpair --thread 40 TSD.check.fasta`;
						foreach my $align(@align)
						{
							my $move_l=0;
							my $move_r=0;
							my @star=();
							if($align=~/\*/)
							{
								@star=$align=~/\*/g;
							}
							elsif($align=~/left/)
							{
								my @temp1=split(" ",$align);
								@left=split("",$temp1[1]);
							}
							elsif($align=~/right/)
							{
								my @temp2=split(" ",$align);
								@right=split("",$temp2[1]);
							}
							if($#star>=7)
							{
								for(my $a=$#left;$a>=$#left-3;$a--)
								{
									unless($left[$a] eq "-")
									{
										if($left[$a] eq $right[$a])
										{
											if($left[$a-1] eq "-" && $left[$a-2] eq "-")
											{
												$move_l++;
											}
											else
											{
												last;
											}
										}
										else
										{
											$move_l++;
										}
									}
								}
								for(my $b=0;$b<=3;$b++)
								{
									unless($right[$b] eq "-")
									{
										if($right[$b] eq $left[$b])
										{
											if($right[$b+1] eq "-" && $right[$b+2] eq "-")
											{
												$move_r++;
											}
											else
											{
												last;
											}
										}
										else
										{
											$move_r++;
										}
									}
								}
								if($check_rt1==1)
								{
									$tir_l=substr($gSeq{$key},$order[$i]-1-$move_l,100);
									$tir_r=substr($gSeq{$key},$pos{$order[$j]}-101+$move_r,100);
									$tir_r=reverse $tir_r;
									$tir_r=~tr/ATCG/TAGC/;
									open(TIR,">>TIR.tag.fasta") or die $!;
									print TIR ">LTIR.$count\n",$tir_l,"\n",">RTIR.$count\n",$tir_r,"\n";
									open(COR,">>MULE.genome.coordinate.txt") or die $!;
									print COR $key,"\t",$order[$i]-1-$move_l,"\t",$pos{$order[$j]}-1+$move_r,"\n";
									$count++
								}
								elsif($check_rt1==2)
								{
									$tir_l=substr($gSeq{$key},$pos{$order[$i]}-1-$move_l,100);
									$tir_r=substr($gSeq{$key},$order[$j]-101+$move_r,100);
									$tir_r=reverse $tir_r;
									$tir_r=~tr/ATCG/TAGC/;
									open(TIR,">>TIR.tag.fasta") or die $!;
									print TIR ">LTIR.$count\n",$tir_l,"\n",">RTIR.$count\n",$tir_r,"\n";
									my @align2=`mafft --quiet --clustalout --maxiterate 300 --localpair --thread 40 TSD.check.fasta`;
								
									open(COR,">>MULE.genome.coordinate.txt") or die $!;
									print COR $key,"\t",$pos{$order[$i]}-1-$move_l,"\t",$order[$j]-1+$move_r,"\n";
									$count++;
								}
							}
						}
					}
				}	
			}
		}
	}
} 


















