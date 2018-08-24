#!/usr/bin/perl
my ($genome_file,$line,$tag_file,$key,$tir_l,$tir_r,$max1,$max2);
my (@tags,@genome,@flank,@pos,@left,@right,@checked_tags);
my (%tags,%gSeq,%pos,%chr,%hconfi,%gflank,%tir_tag);

my $input=shift @ARGV;
open(TAG,"$input") or die $!;
	@tags=<TAG>;

	foreach $line(@tags)
	{
		my @temp=split(" ",$line);
		if ($chr{$temp[0]})
		{
			
			$chr{$temp[0]}=$chr{$temp[0]}+1;
		}
		else
		{
			$chr{$temp[0]}=1;
		}
	}

	foreach $line(keys %chr)
	{
		print $line,"\n" if $chr{$line}>=9
	}
