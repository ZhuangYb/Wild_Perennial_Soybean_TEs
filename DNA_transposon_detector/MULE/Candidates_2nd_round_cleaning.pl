#!/usr/bin/perl
my ($input);
my (@align);
my (%pos);

$input=shift @ARGV;
open(F,"$input") or die $!;
@align=<F>;

for($i=0;$i<$#align-10;$i++)
{
	if($align[$i]=~/^(.+arrow.+)\s+\d+/)
	{
		$pos{$1}=$align[$i+12]
	}
}

open(T1,">>temp1.txt") or die $!;
foreach my $key(keys %pos)
{
	print  T1 $key,"\t",$pos{$key};
}


`bash Sort.sh`;

open(T2,"temp2.txt") or die $!;
@align=<T2>;

open(HIGH,">>confident.list") or die $!;
open(CHECK,">>need2check.list") or die $!;
for($i=1;$i<$#align-1;$i++)
{
	my $check=0;
	my @temp1=split(" ",$align[$i-1]);
	my @temp2=split(" ",$align[$i]);
	my @temp3=split(" ",$align[$i+1]);
	if ($temp1[0] eq $temp2[0])
	{
		if(abs($temp1[5]-$temp2[5])<200)
		{
			print CHECK $align[$i];
			$check=1;
		}
	}

	if ($temp2[0] eq $temp3[0])
	{
		if(abs($temp2[5]-$temp3[5])<200)
		{
			print CHECK $align[$i];
			$check=1;
		}
	}
	if($check==0)
	{
		print HIGH $align[$i];
	}
}

open(F,"need2check.list");
@align=<F>;
chomp @align;

	
for($j=1;$j<$#align-1;$j++)
{
	my $check=0;
	my @temp1=split(" ",$align[$j]);
	my @temp2=split(" ",$align[$j+1]);
	my @temp3=split("\t",$align[$j]);
	my @temp4=split("\t",$align[$j+1]);
	my $s1=substr($temp3[$#temp3],0,20);
	my $s2=substr($temp4[$#temp4],0,20);

	if($temp1[0] eq $temp2[0] && abs($temp1[5]-$temp2[5])<200)
	{

		my @m1=$s1=~/\*/g;
		my @m2=$s2=~/\*/g;
		if($#m1>=$#m2)
		{
			$index{$j+1}=1;
		}
		else
		{
			$index{$j}=1;
		}
	}
}

for($j=0;$j<$#align;$j++)
{
	unless($index{$j})
	{
		print HIGH $align[$j],"\n";
	}
}
close HIGH;
close CHECK;

open(F,"confident.list");
@align=<F>;
open(MULE,">Mutator.identified.list");
foreach my $line(@align)
{
	my @temp=split(" ",$line);
	print MULE join("_",@temp[0..4]),":",$temp[5],"..",$temp[6],"\n";
}
unlink "confident.list";
unlink "need2check.list";
