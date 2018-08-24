#!/usr/bin/perl -w
my($input);
my(%index);
$input=shift @ARGV;
open(F,"$input") or die $!;

my @align=<F>;
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
		print $align[$j],"\n";
	}
}