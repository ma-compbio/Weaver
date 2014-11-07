#!/usr/bin/perl
#
#open(I,"<",shift@ARGV);
while(<STDIN>){
	chomp;
	@m = split(/\s+/);
	if($m[0]  eq "all"){
		next;
	}
	$hash{"$m[0]\t$m[1]\t$m[2]"} += $m[7]*$m[8]/($m[2]-$m[1]+1);
}

foreach $key (keys %hash){
	print $key,"\t",$hash{$key},"\n";
}
