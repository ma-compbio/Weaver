#!/usr/bin/perl
#
#open(I,"<",shift@ARGV);
$L = shift@ARGV;
while(<STDIN>){
	chomp;
	@m = split(/\s+/);
	if($m[0] eq $m[3] && $m[2] ne $m[5] && abs($m[1] - $m[4]) < $L){
	}
	else{
		print $_,"\n";
	}

}

