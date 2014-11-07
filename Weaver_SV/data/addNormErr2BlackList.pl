#!/usr/bin/perl
#
#open(I,"<",shift@ARGV);
while(<STDIN>){
	#10      127190672       +       10      127197228       +       27      7
	chomp;
	@m = split(/\s+/);
	if($m[0] ne $m[3] || ($m[0] eq $m[3] && $m[2] eq $m[5]) || abs($m[1] - $m[4]) > 10000000){
		$m[0] = "chr".$m[0];
		$m[3] = "chr".$m[3];
		$#m = 5;
		print join("\t",@m),"\n";
		#print 
	}

}

