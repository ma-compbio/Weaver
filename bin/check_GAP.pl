#!/usr/bin/perl
#
open(I,"<",shift@ARGV);
while(<I>){
	chomp;
	@m = split(/\s+/);
	if($m[$#m] eq "region"){
		#push @{$region{$m[0]}}, $_;
		$region_b{$m[0]} = $m[1];
		$region_e{$m[0]} = $m[2];
	}
	if($m[3] eq "GAP" || $m[3] eq "Del"){
		push @{$hash_b{$m[0]}}, $m[1];
		push @{$hash_e{$m[0]}}, $m[2];
	}
}

#open(I,"<",shift@ARGV);
while(<STDIN>){
	chomp;
	@m = split(/\s+/);
	$FLAG = 0;
	if(substr($m[0],0,1) ne "c"){
		$m[0] = "chr".$m[0];
	}
	if(!exists $region_b{$m[0]}){
		next;
	}
	if($m[1] >= $region_b{$m[0]} && $m[1] <= $region_e{$m[0]}){
		for $i (0 .. $#{$hash_b{$m[0]}}){
			if($hash_b{$m[0]}[$i]-2 <= $m[1] && $hash_e{$m[0]}[$i]+2 >= $m[1]){
				$FLAG = 1;
				last;
			}
		}
		if($FLAG == 1){
			next;
		}
		if($FLAG == 0){
			print $_,"\n";
		}
	}

}
