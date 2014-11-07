#!/usr/bin/perl
#
$number_flag = shift@ARGV; # 1 for chr 0 for number
open(I,"<",shift@ARGV);
#chr20   62920098        62923146        chr4    120324226       120327271
while(<I>){
	chomp;
	@m = split(/\s+/);
	if($number_flag == 0){
		$m[0] = substr($m[0],3);
		$m[3] = substr($m[3],3);
	}
	push @{$hash{"$m[0]\t$m[3]"}}, [@m];
}

#open(I,"<",shift@ARGV);
#chr2    87983137        +       chr2    112024779       -       3
while(<STDIN>){
	chomp;
	@M = split(/\s+/);
	@m = ();
#	print join("\t",@M),"\n";
	if($M[2] ne "+" && $M[2] ne "-"){
		$m[0] = $M[0];
		$m[3] = $M[4];
		if($M[3] == "+"){
			$m[1] = $M[2];
		}
		else{
			$m[1] = $M[1];
		}
		if($M[7] == "+"){
			$m[4] = $M[6];
		}
		else{
			$m[4] = $M[5];
		}
		$m[2] = $M[3];
		$m[5] = $M[7];
		$m[6] = $M[$#M-1];
		$m[7] = $M[$#M];
	}
	else{
		@m = @M;
	}
	$F = 0;
	if(exists $hash{"$m[0]\t$m[3]"}){
		for $i (0 .. $#{$hash{"$m[0]\t$m[3]"}}){
			if($m[1] <= ${${$hash{"$m[0]\t$m[3]"}}[$i]}[2]+100 && $m[1] >= ${${$hash{"$m[0]\t$m[3]"}}[$i]}[1]-100 && $m[4] <= ${${$hash{"$m[0]\t$m[3]"}}[$i]}[5]+100 && $m[4] >= ${${$hash{"$m[0]\t$m[3]"}}[$i]}[4]-100){
				$F = 1;
				last;
			}
		}
	}
	if($F != 1){
		#print $_,"\n";
		print join("\t",@m),"\n";
	}

}
