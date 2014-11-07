#!/usr/bin/perl
#

#$NUMBER_FLAG = shift@ARGV; ## 1 for chr 0 for number?
$GAPFILE = shift@ARGV;
$SVFILE = shift@ARGV;

open(I,"<",$SVFILE);
$S = <I>;
if(substr($S,0,1) eq "c"){
	$NUMBER_FLAG = 1;
}
else{
	$NUMBER_FLAG = 0;
}
open(I,"<",$GAPFILE);
while(<I>){
	chomp;
	@m = split(/\s+/);
	if($m[$#m] eq "region"){
		if($NUMBER_FLAG == 1){
			push @{$region{$m[0]}}, $_;
		}
		else{
			push @{$region{substr($m[0],3)}}, $_;
		}

	}
	if($m[$#m] eq "GAP"){# || $m[$#m] eq "Del"){
		if($NUMBER_FLAG == 1){
			push @{$GAP{$m[0]}}, $_;
		}
		else{
			push @{$GAP{substr($m[0],3)}}, $_;
		}
	}
}

open(I,"<",$SVFILE);
while(<I>){
	chomp;
	@M = split(/\s+/);
	@m = ();
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
	$F_1 = 0;
	$F_2 = 0;
	if(exists $region{$m[0]}){
		for $i (0 .. $#{$region{$m[0]}}){
			@s = split(/\s+/, $region{$m[0]}[$i]);
#			$s[1] = $s[1] - 5000;# more region ban
#			$s[2] = $s[2] + 5000; #
			if($s[1] < $m[1] && $s[2] > $m[1]){
				$F_1 = 1;
			}
			#if($s[1] < $m[4] && $s[2] > $m[4]){
			#$F_2 = 1;
			#}
		}
	}
	if(exists $region{$m[3]}){
		for $i (0 .. $#{$region{$m[3]}}){
			@s = split(/\s+/, $region{$m[3]}[$i]);
			if($s[1] < $m[4] && $s[2] > $m[4]){
				$F_2 = 1;
			}
		}
	}
	if($F_1 == 1 && $F_2 == 1){
		$K = 0;
		if(exists $GAP{$m[0]}){
			for $i (0 .. $#{$GAP{$m[0]}}){
				@s = split(/\s+/, $GAP{$m[0]}[$i]);
				#	print $s[1],"\t",$m[1],"\t", $s[2],"\n";
				if(abs($s[2] - $s[1]) > 10000){
				$s[1] = $s[1] - 5000;# more region ban
				$s[2] = $s[2] + 5000;
			}
				if($s[1] < $m[1] && $s[2] > $m[1]){
					$K = 1;
					last;
				}
			}
		}
		if($K == 0){
			if(exists $GAP{$m[3]}){
				for $i (0 .. $#{$GAP{$m[3]}}){
					@s = split(/\s+/, $GAP{$m[3]}[$i]);
					if(abs($s[2] - $s[1]) > 10000){
					$s[1] = $s[1] - 5000;
					$s[2] = $s[2] + 5000;
				}
					if($s[1] < $m[4] && $s[2] > $m[4]){
						$K = 1;
						last;
					}
				}
			}
		}
		if($K == 0){
			print $_,"\n";
		}
	}
}
