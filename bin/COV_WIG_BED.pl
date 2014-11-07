#!/usr/bin/perl
#
open(I,"<",shift@ARGV);
$F = 0;
$CC = 0;
while(<I>){
	chomp;
	@m = split(/\s+/);
	if($F == 0){
		if(substr($m[0],0,1) eq "c"){
			$CC = 1;
		}
		$F = 1;
	}
	push @{$hash_b{$m[0]}}, $m[1];
	push @{$hash_e{$m[0]}}, $m[2];
	push @{$hash_k{$m[0]}}, $m[3];
}
$sum = 0;
$chr = "";
$FLL = 0;
open(I,"<",shift@ARGV);
while(<I>){
	chomp;

	if(length($_) > 20){
		$temp_chr  = $chr;
		@m = split(/\s+/);
		$chr = substr($m[1], 9);
		if($CC == 1){
			$chr = "chr".$chr;
		}
		if(!exists $HASH{$chr}){
			$OVFLAG  = 0;
			if($sum != 0){
				$COUNT{$temp_chr}{$i} = $sum;
			}
			$i = 0;
			$sum = 0;
			if($temp_chr ne ""){
				for $ii (0 .. $#{$hash_b{$temp_chr}}){
					if(!exists $COUNT{$temp_chr}){
						$N = 0;
					}
					elsif(!exists $COUNT{$temp_chr}{$ii}){
						$N = 0;
					}
					else{
						$N = $COUNT{$temp_chr}{$ii};
					}
					print $temp_chr,"\t",$hash_b{$temp_chr}[$ii],"\t",$hash_e{$temp_chr}[$ii],"\t",$N/($hash_e{$temp_chr}[$ii] - $hash_b{$temp_chr}[$ii]),"\n";
				}
			}
		}
		$HASH{$chr} = 1;
		$st = substr($m[2], 6);
		$n = -1;
		if(exists $hash_b{$chr}){
			$FLL = 1;
		}
		else{
			$FLL = 0;
		}
		#$sum = 0;
		#$OVFLAG  = 0;
		next;
	}
	if($OVFLAG == 1){
		next;
	}
	$n++;
	if($FLL == 1){
		$pos = $st + $n;
		while($hash_e{$chr}[$i] <= $pos){
			if($sum != 0){
				$COUNT{$chr}{$i} = $sum;
				$sum = 0;
			}
			$i++;
			if($i > $#{$hash_e{$chr}}){
				$OVFLAG = 1;
				last;
			}
		}
		if($OVFLAG == 1){
			next;
		}
		if($hash_b{$chr}[$i] <= $pos && $pos < $hash_e{$chr}[$i]){
			$sum+=$_;
		}
		elsif($hash_b{$chr}[$i] > $pos){
			next;
		}
	}
}
if($sum != 0){
	$COUNT{$chr}{$i} = $sum;
}

for $i (0 .. $#{$hash_b{$chr}}){
	if(!exists $COUNT{$chr}){
		$N = 0;
	}
	elsif(!exists $COUNT{$chr}{$i}){
		$N = 0;
	}
	else{
		$N = $COUNT{$chr}{$i};
	}
	print $chr,"\t",$hash_b{$chr}[$i],"\t",$hash_e{$chr}[$i],"\t",$N/($hash_e{$chr}[$i] - $hash_b{$chr}[$i]),"\n";
}



