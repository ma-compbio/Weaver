#!/usr/bin/perl
#
$CHR_ = shift@ARGV;
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
	if($CHR_ eq $m[0]){
		push @{$hash_b{$m[0]}}, $m[1];
		push @{$hash_e{$m[0]}}, $m[2];
		push @{$hash_k{$m[0]}}, $m[3];
	}
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
			if($temp_chr eq $CHR_){
				last;
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
if($sum != 0 && $chr eq $CHR_){
	$COUNT{$chr}{$i} = $sum;
}

for $i (0 .. $#{$hash_b{$CHR_}}){
	if(!exists $COUNT{$CHR_}){
		$N = 0;
	}
	elsif(!exists $COUNT{$CHR_}{$i}){
		$N = 0;
	}
	else{
		$N = $COUNT{$CHR_}{$i};
	}
	print $CHR_,"\t",$hash_b{$CHR_}[$i],"\t",$hash_e{$CHR_}[$i],"\t",$N/($hash_e{$CHR_}[$i] - $hash_b{$CHR_}[$i]),"\n";
}



