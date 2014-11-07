#!/usr/bin/perl
#
# From SV list to connected regions

open(I,"<",shift@ARGV);
while(<I>){
	chomp;
	@m = split(/\s+/);
	push @{$hash{$m[0]}},$m[1];
	push @{$hash{$m[3]}},$m[4];
	$LINK{$m[0]}{$m[1]} = [$m[3],$m[4],$m[2],$m[5]];
	$LINK{$m[3]}{$m[4]} = [$m[0],$m[1],$m[5],$m[2]];
}

foreach $key (keys %hash){
	@A  =  sort { $a <=> $b } @{$hash{$key}};
	@{$hash{$key}} = @A;
}

%USED;

$N=0;

$RANGE_LIMIT = 5000;
$FRAGMENT_LIMIT = 20;
#grow("11",111028870);
grow("3",146385000);
#if()
#print join("\t",@{$hash{"chr17"}}),"\n";
sub grow{
	if($N > $FRAGMENT_LIMIT){
		return;
	}
	my $CHR = shift @_;
	my $B = shift @_;
	die "chr $CHR does not exist in SV\n" unless exists $LINK{$CHR};
	if(!exists $LINK{$CHR}{$B}){
		die "$B > $hash{$CHR}[$#{$hash{$CHR}}]\n" if ($B > $hash{$CHR}[$#{$hash{$CHR}}]);
		for $i (0 .. $#{$hash{$CHR}}){
			if($hash{$CHR}[$i] > $B){
				if($RANGE_LIMIT < abs($B - $hash{$CHR}[$i])){
					die "$B - $hash{$CHR}[$#{$hash{$CHR}}] > $RANGE_LIMIT\n";
				}
				print $CHR, "\t", $B ,"\t" , $hash{$CHR}[$i],"\t","+\tSTART\n";
				grow($CHR,$hash{$CHR}[$i]);
				return;
			}
		}
	}
	$USED{$CHR}{$B}=1;
	#print $CHR,"\t",$B,"\n";
	my $next_chr = $LINK{$CHR}{$B}[0];
	my $next_b = $LINK{$CHR}{$B}[1];
	my $next_ori = $LINK{$CHR}{$B}[3];
	$USED{$next_chr}{$next_b}=1;
	#print $next_chr,"\t",$next_b,"\t",$next_ori,"\n";
	my $FIND_FLAG = 0;
	my $RANK = -1;
	for my $i (0 .. $#{$hash{$next_chr}}){
		if($next_b == $hash{$next_chr}[$i]){
			$RANK = $i;
			last;
		}
	}
	if($next_ori eq "-"){
		$RANK++;
		while($RANK <= $#{$hash{$next_chr}} && $RANGE_LIMIT > abs($hash{$next_chr}[$RANK] - $next_b)){
			if(!exists $USED{$next_chr}{$hash{$next_chr}[$RANK]}){
				print $next_chr,"\t",$next_b,"\t",$hash{$next_chr}[$RANK],"\t+\n";
				$FIND_FLAG = 1;
				$N++;
				grow($next_chr, $hash{$next_chr}[$RANK]);
				last;
			}
			$RANK++;
		}
	}
	else{
		$RANK--;
		while($RANK >= 0 && $RANGE_LIMIT > abs($hash{$next_chr}[$RANK] - $next_b)){
			if(!exists $USED{$next_chr}{$hash{$next_chr}[$RANK]}){
				print $next_chr,"\t",$hash{$next_chr}[$RANK],"\t",$next_b,"\t-\n";
				$FIND_FLAG = 1;
				$N++;
				grow($next_chr, $hash{$next_chr}[$RANK]);
				last;
			}
			$RANK--;
		}
	}
	if($FIND_FLAG == 0){
		if($next_ori eq "-"){
			print $next_chr,"\t",$next_b,"\t",$next_b+10,"\t+\tEND\n";
		}
		else{
			print $next_chr,"\t",$next_b-10,"\t",$next_b,"\t-\tEND\n";
		}
	}	
}
=cut
foreach $key (keys %USED){
	foreach $B (keys %{$USED{$key}}){
		print $key,"\t",$B,"\n";
	}
}
