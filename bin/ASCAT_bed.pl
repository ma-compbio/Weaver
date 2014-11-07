#!/usr/bin/perl
#430625  1475054 TCGA-13-0727-01 8       3               331     19:24380818..28974637
#open(I,"<",shift@ARGV);
while(<STDIN>){
	chomp;
	@m = split(/\s+/);
	@t = split(/:|\./,$m[$#m]);
	if( $t[0] == 23){
		 $t[0] = "X";
	}
	print $t[0],"\t",$t[1],"\t",$t[3],"\t",$m[3]-$m[4],"\t",$m[4],"\n";
}

