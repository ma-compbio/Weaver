#!/usr/bin/perl
#
use Parallel::ForkManager;
use FindBin qw($Bin);
$FILE = shift@ARGV;
$WIG = shift@ARGV;
$TH = shift@ARGV;
$pm=new Parallel::ForkManager($TH);
open(I,"<",$FILE);
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
	$CHR{$m[0]} = 1;
}
for $chr (keys %CHR){
	my $pid = $pm->start and next;
	print $chr,"\n";
	system("$Bin/SINGLE_COV_WIG_BED.pl $chr $FILE $WIG > tempread_$chr");
	$pm->finish;
}
$pm->wait_all_children;
system("cat tempread_* > tempread; rm tempread_*");
