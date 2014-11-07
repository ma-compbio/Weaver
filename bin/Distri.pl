#!/usr/bin/perl
#
use Parallel::ForkManager;
$pm=new Parallel::ForkManager(7);
$FILE = shift@ARGV;
open(I,"<",$FILE);
while(<I>){
	chomp;
	@m = split(/\s+/);
	if(!exists $hash{$m[0]}){
		push @ALL, $m[0];
		$hash{$m[0]} = 1;
	}

}

for $i (0 .. $#ALL){
	my $pid = $pm->start and next;
	$sb = $ALL[$i];
	open(I,"<",$FILE);
	open(O,">","$sb/$sb.file");
	while(<I>){
		chomp;
		@s = split(/\s+/);
		if($s[0] eq $sb){
			print O $_,"\n";
		}
	}
	#system("cd $sb; coverageBed -split -hist -abam $sb.bam -b $sb.file | /home/yangli9/LUNA/LARGE_RUN/cal.pl | sort -k 1,1 -k 2,2n > $sb.read");
	system("cd $sb; coverageBed -split -hist -abam $sb.Y.bam -b $sb.file | /home/yangli9/LUNA/LARGE_RUN/cal.pl | sort -k 1,1 -k 2,2n > $sb.read");
	$pm->finish;
}
$pm->wait_all_children;
system("rm tempread");
for $i (0 .. $#ALL){
	$sb = $ALL[$i];
	print $sb,"\n";
	system("cat $sb/$sb.read >> tempread");
}

