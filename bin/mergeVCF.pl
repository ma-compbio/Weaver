#!/usr/bin/perl
#
#open(I,"<",shift@ARGV);
$flag = 0;
while(<STDIN>){
	chomp;
	if(substr($_,0,1) eq "#"){
		if($flag == 0){
			print $_,"\n";
		}
		next;
	}
	else{
		$flag = 1;
	}
	print $_,"\n";
	#@m = split(/\s+/);
}

