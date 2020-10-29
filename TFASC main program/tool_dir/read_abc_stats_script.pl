#!/usr/bin/perl

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

my $area = 0;
my $delay = 0;

open(infile, "<$infile");

while(<infile>){
	if(/area\s*\=\s*(\d+\.\d+)\s+/){
		$area = $1;
	}
	if(/delay\s*\=\s*(\d+\.\d+)\s+/){
		$delay = $1;
	}
}	
close(infile);

#print "area = $area\n";
#print "delay = $delay\n";

open(outfile, ">$outfile");
print outfile "${area} ${delay}";
close(outfile);









