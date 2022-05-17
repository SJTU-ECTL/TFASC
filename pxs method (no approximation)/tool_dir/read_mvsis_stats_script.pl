#!/usr/bin/perl

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

my $lit_sop_count = 0;
my $lit_fac_count = 0;
my $cube_count = 0;

open(infile, "<$infile");

while(<infile>){
	if(/lits\s+\=\s+(\d+)\s+/){
		$lit_sop_count = $1;
	}
	if(/lits\(ff\)\s+\=\s+(\d+)\s+/){
		$lit_fac_count = $1;
	}
	if(/cubes\s+\=\s+(\d+)\s+/){
		$cube_count = $1;
	}
}	
close(infile);

#print "lit_sop_count = $lit_sop_count\n";
#print "lit_fac_count = $lit_fac_count\n";
#print "cube_count = $cube_count\n";

open(outfile, ">$outfile");
print outfile "${lit_sop_count} ${lit_fac_count} ${cube_count}";
close(outfile);









