#!/usr/bin/perl

my $dir = $ARGV[0];

if(-e $dir){
	`rm $dir*`;	
}
else{
	`mkdir $dir`;
}








