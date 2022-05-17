#!/usr/bin/perl

my $main_output_dir = $ARGV[0]; # like "./output_dir/demo_results"
my $hardware_cost_summary_filename_in = $ARGV[1];
my $runtime_summary_filename_in = $ARGV[2];

my $hardware_cost_summary_filename = "${main_output_dir}/$hardware_cost_summary_filename_in";
my $runtime_summary_filename = "${main_output_dir}/$runtime_summary_filename_in";

open(outfile1, ">$runtime_summary_filename") or die "Cannot open output file $runtime_summary_filename\n";
open(outfile2, ">$hardware_cost_summary_filename") or die "Cannot open output file $hardware_cost_summary_filename\n";

my $bm_count = 0;
my $runtime_in_seconds_sum = 0;

#for(my $group_id = 1; $group_id <= 4; $group_id++){
for(my $group_id = 1; $group_id <= 1; $group_id++){
	for(my $bm_id = 1; $bm_id <= 12; $bm_id++){
	
		my $bm_name = "bm${bm_id}.${group_id}";
		if($bm_name eq "bm7.1" or $bm_name eq "bm7.2"){
			next;
		}
		
		$bm_count++;

		my $output_subdir = "${main_output_dir}/${bm_name}";
		
		my $input_runtime_report_filename = "${output_subdir}/runtime_report.txt";
		my $input_hardware_cost_report_filename = "${output_subdir}/bestSol_summary_detailed_ABC.txt";

		#open(FILE, "<./runtime_report.txt") or die "Can't open runtime_report.txt\n"; 
		open(infile1, "<$input_runtime_report_filename") or die "Can't open input_runtime_report for ${bm_name}\n";

		my $runtime_in_seconds;
		while(<infile1>){
			chomp;
			if(/chrono.*\s=\s(\S+)/){
				print "$_\n";
				$runtime_in_seconds = $1;
				print "$runtime_in_seconds\n";
				print outfile1 "$bm_name $runtime_in_seconds\n";
				$runtime_in_seconds_sum = $runtime_in_seconds_sum + $runtime_in_seconds;
			}
		}

		close(infile1);


		open(infile2, "<$input_hardware_cost_report_filename") or die "Can't open input_hardware_cost_report for ${bm_name}\n";

		my $area, $delay, $ADP;
		while(<infile2>){
			chomp;
			if(/area \(ABC\) = \s+(\S+)/){
				print "$_\n";
				$area = $1;
				print "$area\n";
				print outfile2 "$bm_name $area ";
			}
			if(/delay \(ABC\) = \s+(\S+)/){
				print "$_\n";
				$delay = $1;
				print "$delay\n";
				print outfile2 "$delay ";
			}
			if(/area-delay-product \(ABC\) = \s+(\S+)/){
				print "$_\n";
				$ADP = $1;
				print "$ADP\n";
				print outfile2 "$ADP\n";
			}
		}

		close(infile2);
	}
}

my $runtime_in_seconds_avg = $runtime_in_seconds_sum / $bm_count;
print outfile1 "\nruntime_avg(seconds): $runtime_in_seconds_avg\n";

close(outfile1);
close(outfile2);







