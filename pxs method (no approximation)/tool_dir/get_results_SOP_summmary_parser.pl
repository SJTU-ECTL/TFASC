#!/usr/bin/perl

my $main_output_dir = $ARGV[0]; # like "./output_dir/demo_results"
my $hardware_cost_summary_filename_in = $ARGV[1];

my $hardware_cost_summary_filename = "${main_output_dir}/$hardware_cost_summary_filename_in";

open(outfile, ">$hardware_cost_summary_filename") or die "Cannot open output file $hardware_cost_summary_filename\n";

my $bm_count = 0;
my $runtime_in_seconds_sum = 0;


for(my $group_id = 1; $group_id <= 4; $group_id++){
	for(my $bm_id = 1; $bm_id <= 12; $bm_id++){
	
		my $bm_name = "bm${bm_id}.${group_id}";
		#if($bm_name eq "bm7.1" or $bm_name eq "bm7.2"){
		#	next;
		#}
		
		$bm_count++;

		my $output_subdir = "${main_output_dir}/${bm_name}";
		
		my $input_hardware_cost_report_filename = "${output_subdir}/bestSol_summary_detailed_ESPRESSO.txt";

		open(infile, "<$input_hardware_cost_report_filename") or die "Can't open input_hardware_cost_report for ${bm_name}\n";

		my $literalCount, $cubeCount;
		while(<infile>){
			chomp;
			if(/SOP literal count \(Espresso\) = (\d+)/){
				print "$_\n";
				$literalCount = $1;
				print "$literalCount\n";
				print outfile "$bm_name $literalCount ";
			}
			if(/cube count \(Espresso\) = (\d+)/){
				print "$_\n";
				$cubeCount = $1;
				print "$cubeCount\n";
				print outfile "$cubeCount\n";
			}
		}
		close(infile);
	}
}

close(outfile);







