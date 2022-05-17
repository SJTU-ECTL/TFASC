#!/usr/bin/perl

my $main_output_dir = $ARGV[0]; # like "./output_dir/demo_results"
my $calling_number_filename_in = $ARGV[1];

my $calling_number_summary_filename = "${main_output_dir}/$calling_number_filename_in";

open(outfile, ">$calling_number_summary_filename") or die "Cannot open output file $calling_number_summary_filename\n";


#for(my $group_id = 1; $group_id <= 4; $group_id++){
for(my $group_id = 1; $group_id <= 1; $group_id++){
	for(my $bm_id = 1; $bm_id <= 12; $bm_id++){
	
		my $bm_name = "bm${bm_id}.${group_id}";
		if($bm_name eq "bm7.1" or $bm_name eq "bm7.2"){
			next;
		}
		
		my $output_subdir = "${main_output_dir}/${bm_name}";
		
		my $report_filename = "${output_subdir}/ZhaoMethod_LUT_calling_number.txt";

		open(infile, "<$report_filename") or die "Can't open ZhaoMethod_LUT_calling_number.txt for ${bm_name}\n";

		my $zhao_num, $LUT_num;
		while(<infile>){
			chomp;
			if(/Zhao method calling number = (\d+)/){
				print "$_\n";
				$zhao_num = $1;
				print "$zhao_num\n";
				print outfile "$bm_name zhao_number = $zhao_num, ";
			}
			if(/LUT entry found number = (\d+)/){
				print "$_\n";
				$LUT_num = $1;
				print "$LUT_num\n";
				print outfile "LUT entry found number = $LUT_num\n";
			}
		}
		close(infile);
	}
}

close(outfile);







