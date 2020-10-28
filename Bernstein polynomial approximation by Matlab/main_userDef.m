clc;
clear all;

%=======================================================================
% here are the user-input parameters: n for degree, and m for precision
n = 6;
m = 8;
%=======================================================================

current_dir = '.\\';
addpath([current_dir, 'common_functions\\']);

demo_output_summary_filename = [current_dir, 'output_dir\\user_results\\summary.txt'];
demo_output_for_TFASC_input_filename = [current_dir, 'output_dir\\user_results\\input_files_for_TFASC\\TFASC_input_file.txt'];

fid_summary = fopen(demo_output_summary_filename,'w');
fid_TFASC = fopen(demo_output_for_TFASC_input_filename, 'w');

degree = n;
[coefs, obj, H, b, c] = Bern_appr(@target_function_user_defined, degree);
		
fprintf(fid_summary, '%s%d%s%d\r\n', 'n = ', n, ', m = ', m);
Bernstein_coefficient_vec = coefs'
fprintf(fid_summary, '%s ', 'Bernstein coefficients: ');
fid_summary = print_double_vec_to_file(fid_summary, Bernstein_coefficient_vec);
fprintf(fid_summary, '\r\n');
		
feature_vector = BernCoefVec_to_featureVector_converter(Bernstein_coefficient_vec, m)
fprintf(fid_summary, '%s ', 'feature vector: ');
fid_summary = print_int_vec_to_file(fid_summary, feature_vector);
fprintf(fid_summary, '\r\n\r\n');
        
% print to generated input file for TFASC
fprintf(fid_TFASC, '%d\r\n%d\r\n', n, m);
fid_TFASC = print_int_vec_to_file(fid_TFASC, feature_vector);
        
fclose(fid_TFASC);	
fclose(fid_summary);



