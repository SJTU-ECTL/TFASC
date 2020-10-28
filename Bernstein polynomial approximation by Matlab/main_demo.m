clc;
close all;
clear all;

current_dir = '.\\';
addpath([current_dir, 'demo_target_functions\\']);
addpath([current_dir, 'common_functions\\']);

demo_output_summary_filename = [current_dir, 'output_dir\\demo_results\\summary.txt'];
demo_output_for_TFASC_dir = [current_dir, 'output_dir\\demo_results\\input_files_for_TFASC\\'];

demo_output_summary_filename
fid_summary = fopen(demo_output_summary_filename,'w');

bm_id_start = 1;
bm_id_end = 12;

for test_group_id = 1:4
	if(test_group_id == 1)
		n = 4;
		m = 4;
	elseif(test_group_id == 2)
		n = 4;
		m = 8;
	elseif(test_group_id == 3)
		n = 6;
		m = 4;
    elseif(test_group_id == 4)
		n = 6;
		m = 8;
	end
	
	for bm_id = bm_id_start:1:bm_id_end
    
        TFASC_input_filename = [demo_output_for_TFASC_dir,'bm',num2str(bm_id),'.',num2str(test_group_id),'.txt'];
        fid_TFASC = fopen(TFASC_input_filename,'w');
                
        bm_name = ['bm',num2str(bm_id),'.',num2str(test_group_id)]
        
		if(bm_id == 1)
			target_func_test = @func1;
			fprintf(fid_summary, '%s\r\n', [bm_name,': sin(x)']);
		elseif(bm_id == 2)
			target_func_test = @func2;
			fprintf(fid_summary, '%s\r\n', [bm_name,': cos(x)']);
		elseif(bm_id == 3)
			target_func_test = @func3;
			fprintf(fid_summary, '%s\r\n', [bm_name,': exp(-x)']);
		elseif(bm_id == 4)
			target_func_test = @func4;
			fprintf(fid_summary, '%s\r\n', [bm_name,': log(x+1)']);
		elseif(bm_id == 5)
			target_func_test = @func5;
			fprintf(fid_summary, '%s\r\n', [bm_name,': sin(pi*x)/pi']);
		elseif(bm_id == 6)
			target_func_test = @func6;
			fprintf(fid_summary, '%s\r\n', [bm_name,': tanh(x)']);
		elseif(bm_id == 7)
			target_func_test = @func7;
			fprintf(fid_summary, '%s\r\n', [bm_name,': tanh(4x)']);
		elseif(bm_id == 8)
			target_func_test = @func8;
			fprintf(fid_summary, '%s\r\n', [bm_name,': x^0.45']);
		elseif(bm_id == 9)
			target_func_test = @func9;
			fprintf(fid_summary, '%s\r\n', [bm_name,': exp(-2x)']);
		elseif(bm_id == 10)
			target_func_test = @func10;
			fprintf(fid_summary, '%s\r\n', [bm_name,': sigmoid(x)']);
		elseif(bm_id == 11)
			target_func_test = @func11;
			fprintf(fid_summary, '%s\r\n', [bm_name,': x^2.2']);
		elseif(bm_id == 12)
			target_func_test = @func12;
			fprintf(fid_summary, '%s\r\n', [bm_name,': 0.5*cos(pi*x)+0.5']);
		end
		
		degree = n;
		[coefs, obj, H, b, c] = Bern_appr(target_func_test, degree);
		
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
    end    
end

fclose(fid_summary);


