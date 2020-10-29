// Benchmark "/home/wch/wangchen_research/approxSC_upload/DA_programs/linux_program/output_dir/user_results/test1/test1-solution-000" written by ABC on Thu Oct 29 11:10:53 2020

module \/home/wch/wangchen_research/approxSC_upload/DA_programs/linux_program/output_dir/user_results/test1/test1-solution-000  ( 
    x0, x1, x2, x3, x4, x5, x6, x7,
    z0  );
  input  x0, x1, x2, x3, x4, x5, x6, x7;
  output z0;
  wire new_n10_;
  nand4  g0(.a(x6), .b(x5), .c(x4), .d(x0), .O(new_n10_));
  and2   g1(.a(new_n10_), .b(x3), .O(z0));
endmodule


