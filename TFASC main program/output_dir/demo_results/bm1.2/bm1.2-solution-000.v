// Benchmark "/home/wch/wangchen_research/approxSC_upload/DA_programs/linux_program/output_dir/demo_results/bm1.2/bm1.2-solution-000" written by ABC on Thu Oct 29 11:11:05 2020

module \/home/wch/wangchen_research/approxSC_upload/DA_programs/linux_program/output_dir/demo_results/bm1.2/bm1.2-solution-000  ( 
    x00, x01, x02, x03, x04, x05, x06, x07, x08, x09, x10, x11,
    z0  );
  input  x00, x01, x02, x03, x04, x05, x06, x07, x08, x09, x10, x11;
  output z0;
  wire new_n14_, new_n15_, new_n16_, new_n17_, new_n18_, new_n19_, new_n20_,
    new_n21_, new_n22_, new_n23_, new_n24_, new_n25_, new_n26_, new_n27_,
    new_n28_, new_n29_, new_n30_, new_n31_, new_n32_, new_n33_;
  inv1   g00(.a(x01), .O(new_n14_));
  inv1   g01(.a(x04), .O(new_n15_));
  inv1   g02(.a(x05), .O(new_n16_));
  nand2  g03(.a(new_n16_), .b(new_n15_), .O(new_n17_));
  inv1   g04(.a(x06), .O(new_n18_));
  inv1   g05(.a(x08), .O(new_n19_));
  inv1   g06(.a(x09), .O(new_n20_));
  nand4  g07(.a(new_n20_), .b(new_n19_), .c(new_n18_), .d(x03), .O(new_n21_));
  aoi21  g08(.a(new_n21_), .b(new_n17_), .c(new_n14_), .O(new_n22_));
  inv1   g09(.a(x07), .O(new_n23_));
  nand3  g10(.a(new_n23_), .b(new_n18_), .c(x03), .O(new_n24_));
  inv1   g11(.a(new_n24_), .O(new_n25_));
  oai21  g12(.a(new_n25_), .b(new_n22_), .c(x02), .O(new_n26_));
  inv1   g13(.a(x10), .O(new_n27_));
  nand2  g14(.a(new_n27_), .b(x06), .O(new_n28_));
  nand2  g15(.a(x11), .b(x07), .O(new_n29_));
  nand3  g16(.a(new_n29_), .b(new_n20_), .c(new_n19_), .O(new_n30_));
  nor2   g17(.a(new_n30_), .b(new_n28_), .O(new_n31_));
  nand3  g18(.a(x05), .b(x04), .c(x00), .O(new_n32_));
  oai21  g19(.a(new_n32_), .b(new_n31_), .c(x03), .O(new_n33_));
  nand2  g20(.a(new_n33_), .b(new_n26_), .O(z0));
endmodule


