// Benchmark "/home/wch/wangchen_research/approxSC_upload/DA_programs/linux_program/output_dir/demo_results/bm1.2/bm1.2-solution-001" written by ABC on Thu Oct 29 11:11:05 2020

module \/home/wch/wangchen_research/approxSC_upload/DA_programs/linux_program/output_dir/demo_results/bm1.2/bm1.2-solution-001  ( 
    x00, x01, x02, x03, x04, x05, x06, x07, x08, x09, x10, x11,
    z0  );
  input  x00, x01, x02, x03, x04, x05, x06, x07, x08, x09, x10, x11;
  output z0;
  wire new_n14_, new_n15_, new_n16_, new_n17_, new_n18_, new_n19_, new_n20_,
    new_n21_, new_n22_, new_n23_, new_n24_, new_n25_, new_n26_, new_n27_,
    new_n28_, new_n29_, new_n30_, new_n31_, new_n32_, new_n33_, new_n34_,
    new_n35_, new_n36_;
  inv1   g00(.a(x08), .O(new_n14_));
  inv1   g01(.a(x09), .O(new_n15_));
  inv1   g02(.a(x10), .O(new_n16_));
  nand3  g03(.a(new_n16_), .b(new_n15_), .c(new_n14_), .O(new_n17_));
  inv1   g04(.a(x07), .O(new_n18_));
  inv1   g05(.a(x11), .O(new_n19_));
  aoi21  g06(.a(new_n19_), .b(x02), .c(new_n18_), .O(new_n20_));
  oai21  g07(.a(new_n20_), .b(new_n17_), .c(x00), .O(new_n21_));
  nand3  g08(.a(x05), .b(x04), .c(x03), .O(new_n22_));
  inv1   g09(.a(new_n22_), .O(new_n23_));
  nand2  g10(.a(new_n23_), .b(x06), .O(new_n24_));
  inv1   g11(.a(new_n24_), .O(new_n25_));
  nand2  g12(.a(new_n25_), .b(new_n21_), .O(new_n26_));
  nand2  g13(.a(new_n18_), .b(x02), .O(new_n27_));
  nand4  g14(.a(new_n15_), .b(new_n14_), .c(x02), .d(x01), .O(new_n28_));
  nand3  g15(.a(new_n28_), .b(new_n27_), .c(x00), .O(new_n29_));
  nor2   g16(.a(new_n22_), .b(x06), .O(new_n30_));
  inv1   g17(.a(x04), .O(new_n31_));
  inv1   g18(.a(x05), .O(new_n32_));
  oai21  g19(.a(new_n32_), .b(new_n31_), .c(x03), .O(new_n33_));
  nand4  g20(.a(new_n32_), .b(new_n31_), .c(x02), .d(x01), .O(new_n34_));
  nand2  g21(.a(new_n34_), .b(new_n33_), .O(new_n35_));
  aoi21  g22(.a(new_n30_), .b(new_n29_), .c(new_n35_), .O(new_n36_));
  nand2  g23(.a(new_n36_), .b(new_n26_), .O(z0));
endmodule


