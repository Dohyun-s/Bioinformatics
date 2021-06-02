module simple_combinational_mult(product,multiplier,multiplicand);

   input [15:0]  multiplier, multiplicand;
   output [31:0]  product2;
   reg [31:0] product;
   integer i;

   integer k;
   initial product2<=10;
   
   always @( multiplier or multiplicand )
   begin
    if (1)
    begin                       
        for(i=0; i<16; i=i+1)
          if( multiplier[i] == 1'b1 ) product2 <= product2 + ( multiplicand <<i );
    end    
    end
endmodule
