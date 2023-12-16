`timescale 1ns / 1ps

module dtage_tb;

    // Inputs
    reg clk;
    reg rst_n;
    reg [15:0] pc;
    reg [15:0] h [0:3];

    // Output
    wire [15:0] prediction;

    dtage uut (
        .clk(clk), 
        .rst_n(rst_n), 
        .pc(pc), 
        .h(h), 
        .prediction(prediction)
    );

    always #10 clk = ~clk; // Clock with period of 20ns

    // Test sequence
    initial begin
        // Initialize Inputs
		
		#1;
        clk = 0;
        rst_n = 0;
        pc = 0;
        h[0] = 16'h0;
        h[1] = 16'h0;
        h[2] = 16'h0;
        h[3] = 16'h0;
		
		#1;
		rst_n = 1;
	// waiting for rest
        #100;
       

        // Apply test stimulus

        #20 pc = 16'h0001; h[0] = 16'hAAAA;
        #20 pc = 16'h0002; h[1] = 16'h5555;
        #20 pc = 16'h0003; h[2] = 16'hAAAA;
        #20 pc = 16'h0004; h[3] = 16'h5555;

        #100 $finish;
    end
    
    always @(posedge clk) begin
        if (!rst_n) begin // Only display after reset is released
            $display("Time=%t | PC=%h | Prediction=%h", $time, pc, prediction);
        end
    end

endmodule