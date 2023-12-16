module dtage #(parameter NUM_TABLES = 2, parameter NUM_TILES = 4)(
	input clk,
	input rst_n,
	input [15:0] pc,
	input [15:0] h [0:NUM_TABLES-1],
	output logic [15:0] prediction
);
	// 2 HISTORY TABLES, each can have a max of 4 tiles
	logic [2:0] config_vector [2];
	assign config_vector[0] = 1;
	assign config_vector[1] = 1;
	
	logic  history_vector[2];
	assign history_vector[0] = 0;
	assign history_vector[1] = 0;
	// pre
	genvar i;
	wire [15:0] base_prediction;
	assign base_prediction = '0;
	logic [2:0] [0:NUM_TABLES-1] scatter_index;
	logic [0:NUM_TABLES-1] tag;
	// inbetween
	typedef struct packed{
		reg [15:0] prediction; 
		reg [5:0] tag; 
	} line;
	typedef struct packed{
		reg [15:0] input_index;
		line line0;
		line line1;
		reg buff_sel;
	} tile;
	
	// possible parametrization?
	tile tiles[4];

	// post
	reg [15:0] tile0_pred_out;
	reg [15:0] tile1_pred_out;
	reg [15:0] tile2_pred_out;
	reg [15:0] tile3_pred_out;
	reg [5:0] tile0_tag_out;
	reg [5:0] tile1_tag_out;
	reg [5:0] tile2_tag_out;
	reg [5:0] tile3_tag_out;
	
	reg [15:0] table0_pred_out;
	reg [15:0] table1_pred_out;
	reg [5:0] table0_tag_out;
	reg [5:0] table1_tag_out;
	
	reg [15:0] mux0, mux1;
	
	// data flow
	for(i = 0; i < NUM_TABLES; i++) begin
		assign scatter_index[i] = pc[0:2] ^ h[i][0:2];
		assign tag[i] = pc ^ h[i];
	end

	always @ (*) begin 
		tiles[0].input_index = scatter_index[0];
			
		if (scatter_index[0][2:1] == 0) begin
			tiles[0].buff_sel = 1;
		end else
			tiles[0].buff_sel = 0;
			
		tiles[1].input_index = scatter_index[1];
			
		if (scatter_index[1][2:1] == 1) begin
			tiles[1].buff_sel = 1;
		end else
			 tiles[1].buff_sel = 0;
			
		tiles[2].input_index = scatter_index[2];
			
		if (scatter_index[2][2:1] == 2) begin
			tiles[2].buff_sel = 1;
		end else
			tiles[2].buff_sel = 0;
			
		tiles[3].input_index = scatter_index[3];
			
		if (scatter_index[3][2:1] == 3) begin
			tiles[3].buff_sel = 1;
		end else
			tiles[3].buff_sel = 0;
			
		tile0_pred_out = (tiles[0].buff_sel)? {(scatter_index[0][0])? tiles[0].line1.prediction: tiles[0].line0.prediction}: 'x;
		tile1_pred_out = (tiles[1].buff_sel)? {(scatter_index[1][0])? tiles[1].line1.prediction: tiles[1].line0.prediction}: 'x;
		tile2_pred_out = (tiles[2].buff_sel)? {(scatter_index[2][0])? tiles[2].line1.prediction: tiles[2].line0.prediction}: 'x;
		tile3_pred_out = (tiles[3].buff_sel)? {(scatter_index[3][0])? tiles[3].line1.prediction: tiles[3].line0.prediction}: 'x;
		
		tile0_tag_out = (tiles[0].buff_sel)? {(scatter_index[0][0])? tiles[0].line1.tag: tiles[0].line0.tag}: 'x;
		tile1_tag_out = (tiles[1].buff_sel)? {(scatter_index[1][0])? tiles[1].line1.tag: tiles[1].line0.tag}: 'x;
		tile2_tag_out = (tiles[2].buff_sel)? {(scatter_index[2][0])? tiles[2].line1.tag: tiles[2].line0.tag}: 'x;
		tile3_tag_out = (tiles[3].buff_sel)? {(scatter_index[3][0])? tiles[3].line1.tag: tiles[3].line0.tag}: 'x;
		
		table0_pred_out = (history_vector[0])? tile1_pred_out: tile0_pred_out;
		table0_tag_out = (history_vector[0])? tile1_tag_out: tile0_tag_out;
		table1_pred_out = (history_vector[1])? tile3_pred_out: tile2_pred_out;
		table1_tag_out = (history_vector[1])? tile3_tag_out: tile2_tag_out;
		
		assign mux0 = (tag[0] == table0_tag_out)? table0_pred_out: base_prediction;
		assign mux1 = (tag[1] == table1_tag_out)? table1_pred_out: mux0;
		assign prediction = mux1;
	end
		
	
	
endmodule