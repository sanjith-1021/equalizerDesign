function read_order = int_indx_gen(depth)
    row_indices = zeros(40*depth, 1);
    col_indices = zeros(40*depth, 1);
    row = 1; col = 1;
    last_col_when_row_max = 1;
    for i = 1:40*depth
        row_indices(i) = row;
        col_indices(i) = col;
        if row == 40
            row = 1;
            col = last_col_when_row_max + 1;
            last_col_when_row_max = last_col_when_row_max + 1;
        else
            row = row + 1;
            col = mod(col - 18, depth) + 1;
        end
    end
    read_order = 40*(col_indices - 1) + (row_indices - 1) + 1;
end