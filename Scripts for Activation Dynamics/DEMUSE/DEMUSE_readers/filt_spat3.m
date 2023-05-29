function lp_mtx = filt_spat3(sdl_mtx,sig);

% DD signal
Row = 11;
Col = 5;


lp_mtx = cell(11,3);

for row = 1 : Row
    if mod(row,2) == 0  
        for column = 2 : Col-1
            %disp(['row: ' num2str(row) ' column: ' num2str(column) ' index: ' num2str((Col - 1) + (row)*Col - column-1) ' ' num2str((Col - 1) + (row)*Col - column)]);
            lp_mtx{row,column-1} = sdl_mtx{row,column}  - sdl_mtx{row + 1,column};
            lp_mtx{row,column-1} = lp_mtx{row,column-1} + sig((Col - 1) + (row)*Col - column-1, :);
            lp_mtx{row,column-1} = lp_mtx{row,column-1} - sig((Col - 1) + (row)*Col - column, :);
        end
    elseif mod(row,2) == 1  
        for column = 2 : Col-1
            %disp(['row: ' num2str(row) ' column: ' num2str(column) ' index: ' num2str((Col - 2) + (row-1)*Col + column-1) ' ' num2str((Col - 2) + (row - 1)*Col + column)]);
            lp_mtx{row,column-1} = sdl_mtx{row,column}  - sdl_mtx{row + 1,column};
            lp_mtx{row,column-1} = lp_mtx{row,column-1} + sig((Col - 2) + (row - 1)*Col + column-1, :);
            lp_mtx{row,column-1} = lp_mtx{row,column-1} - sig((Col - 2) + (row - 1)*Col + column, :);
        end
    end
end