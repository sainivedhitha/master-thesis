function sdl_mtx = Resystem1(sig_sdl);

%global sdl_mtx;

temp = sig_sdl;
Row = 12;
Col = 5;

sdl_mtx = cell(Row,Col);

k = 1;
for column = 2 : Col - 1
    sdl_mtx{1,column} = temp(4-k, :);
    %disp(['row: ' num2str(1) ' column: ' num2str(column) ' index: ' num2str(4-k)]);
    k = k + 1;
end

for row = 2 : Row - 1
   k = 0;
   if mod(row,2) == 0      						 % siamo sulla riga pari
      for column = 1 : Col
         sdl_mtx{row,column} = temp((Col - 1) + (row - 2)*Col + k, :);
         %disp(['row: ' num2str(row) ' column: ' num2str(column) ' index: ' num2str((Col - 1) + (row - 2)*Col + k)]);
         k = k + 1;
      end
   elseif mod(row,2) == 1                     % siamo sulla riga dispari
      for column = 1 : Col
         sdl_mtx{row,column} = temp((Col - 2) + (row - 1)*Col - k, :);
          %disp(['row: ' num2str(row) ' column: ' num2str(column) ' index: ' num2str((Col - 2) + (row - 1)*Col - k)]);
         k = k + 1;
      end 
   end
end

k = 0;
for column = 2 : Col - 1
   sdl_mtx{12,column} = temp(54 + k, :);
   %disp(['row: ' num2str(12) ' column: ' num2str(column) ' index: ' num2str(54+k)]);
   k = k + 1;
end
