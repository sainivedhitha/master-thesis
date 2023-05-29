function dd_mtx = filt_spat2(sdl_mtx);

% DD signal
Row = 11;
Col = 5;

dd_mtx = cell(11,5);

for column = 2 : Col - 1
   dd_mtx{1,column} = sdl_mtx{1,column} - sdl_mtx{2,column};
end

for row = 2 : Row - 1
   for column = 1 : Col
      dd_mtx{row,column} = sdl_mtx{row,column} - sdl_mtx{row + 1,column};
   end
end

for column = 2 : Col - 1
   dd_mtx{11,column} = sdl_mtx{11,column} - sdl_mtx{12,column};
end