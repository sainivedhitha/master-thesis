function final_array = first_elements(bigarray)

   ncells = length(bigarray);

    for n27 = 1 : ncells
        array1 = bigarray{1,n27};
        
        array1_len = length(array1);
        temp_array = array1(1,1);

        for i27=1:array1_len-1
            difference = - array1(1,i27) + array1(1,i27+1);
            if difference > 1000
               temp_array = [temp_array array1(1,i27+1)];
            end
        end
        
        final_array{1,n27} = temp_array;

    end
end