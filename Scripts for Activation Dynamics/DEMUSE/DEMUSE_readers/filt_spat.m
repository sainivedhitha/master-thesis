function [sig_sdl] = filt_spat(sig);

Row = 12;
Col = 5;

% SDL signal
for row = 1 : Row
   k = 0;
   if row == 1
      for ch = 1 : Col - 2                      % da 1 a 3  
%          ch : (6 - k)
         sig_sdl(ch,:) = sum(sig(ch : (6 - k),:));
         k = k + 1;
      end
   elseif row == 2
      for ch = Col - 1 : 2*Col - 3              % da 4 a 7
%          ch : (2*Col + 2) - k
         sig_sdl(ch,:) = sum(sig(ch : (2*Col + 2) - k,:));
         k = k + 1;
      end
      sig_sdl(8,:) = sig(8,:);
   elseif row == 3
      for ch =  2*Col - 1 : 3*Col - 3           % da 9 a 12
%          ch : (3*Col + 2) - k
         sig_sdl(ch,:) = sum(sig(ch : (3*Col + 2) - k,:));
         k = k + 1;
      end
      sig_sdl(13,:) = sig(13,:);
   elseif row == 4
      for ch = 3*Col - 1 : 4*Col - 3            % da 14 a 17
%          ch : (4*Col + 2) - k
         sig_sdl(ch,:) = sum(sig(ch : (4*Col + 2) - k,:));
         k = k + 1;
      end
      sig_sdl(18,:) = sig(18,:);
   elseif row == 5
      for ch = 4*Col - 1 : 5*Col - 3            % da 19 a 22
%          ch : (5*Col + 2) - k
         sig_sdl(ch,:) = sum(sig(ch : (5*Col + 2) - k,:));
         k = k + 1;
      end
      sig_sdl(23,:) = sig(23,:);
   elseif row == 6
      for ch = 5*Col - 1 : 6*Col - 3            % da 24 a 27
%          ch : (6*Col + 2) - k
         sig_sdl(ch,:) = sum(sig(ch : (6*Col + 2) - k,:));
         k = k + 1;
      end
      sig_sdl(28,:) = sig(28,:);
   elseif row == 7
      for ch = 6*Col - 1 : 7*Col - 3            % da 29 a 32
%          ch : (7*Col + 2) - k
         sig_sdl(ch,:) = sum(sig(ch : (7*Col + 2) - k,:));
         k = k + 1;
      end
      sig_sdl(33,:) = sig(33,:);
   elseif row == 8
      for ch = 7*Col - 1 : 8*Col - 3            % da 34 a 37
%          ch : (8*Col + 2) - k
         sig_sdl(ch,:) = sum(sig(ch : (8*Col + 2) - k,:));
         k = k + 1;
      end
      sig_sdl(38,:) = sig(38,:);
   elseif row == 9
      for ch = 8*Col - 1 : 9*Col - 3            % da 39 a 42
%          ch : (9*Col + 2) - k
         sig_sdl(ch,:) = sum(sig(ch : (9*Col + 2) - k,:));
         k = k + 1;
      end
      sig_sdl(43,:) = sig(43,:);
   elseif row == 10
      for ch = 9*Col - 1 : 10*Col - 3           % da 44 a 47
%          ch : (10*Col + 2) - k
         sig_sdl(ch,:) = sum(sig(ch : (10*Col + 2) - k,:));
         k = k + 1;
      end
      sig_sdl(48,:) = sig(48,:);
   elseif row == 11
      for ch = 10*Col - 1 : 11*Col - 3          % da 49 a 52
%          ch : (11*Col + 2) - k
         sig_sdl(ch,:) = sum(sig(ch : (11*Col + 2) - k,:));
         k = k + 1;
      end
      sig_sdl(53,:) = sig(53,:);
   elseif row == 12
      for ch = 11*Col - 1 : 12*Col - 4          % da 54 a 56
%          (ch + 1) : 12*Col - k
         sig_sdl(ch,:) = sum(sig((ch + 1) : 12*Col - k,:));
         k = k + 1;
      end
   end
end
