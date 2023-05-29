function Sig_OK = Matrix_reconstr(Electrode_Map,Signal,Acquisition_mode)

%% ------------------------------------------------------------------------------------ %
%                 Matrix reconstruction  -  Single differential signals                 %
% ------------------------------------------------------------------------------------- %
clear Sig_OK
if strcmpi(Acquisition_mode,'Differential')
    % ricostruisco i canali uno per uno
    for columna = 1:size(Electrode_Map,1) % ciclo sulle colonne
        for rowa = 1:size(Electrode_Map,2)-1
            if not(isnan(Electrode_Map(columna,rowa)-Electrode_Map(columna,rowa+1)))
                if Electrode_Map(columna,rowa)<Electrode_Map(columna,rowa+1) % caso di somme negative
                    vett_somme{columna}{rowa} = -Electrode_Map(columna,rowa);
                    for ind_sum = 1:Electrode_Map(columna,rowa+1)-Electrode_Map(columna,rowa)-1;
                        vett_somme{columna}{rowa} = [vett_somme{columna}{rowa} -Electrode_Map(columna,rowa)-ind_sum];
                    end
                else
                    vett_somme{columna}{rowa} = Electrode_Map(columna,rowa+1);
                    for ind_sum = 1:Electrode_Map(columna,rowa)-Electrode_Map(columna,rowa+1)-1;
                        vett_somme{columna}{rowa} = [vett_somme{columna}{rowa} Electrode_Map(columna,rowa+1)+ind_sum];
                    end
                end
                Sig_OK{columna}(rowa,:) = zeros(1,size(Signal,2));
                if sum(vett_somme{columna}{rowa})<inf
                    for ind_sum = 1: size(vett_somme{columna}{rowa},2)
                        Sig_OK{columna}(rowa,:) = Sig_OK{columna}(rowa,:) + sign(vett_somme{columna}{rowa}(ind_sum))*Signal(abs(vett_somme{columna}{rowa}(ind_sum)),:);
                    end
                else
                    Sig_OK{columna}(rowa,:) = zeros(1,size(Signal,2))*nan;
                end
            end
        end
    end
else % 'Monopolar'
    %% ----------------------------------------------------------------------------- %
    %                 Matrix reconstruction  -  Monopolar  signals                   %
    % ------------------------------------------------------------------------------ %

    for columna = 1:size(Electrode_Map,1) % ciclo sulle colonne
        for rowa = 1:size(Electrode_Map,2)-1
            if not(isnan(Electrode_Map(columna,rowa)-Electrode_Map(columna,rowa+1)))
                Sig_OK{columna}(rowa,:) = Signal(Electrode_Map(columna,rowa+1),:) - Signal(Electrode_Map(columna,rowa),:);
            end
        end
    end
    
end

