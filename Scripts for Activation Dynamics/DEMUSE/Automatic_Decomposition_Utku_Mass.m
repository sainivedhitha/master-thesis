clc
clear all;
close all;

[Name, path] = uigetfile;
dirlist2 = dir(path);

if ~isempty(dir([path,'abstract.xml'])),
    abs=xmlread([path, 'abstract']);
    abstract=xmlwrite([path, 'abstract']);
    ind1=strfind(abstract,'<protocol_code>');
    ind2=strfind(abstract,'</protocol_code>');
    newfilename=abstract(ind1+length(['<protocol_code>']):ind2-1);
else
    disp('default sample rate and gain was used. please check!')
    fsamp=2048;
    gain=1000;
end

exFactor = 10; 
LowFreq = 50;
HighFreq = 500;
DecompRuns = 70;
signalQuality = 5;
DiffModeChkbox = 0;

muscles={'GASlat','GASmed','PER','SOL','TA'};


Finterf =50;
%fsamp=2048;

[BF,AF] = butter(4,[LowFreq HighFreq]*2/fsamp);


flag=0;


for dirindex2 = 1 : size(dirlist2,1),
    
         if ~isempty(strfind(dirlist2(dirindex2).name,'.sig')), %&& isempty(strfind(dirlist2(dirindex2).name,'_decomposed')),
             flag=1;

             [path,dirlist2(dirindex2).name]

             h = fopen([path,dirlist2(dirindex2).name]);
             SIG = fread(h, [258 fsamp*inf], 'short');
             fclose(h);
             %gain=1000;
             DINAMICA = 10;			%	-5V +5V
             PRECISION = 12; 		%	number of bits of the A/D converter
             SIG = SIG.*DINAMICA/((2^PRECISION)*gain).*1e6; % The data are now converted in microvolt

             
             force = SIG(257,:);

             for i=1:length(muscles),
                 
             muscles(i)    
                 if i==1,EMG=SIG(33:64,:);end    %GASlat
                 if i==2,EMG=SIG(129:192,:);end  %GASmed
                 if i==3,EMG=SIG(1:32,:);end     %PER
                 if i==4,EMG=SIG(65:128,:);end   %SOL
                 if i==5,EMG=SIG(193:256,:);end  %TA
                 
                Ni=find(dirlist2(dirindex2).name=='.');
                Pos_MVC=dirlist2(dirindex2).name(1:Ni-1)
                
                switch Pos_MVC
                    case 'Anatomy-30'
                        epoch=[16,29,30,44,45,57,58,72];
                    case 'Anatomy-50'
                        epoch=[12,25,26,39,40,52,53,67];
                    case 'Anatomy-70'
                        epoch=[14,27,28,41,42,54,55,69];
                    case 'Anatomy-90'
                        epoch=[17,30,31,44,45,59,60,73];
                    case 'Dorsi-30'
                        epoch=[24,38,39,52,53,66,67,81,82,94,95,108];
                    case 'Dorsi-50'
                        epoch=[11,24,25,38,39,52,53,66];
                    case 'Dorsi-70'
                        epoch=[11,24,25,38,39,52,53,66];
                    case 'Dorsi-90'
                        epoch=[13,26,27,40,41,54,55,68];
                    case 'Plantar-30'
                        epoch=[18,31,32,45,46,59,60,73];
                    case 'Plantar-50'
                        epoch=[13,26,27,40,41,54,55,68];
                    case 'Plantar-70'
                        epoch=[11,24,25,38,39,52,53,67];
                    case 'Plantar-90'
                        epoch=[13,26,27,40,41,54,55,68];
                end

                
                if size(EMG,2)<size(EMG,1), EMG = EMG';end
                
                sigLength=length(EMG);
                Frad=round(fsamp/200);
                fInt=round(0:Finterf*sigLength/fsamp:(sigLength/2));
                for ch = 1 : size(EMG,1)                   
                    EMG(ch,:)=real(subtractFreq2(EMG(ch,:),fInt,Frad,0));
                end
                EMG=filtfilt(BF,AF,EMG')';
                disp('Filtering DONE!!!');
                tic
                step=5; %in second
                flag2=0;
                for e=1:2:length(epoch)-1;              %1:step*fsamp:length(EMG)-(step*fsamp),
                    flag2=flag2+1;
                    
                    z={dirlist2(:).name};
                    chkfolder=strfind(z,[muscles{i},dirlist2(dirindex2).name(1:end-4),'_decomposed',num2str(flag2)]);
                    if isempty(cell2mat(chkfolder)),

                    
                    [ActIndWhole,ActIndResid,MUPulses,Cost,MUPulses2,Cost2,IPTs,IPTs2,ProcTime] = ...
                        GUICKCgrad5(EMG(:,epoch(e)*fsamp:epoch(e+1)*fsamp),exFactor,fsamp,DecompRuns,5,0,0); % run the decomposition
                                                    
                    [MUPulses,newInd]=eliminateDuplicateIPT(MUPulses,fsamp,1,0.3);
                    Cost = {Cost{newInd}};
                    toc
                    
                    IPTsAll = [IPTs IPTs2];
                    MUFirings = [MUPulses MUPulses2];
                    CostAll = [Cost Cost2];
                    
                    save([path,muscles{i},dirlist2(dirindex2).name(1:end-4),'_decomposed',num2str(flag2)],'MUFirings','MUPulses','MUPulses2','Cost','Cost2','CostAll','IPTs','IPTs2','IPTsAll','EMG','force');
                    clear MUFirings MUCost Cost Cost2 CostAll MUPulses MUPulses2 IPTsAll IPTs IPs2
                    [path,muscles{i},dirlist2(dirindex2).name(1:end-4),'_decomposed',num2str(flag2)]
                    else
                    disp([muscles{i},dirlist2(dirindex2).name(1:end-4),'_decomposed',num2str(flag2),' ALREADY DECOMPOSED!'])
                    end
                end
                clear EMG 

             end    
        end
        flag=0;
end
  

% if flag ==1, disp('Decomposed DONE!!!');
% else disp('Nothing NEW!!!!!!')   
% end
