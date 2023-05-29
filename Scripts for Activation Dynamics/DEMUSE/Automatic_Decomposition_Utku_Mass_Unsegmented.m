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
LowF*req = 50;
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
    
         if ~isempty(strfind(dirlist2(dirindex2).name,'.sig')) && isempty(strfind(dirlist2(dirindex2).name,'_decomposed')),
             flag=1;

             [path,dirlist2(dirindex2).name]

             h = fopen([path,dirlist2(dirindex2).name]);
             SIG = fread(h, [258 fsamp*inf], 'short');
             fclose(h);
             %gain=1000;
             DINAMICA = 10;			%	-5V +5V
             PRECISION = 12; 		%	number of bits of the A/D converter
             SIGall = SIG.*DINAMICA/((2^PRECISION)*gain).*1e6; % The data are now converted in microvolt
             %load([path,'\',dirlist(dirindex).name,'\',dirlist2(dirindex2).name]);
             
             force = SIG(257,:);

             for i=1:length(muscles),

                 if i==1,EMG=SIGall(33:64,:);end    %GASLat
                 if i==2,EMG=SIGall(129:192,:);end  %GASMed
                 if i==3,EMG=SIGall(1:32,:);end     %PER
                 if i==4,EMG=SIGall(65:128,:);end   %SOL
                 if i==5,EMG=SIGall(193:256,:);end  %TA

                
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
                
                [ActIndWhole,ActIndResid,MUPulses,Cost,MUPulses2,Cost2,IPTs,IPTs2,ProcTime] = ...
                    GUICKCgrad5(EMG,exFactor,fsamp,DecompRuns,5,0,0); % run the decomposition
                
                [MUPulses,newInd]=eliminateDuplicateIPT(MUPulses,fsamp,1,0.3);
                Cost = {Cost{newInd}};
                
                MUFirings = [MUPulses MUPulses2];
                CostAll = [Cost Cost2];
                IPTsAll=[IPTs IPTs2];
                
                toc
                %MUFirings = MUPulses;
                %MUCost = Cost;

            %MUFiring

                %save([path,muscles{i},dirlist2(dirindex2).name(1:end-4),'_decomposed'],'MUFiring','MUCost','EMG','force');
                save([path,muscles{i},dirlist2(dirindex2).name(1:end-4),'_decomposed'],'MUFirings','MUPulses','MUPulses2','Cost','Cost2','CostAll','IPTs','IPTs2','IPTsAll','EMG','force');
                %delete([path,'/',dirlist(dirindex).name,'/',dirlist2(dirindex2).name]);

                clear MUFirings MUCost Cost Cost2 CostAll MUPulses MUPulses2 IPTsAll IPTs IPs2 
                
                [path,muscles{i},dirlist2(dirindex2).name(1:end-4),'_decomposed']

             end    
        end
        flag=0;
end
  

% if flag ==1, disp('Decomposed DONE!!!');
% else disp('Nothing NEW!!!!!!')   
% end
