function [spikeTrains,dischargeRates,averageDR,stdDR,IDI,averageIDI,ax] = getSpikeTrains(muPulses, tVec, fs, force, plotFlag)


% Graph settings
color = {'#A2142F','#D95319','#EDB120','#4DBEEE','#0072BD'};
colSeq=1; % initializing color sequence
linWidth = 0.5;

% if plotFlag
%     figure;
% end

% Initializing variables
spikeTrains = zeros(length(tVec),length(muPulses));
dischargeRates = NaN(length(tVec),length(muPulses));
%% For every motor unit, plots the spike trains and discharge rates
for muCount = 1:length(muPulses)
    spikeTrains(muPulses{muCount},muCount) = 1;
    dischargeRates(muPulses{muCount}(2:end),muCount) = fs./diff(muPulses{muCount});
    dischargeRates(dischargeRates<3 | dischargeRates>50) = NaN; % values below 3 hz may be due to big pauses (non physiological 
    spikeTrains(dischargeRates<3 | dischargeRates>50) = 0; % same notice this eliminates some spikes
    if plotFlag
        % to plot spikes neatly, we need to give two values per instance
        % (sample), i.e. we need to go from 0 to 1 in a single sample,
        % hence the elements (samples) of t2 are repeated and the spike 
        % locations are updated 
        t2= repelem(tVec,2); % duplicating the samples of t vector  
        spikePos = (muPulses{muCount})*2; % updating spike locations
        % Plot variable is full of NaNs to leave blank the samples without
        % spikes
        plotVar = NaN(size(t2)); % creating plot variable
        plotVar(spikePos) = muCount-0.4;
        plotVar(spikePos+1) = muCount+0.4;
        plotVar = plotVar(1:length(t2)); % creating plot variable       
        ax1=subplot(1,2,1);
        plot(t2,plotVar,'Color', color{colSeq},'LineWidth',linWidth); %,
        title('MU spike trains')
        hold on
        % to check innervation pulse trains(IPTs1 needs to be added as inpu
        %         plot(tVec(1:end-1),0.8*IPTs1{trialCount}/max(IPTs1{trialCount})+trialCount-0.4,'Color', color{c+1})
        %         hold on
        ax2=subplot(1,2,2);
        plot(tVec, dischargeRates(:,muCount)/30+muCount-0.5, 'o', 'MarkerSize',5,  'MarkerEdgeColor',color{colSeq},...
        'MarkerFaceColor',color{colSeq}); % normalized to 30 pps
        title('MU discharge rates')
        hold on
        if colSeq < 5
            colSeq = colSeq + 1 ;
        else
            colSeq=1;
        end
    end
end 
averageDR=mean(dischargeRates,'omitnan');
stdDR=std(dischargeRates,'omitnan');
IDI=1./dischargeRates;
averageIDI=mean(IDI,'omitnan');

if plotFlag
    ax= [ax1,ax2];

    subplot(1,2,1)
    xlim([-1, max(tVec)+1]);
    ylim([0 length(muPulses)+0.5]);
    yticks(0:1:length(muPulses)+1)
    xlabel('Time (s)')
    ylabel('Motor unit (#)')
    subplot(1,2,2)
    xlim([-1, max(tVec)+1]);
    ylim([0 length(muPulses)+0.5]);
    yticks(0:1:length(muPulses)+1)
    yticklabels({'',num2str(averageDR(:),'%.1f')});
    xlabel('Time (s)')
    ylabel('Mean discharge rate (pps)')
    if ~isempty(force)
        for i=1:2
            yyaxis(ax(i),'right')
            ax(i).YAxis(2).Color = 'k';
            hold on, plot(ax(i),tVec,force,'Color',[0 0 0]);
            ylabel(ax(i),'Force (N)')
        end
  %      to do it manually, only if the force is desired to be in the back
% 
%         ax(1).Color= "none";
%         ax1_pos = ax(1).Position; % position of first axes
%         ax2 = axes('Position',ax1_pos,...
%             'YAxisLocation','right',...
%             'Color','none');
%         uistack(ax1, 'bottom')
%         linkaxes([ax(1),ax2],'x')
%         set(gca,'xtick',[])
%         ax2.Color = 'white';
%         
    end

    linkaxes(ax)
else
    ax=[];
end
