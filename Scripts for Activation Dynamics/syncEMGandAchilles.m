function [ EMG ] = syncEMGandAchilles( EMG, AchData)
%Function that synchronizes Achilles and EMG data

syncSignal = EMG.raw(:,132);
% [b,a] = butter(2,10/(2048/2),'low');
% syncSignal = filtfilt(b,a,syncSignal);
syncSignal(syncSignal == 0) = [];
syncSignal = syncSignal/max(syncSignal);
AchData = AchData/max(AchData);
% delay = finddelay(detrend(AchData/max(AchData)),detrend(syncSignal(1:end-maxIndex+50)/max(syncSignal(1:end-maxIndex+50))));       % Be aware that the index is not the best option
% delay = finddelay(AchData/max(AchData),syncSignal(1:end-maxIndex+50)/max(syncSignal(1:end-maxIndex+50)));       % Be aware that the index is not the best option
delay = finddelay(detrend(AchData),detrend(syncSignal(1:end)));
% delay = finddelay(AchData,syncSignal);

if (delay < 1) delay = 1; end
  figure, plot( 1000000*(AchData- mean(AchData))), hold on, plot(syncSignal(delay:(delay+length(AchData)-1))-mean(syncSignal(delay:(delay+length(AchData)-1))));
%    figure, plot( (AchData- mean(AchData))), hold on, plot(syncSignal(delay:(delay+length(AchData)-1))-mean(syncSignal(delay:(delay+length(AchData)-1))));

  legend('torque','sync')
  title('Sync. Achilles-EMG')
EMG.sync = EMG.raw(delay:(delay+length(AchData)-1),:);

end
