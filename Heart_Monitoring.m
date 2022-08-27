clear all
clc
close all

a=arduino();
ecg=0;

ubidotsURL = 'http://industrial.api.ubidots.com/api/v1.6/devices/heart-rate-monitor?token=BBFF-viydUhG7EoRjPx0gftf5jJb7F4eJbr';
options = weboptions('MediaType','application/json');

for k=1:500
    b=readVoltage(a,'A0');
    writeDigitalPin(a, 'D13', 1);
    ecg=[ecg,b]; //concatenating
    plot(ecg);
    %data = struct('ECG', ecg);
    %response = webwrite(ubidotsURL, data, options);
    grid on;
    drawnow;
end
writeDigitalPin(a, 'D13', 0);

%THINGSPEAK SETUP------------------------------
channelId = 1792534;
writeKey = 'MEDW4EK2HIVC8XYO';

%Creating timestamp----------------------------
tStamps = [datetime('now')-seconds(500):seconds(1):datetime('now')]';
dataTable = timetable(tStamps, ecg');
%thingSpeakWrite(channelId,dataTable,'WriteKey',writeKey);

ecg = ecg';
ecg=ecg-mean(ecg);
ecg=ecg/max(abs(ecg));

Fs = 150;
slen=length(ecg);
t=(1:slen)/Fs;

%% Low pass filtering %%

b=[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a=[1 -2 1] * 32;

%% filtering ECG signal using LPF %%

ecg_out1=filter(b,a,ecg);
ecg_out1=ecg_out1 - mean(ecg_out1);
ecg_out1=ecg_out1/max(abs(ecg_out1));

%% High pass filtering %%

b2=[-1,zeros(1,15),32,-32,zeros(1,14),1];
a2=[1,-1] * 32;

%% filtering the ECG signal from HPF

ecg_out2 =filter(b2,a2,ecg_out1);
ecg_out2= [ecg_out2(1:40)*0.25;ecg_out2(41:end)];
ecg_out2 = ecg_out2/max(abs(ecg_out2));

%% derivative operator %% 

b3= [2 1 0 -1 -2];
a3 = [1] *8;

%% filtering the ECG signal from derivative operator %%

ecg_out3 = filter(b3,a3,ecg_out2);
%% Normalization %%
ecg_out3=ecg_out3 - mean(ecg_out3);
ecg_out3 = ecg_out3/max(abs(ecg_out3));

%% squaring operation %%

ecg_out4 = ecg_out3.^2;
ecg_out4 = ecg_out4/max(abs(ecg_out4));

%% moving window integration operation %%

ecg_out4pad = [zeros(1,29) ecg_out4' zeros(1,29)];

for i=30:length(ecg_out4pad)-29
    ecg_int(i-29) = sum(ecg_out4pad(i-29:i))/30;
end

ecg5 = ecg_int';
ecg5 =ecg5/max(abs(ecg5));

%% Thresholding operation %%
TH = mean(ecg5)*max(ecg5); % Set threshold
w=(ecg5>(TH));

x=find(diff([0 w']) == 1); % Finding location of 0 to 1 transition
y=find(diff([w' 0]) == -1); % Finding location of 1 to 0 transition

%% cancelling the delay due to LOW PASS FILTER and HIGH PASS FILTER %%
x=x-(6+16); % 6 DELAY BY LPF & 16 DELAY BY HPF
y=y-(6+16);

%% Detect Q,R,S points %%
for i=1:length(x)
    %% R Locations %%
    [R_val(i),R_loc(i)]=max(ecg(x(i):y(i)));
    R_loc(i) = R_loc(i)-1 + x(i); % adding offset
    
    %% Q Locations %%
    [Q_val(i),Q_loc(i)]=min(ecg(R_loc(i):-1:R_loc(i)-8));
    Q_loc(i) = R_loc(i)-Q_loc(i)+1; % adding offset
    
    %% S Locations %%
    [S_val(i),S_loc(i)]=min(ecg(R_loc(i):R_loc(i)+10));
    S_loc(i) = R_loc(i)+S_loc(i)-1; % adding offset
end

figure(2)
plot(t,ecg);
xlim([0,10]);
grid on;
grid minor;
hold on;
plot(t(Q_loc),Q_val,'b+');
plot(t(R_loc),R_val,'r*');
plot(t(S_loc),S_val,'go');
xlabel('Time (seconds)')
ylabel('Signal Amplitude');
title('ECG Signal with QRS Locations after Thresholding');
legend('ECG','Q','R','S');
hold off;
%set(handles.edit3,'string',TH);
disp(['Threshold Value = ' num2str(TH) ' mV']);

%% Calculation of HEART RATE %%
HT=ceil((length(R_loc)*60)/t(end)); % calculate Heart rate
%% Calculation of QRS Duration %%
ecg6 = zeros(slen,1);
ecg6(w)=1;
ecg7 = diff([ecg6;0]);

x1= find(ecg7>0);
y1 = find(ecg7<0);

z1=y1-x1;
dur=mean(z1)*(1/Fs); % Calculate QRS duration

%% Calculation of RR interval %%
%RR_loc=[diff(R_loc),0]; % R-R Duration in Samples
%RR=RR_loc/Fs; % R-R Duration in Seconds
%R_val
RR=diff(R_val);
RR_sq=RR.^2;

HRV = sqrt(mean(RR_sq));

disp(['Heart Rate for ecg = ' num2str(HT) ' beats/min']);

disp(['Heart rate variability = ' num2str(HRV) ' s']);

disp(['QRS Duration for ecg = ' num2str(dur) ' sec']);

data = struct('HeartRate', HT, 'HRV', HRV, 'QRSduration', dur);
response = webwrite(ubidotsURL, data, options);