%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 362 Project Template File
% student name:
% student number:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all, clear, clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import audio data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = importdata('ENGR_362_guitar_Fs_is_48000_Hz.txt');
Fs = 48000;                             % sampling freq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Play sound of raw audio data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% soundsc(y,Fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graph in time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ts = 1/Fs;                              % sampling period       
Length_y = length(y(:,1));              % length of signal
time = (0:Length_y-1)*Ts;               % time vector
figure,plot(time,y), axis tight
title('y(t) vs t')                      % labels
xlabel('t (s)')                         % labels
ylabel('y(t)')                          % labels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graph with DFT/FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = fft(y);                             % Discrete Fourier transform
F1 = abs(Y/Length_y);                   % frequency
F2 = F1(1:Length_y/2+1);                % half of frequency
F2(2:end-1) = 2*F2(2:end-1);            % Discrete Fourier transform
f = Fs*(0:(Length_y/2))/Length_y;       % freq vector [Hz]
f_kHz = f/1000;                         % freq vector [kHz]
figure,plot(f_kHz,F2)                   % plot 
axis([0  max(f_kHz)/25 0 max(F2)])         % axis details
title('Y(F) vs F')                           % labels
xlabel('F (kHz)')                       % labels
ylabel('Y(F)')                          % labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D major chord frequencies for notes D3, A3, D4, F#4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D3 = 146.83;                            % freq of note D3 [Hz]
D3_int = round(D3/max(f)*length(f))     % associated integer to above freq
A3 = 220.00;                            % freq of note A3 [Hz]
A3_int = round(A3/max(f)*length(f))     % associated integer to above freq
D4 = 293.66;                            % freq of note D4 [Hz]
D4_int = round(D4/max(f)*length(f))     % associated integer to above freq
F_sharp_4 = 369.99;                     % freq of note F#4 [Hz]
F_sharp_4_int = ...
    round(F_sharp_4/max(f)*length(f))   % associated integer to above freq

note_freq = [D3 A3 D4 F_sharp_4];       % vector of all note freqs
note_freq_int = ...
  [D3_int A3_int D4_int F_sharp_4_int]; % vector of all int note freqs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop to apply filter bank

N = 7; % define the order
R = 2; % define the ripple
Limit = 5; % define the upper or lower limit
Y_filtered = []; % create a variable to save filtered siganals
F2_s = []; % create a  variable to save discrete Fourier transform of each filtered signal to calculate amplitudes

for fre = [146.83 220 293.66 369.99]
    [b1,a1] = cheby1(N,R,(fre+Limit)/Fs*2,'low');  %low-pass filter
    y1_filtered = filter(b1, a1, y);
    figure,freqz(b1,a1),hold on  %Generate a graph
    [b2,a2] = cheby1(N,R,(fre-Limit)/Fs*2,'high'); %high-pass filter
    y2_filtered = filter(b2, a2, y1_filtered);
    freqz(b2,a2),hold off

    soundsc(y2_filtered,Fs) %play the audio of the filtered signal

    Y2_filtered = fft(y2_filtered);         % Discrete Fourier transform
    Y_filtered = [Y_filtered,Y2_filtered ]; % Save the filtered signal
    F1 = abs(Y2_filtered/Length_y);         % frequency
    F2 = F1(1:floor(Length_y/2+1));         % half of frequency
    F2(2:end-1) = 2*F2(2:end-1);            % Discrete Fourier transform
    F2_s = [F2_s, F2];                      % Save the above result of each filtered signal
    f = Fs*(0:(Length_y/2))/Length_y;       % freq vector [Hz]
    figure,plot(f,F2)                   % plot 
    axis([0  max(f)/50 0 max(F2)])         % axis details

    dec1 = fre;
    str1 = 'Y(F) for ';
    str2 = num2str(dec1);
    str3 = 'Hz';
    str4 = ' vs F';
    Title = [str1,str2,str3,str4];  % Create a title for each signal
    Y_label = [str1,str2,str3];     % Create a y-axis label for each signal
    title(Title);                           % labels
    xlabel('F (Hz)');                      % labels
    ylabel(Y_label);                       %labels
    
end


%Retrieve Amplitudes of Each Filtered Signal

F2_D3 = F2_s(:,1);
amplitude_D3 = abs(max(abs(F2_D3(D3_int-10:D3_int+10))));

F2_A3 = F2_s(:,2);
amplitude_A3 = abs(max(abs(F2_A3(A3_int-10:A3_int+10))));

F2_D4 = F2_s(:,3);
amplitude_D4 = abs(abs(max(F2_D4(D4_int-10:D4_int+10))));

F2_F4 = F2_s(:,4);
amplitude_F4 = abs(abs(max(F2_F4(F_sharp_4_int-10:F_sharp_4_int+10))));


%Balance amplitudes to the max amplitude among the filtered signals and Recombine them;
max_amplitude = max([amplitude_A3,amplitude_D3,amplitude_D4,amplitude_F4]);
Y_merged = Y_filtered(:,1)/amplitude_D3.*max_amplitude + Y_filtered(:,2)/amplitude_A3.*max_amplitude + Y_filtered(:,3)/amplitude_D4.*max_amplitude + Y_filtered(:,4)/amplitude_F4.*max_amplitude ;

%Generate a graph of the final compressed signal in the frequency domain

F1 = abs(Y_merged/Length_y);             % frequency
F2 = F1(1:floor(Length_y/2+1));         % half of frequency
F2(2:end-1) = 2*F2(2:end-1);            % Discrete Fourier transform
f = Fs*(0:(Length_y/2))/Length_y;       % freq vector [Hz]
figure,plot(f,F2)                       % plot 
axis([0  max(f)/50 0 max(F2)])          % axis details
title('Y merged vs F')                  % labels
xlabel('F (Hz)')                        % labels
ylabel('Y merged')                      % labels



