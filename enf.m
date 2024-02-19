function [y1, y2] = enf(x, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency)
%% Parameters
%sampleSizeSeconds = 16; %Seconds per Sample
%padRat = 2; %Padding Ratio
%[x,fs] = audioread('groundtruth1.wav'); %.wav file, Vals in x, fs is sampling freq of .wav



k1 = Frequency-1; %frequency min
k2 = Frequency+1; %freq max

%overlap_rate = 0.5; %Overlap rate between [0-0.99]
%% Input Segmentation 
%segment_size = sampleSizeSeconds*fs; %Window in #of Samples

padL = ZeroPad; %Zero Padding based on raio of window size

d = round(BlockSize*(1-Overlap)); %Overlap index
idx=(1:BlockSize)'+(0:d:(numel(x)-BlockSize));%Overlapped index for window (not padded)

matrix = x(idx); %Values of .wav file @fs in overlapped window

%% Windowing Functions


% Gas= gausswin(BlockSize);
% Ham= hamming(segment_size,"periodic");
% Bm = blackmanharris(segment_size, "periodic"); %Chose BMH because of balanced shape

Bmmatrix = Window.*matrix; %Applying Window Function


%% Zero Padding
PadMat = padarray(Bmmatrix, padL, 0, 'post'); %Zero Padding Windowed Segmented Matrix
%% FFT
FtMat = fft(PadMat);%FFT
%AbsFtMat = abs(fft(PadMat));%Abs FFT

[k,m] = size(FtMat); %idk
Eidx= abs(FtMat); %ENF

[K,M] = size(Eidx); %Compute Row, Col Size after Padding
M = 1:M; %Turn M into Range
t = 1:k; %Turn K into Range t
%% Ext Freq of Interest
Fidx(t) = Fs*(t-1)/(BlockSize + padL); %Frequencies of entire Padded FFT @fs

Fint = find((k1<Fidx)&(Fidx<k2)); %Indexes of Fidx @ Freq of Interest [k1,k2]
Fval = Fidx((k1<Fidx)&(Fidx<k2)); %Frequency Values in Region of Interest

var = Eidx(Fint,M); %Energy @ frequencies of interest using index value of Fidx
%% Max and Weighted Energy
[~,Midx] = max(var,[], 1); %Max energy value(unused), and index across all Freq Of Int
MFreqs = Fval(Midx); %Freq Val @ max energy per block

[Freq,blk]= meshgrid(Fval,M); %Freq Grid between K1 and K2 mapped to Blocks


vsum = sum(var); %Sum of Energies across Column
Nsum = sum(var.*(Freq.')); %Product of Energy and Frequencies Across Column(Blocks)

WENF = Nsum./vsum; %Weighted Energy Calculation

%% Plots 
%Spectrogram
figure 
surf(Freq, blk, var.'); %Freq Values, Windows, and Energy
title(['ENF Spectrogram @[',num2str(k1),'Hz-',num2str(k2),'Hz]']);
ylabel('Block Number') ;
zlabel('Energy (Joules)') ;
xlabel('ENF Freq (Hz)') ;
figure

%Max ENF
subplot(2,1,1);
plot(blk,MFreqs); %block and MENF
title(['ENF MaxFreq per Window @[',num2str(k1),'Hz-',num2str(k2),'Hz]']);
grid on;
xlabel('Block Number');
ylabel('Max ENF (Hz)'); 

%Weighted ENF
subplot(2,1,2);
plot(blk,WENF); %block and WENF
title(['Weighted Avg of @[',num2str(k1),'Hz-',num2str(k2),'Hz]']);
grid on;
xlabel('Block Number'); 
ylabel('Weighted ENF (Hz)'); 

y1 = MFreqs;
y2 = WENF;

end