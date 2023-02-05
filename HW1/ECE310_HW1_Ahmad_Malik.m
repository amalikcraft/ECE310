%Ahmad Malik
%ECE310 HW1

clc;
close all;
clear;

%% Question 2
fprintf('Question 2\n');
%a)
f = 10e3;
fs = 50e3;
N = 256;
binSpace = fs/N;
fprintf('\ta)\n\t\tBin Spacing: %.2f Hz\n', binSpace);

%b)
k1 = round(f/binSpace);
k2 = round((fs-f)/binSpace);
fprintf('\tb)\n\t\tK1 = %d\n\t\tK2 = %d\n', k1, k2); 

%c)
%Offset for first k index
offset = 2*pi*(f - k1*binSpace)/fs;
%Straddle Loss
straddleLoss = abs(20*log10(diric(offset,250)));
fprintf('\tc)\n\t\tStraddle Loss: %.4f\n',straddleLoss);

%d)
%hamming windows padded with Zeros at end
Hamming = [hamming(250); zeros(6, 1)].';
% W(w')
w_prime = sum(exp(-1i*offset*(1:N)).*Hamming);
%|W(0)|
w0 = abs(sum(exp(0).*Hamming));
hamming_straddleLoss = abs(20*log10(abs(w_prime)/w0));
fprintf('\td)\n\t\tStraddle Loss using Hamming Window: %.4f\n',hamming_straddleLoss);

%e)
fprintf('\te)\n\t\t1024-point DTFT:\n');
N2 = 1024;
binSpace2 = fs/N2;
k = round(f/binSpace2);
%Offset for a k index
offset2 = 2*pi*(f - k*binSpace2)/fs;
%Straddle Loss
straddleLoss = abs(20*log10(diric(offset2,250)));
fprintf('\t\tStraddle Loss : %.4f',straddleLoss);
% W(w')
w_prime = sum(exp(-1i*offset2*(1:N)).*Hamming);
%|W(0)|
w0 = abs(sum(exp(0).*Hamming));
hamming_straddleLoss2 = abs(20*log10(abs(w_prime)/w0));
fprintf('\t\n\t\tStraddle Loss using Hamming Window : %.4f\n',hamming_straddleLoss2);

%% Question 4

%b)
%rectangular window:
N = 200;
Rec_Wind = ones(1,200);
%subsampling window:
Sub_Wind = mod((0:1:N-1),4)==0;
%Computing DFTS
DFT_RW = fft(Rec_Wind, N);
DFT_SW = fft(Sub_Wind, N);
%Stem Plots
figure;
subplot(2,1,1);
stem(1:N, abs(DFT_RW));
title('DFT of Rectangular Window')
xlabel('Frequency(Hz)');
subplot(2,1,2);
stem(1:N, abs(DFT_SW));
title('DFT of Subsampling Window')
xlabel('Frequency(Hz)');

%c)
N0 = 2^(nextpow2(16*N)); %useful function that finds next power of 2 for ffts
DFT_REC2 = fftshift(fft(Rec_Wind, N0));
DFT_SUB2 = fftshift(fft(Sub_Wind, N0));
radAxis = linspace(-pi,pi,N0);
%Magnitude Spectra Plots
figure;
subplot(2,1,1);
plot(radAxis, abs(DFT_REC2/N)/max(DFT_REC2/N));
title('Magnitude Spectra:  Rectangular Window')
xlabel('Frequency(Radians)');
subplot(2,1,2);
xlim([-pi,pi]);
plot(radAxis, abs(DFT_SUB2/N)/max(DFT_SUB2/N));
title('Magnitude Spectra:  SubSampling Window')
xlabel('Frequency(Radians)');
xlim([-pi,pi]);

%d)
A = 1;
Wo = 0.4;
N1 = 0:N-1;
%X samples
x = A*cos(Wo*N1); 
%w[n] and w0[n] of x[n]
X_Sub = x.*Sub_Wind;
X_Rec = x.*Rec_Wind;
%DFTs
X_DFT_Sub = abs(fft(X_Sub, N0));
X_DFT_Rec = abs(fft(X_Rec, N0));

omega = linspace(0,pi,N0);
figure;
plot(omega,abs(X_DFT_Rec)/N);  
hold on 
plot(omega,abs(X_DFT_Sub)/N);
xlim([0,pi]);
title('4096 - Point DFT of X windowed by w[n] and w0[n]');
xlabel('Frequency(Radians)');
ylabel('Magnitude Spectra');
legend('w[n]','W0[n]');

%e)
%Qualitatively, there is a difference between the spectral resolution of both graphs.

%% Question 5

%a)
load handel.mat;
t = 8; %seconds
y1 = y(1:Fs*t);
newY = reshape(y1, Fs, t);
%Graph of the DFT
figure;
for i = 1:8
    hold on
    data = abs(fft(newY(:,i), Fs)/Fs);
    plot(data(1:(Fs/2)));
end
title('DFT of the block window');
xlabel('Frequency (Hz)');
%peaks around integer multiples of 700Hz

%b)
%Creating the hamming window
f=0:1/2:Fs/2;
%Welsh periodogram 
p = 20*log10(abs(pwelch(y, Fs/2, Fs/4,  Fs*2)));
figure;
plot(f, p);
title('Welsh periodogram');
ylabel('decibels');
xlabel('Frequency (Hz)');

%c)
%Smiliar to the DFT of the individual segments, the Periodogram also has 
%peaks close to integer multiples of 700hz.

%d)
%spectrogram calculation
figure;
spectrogram(y, hamming(Fs/2), Fs/8, Fs);
title('Spectrogram')
