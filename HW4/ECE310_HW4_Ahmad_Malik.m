%Ahmad Malik
%ECE 310 - Fall 2022
%Fontaine
%Problem Set 4

clc
clear 
close all

%% Question 1

%a)
B = 800;
fso = 2e3;
M = 16;
fs = M*fso;

% Generate bandlimited gaussian data 
N = 2^(nextpow2(fs)); % Will use N as number of samples to generate
Nfft= 2*N;

%-what is the formula for the k0 DFT index that is closest to B, given radix Nfft DFT?
k0 = round(B*Nfft/fs);

xf = randn(1,k0+1); %% Generates iid N(0,1) samples, which will occupy DFT indices 0<=k<=k0
x = real(ifft(xf,Nfft)); %% Indices k0+1<=k<=Nfft-1 are all 0
x = x(1:N);
x = x/std(x); % Normalize to unit variance

%b)
%Q: Why is the signal in the frequency domain (xf extended with 0's to length N) not stationary? 
%A: Because the variance is not constant throughout the data.

%Q: Why is x stationary?
%A: Because its is the Fourier Transform of IID Gaussian samples, which is also stationary.

%Q: Why is x Gaussian and 0-mean?
%A: Because the fourier transform is linear and X consists of jointly Gaussian Random Variables.

%c)
[psd,w] = pwelch(x);
f = w*(fs/(2*pi));
figure;
plot(f, psd);
xline(B,'-r',{'B = 800Hz'});
title('Estimated Power Spectrum of X');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 1500]);

%d)
%Q: Why does ifft(xf,Nfft) produce a complex but symmetric result?
%A: From the IDFT formula, there is conjugate complex symmetry where X_N-n = X*n, 
%   therefore it is possible that for even real values of X, the output can be complex.
%   So the ifft of ra real signal in the frequency domain is a complex symmetric signal in the time domain.

%e)
%Q: Why do I say the operation x = x/std(x) normalizes to unit variance?
%A: Assume Gaussian, and it is unit variance because we normalize it (x/std(x) by ergodicity.

%f)
%Q: If time-domain stationary Gaussian samples occupy a limited bandwidth (here up to B only), 
%   why are they necessarily correlated?

%A: We can assume it is correlated because the psd is not flat within the 800hz bandwidth

%% Question 2

Ns = 20;
Qbits = 3;
scale = 3;
x1 = 5*(2*rand(1,Ns)-1);
xq = quantop(x1,Qbits,scale);

figure;
hold on;
stem(x1);
stem(xq);
xlabel("Samples");
ylabel("Magnitude");
title("Comparing Original and Quantized Values");
legend("Original ","Quantized");
hold off;

%Plot shows that the quantized values are comparably smaller in magnitude
%than the original values. The maximum value is 2.25.

%The reason we add 1 then subtract 1 is because after scaling X, we want to
%make x values positive for rounding since negative values are rounded with bias. 
%After rounding x to between [0,2), we move x to [-1, 1] by subtracting 1,
%and reverse scaling.



%% Question 3

%a)

%Wp = (7/8)Wc = (7/8)(pi/M) = (7/8)(pi/16) = (7/128)*pi
%                         Wp = 2*pi*(fp/fs) = (7/128)*pi
%                             2*pi*(fp/32khz) = (7/128)*pi
%                                 (fp = 875Hz) > (B = 800Hz)
fprintf("(fp = %fHz)  >  (B = %fHz)\n\n", 875, B);                            
%b)
% f_crit = [0, 7/8 , 1/16 , 1/16 , 1 ]
% acrit = [ 1 , 1 , 0 , 0 ]

%c)
M = 16;
wc = 1/M;
fcrit = [0 , (7/8)*wc , wc , 1];    
acrit = [1 , 1  , 0  , 0]; 

b0 = firpm(221,fcrit,acrit);
b0 = b0/norm(b0);
a0 = 1;  % normalize filter to unit power

[h,w] = freqz(b0,a0);
phase = unwrap(angle(h))*180/pi;
f2 = (w*fs)/(2*pi); 
 
figure;
subplot(1,2,1);
plot(f2,abs(h));
xline(B,'-r',{'B = 800Hz'});
title("Decimation Filter Magnitude Response");
xlabel("Frequency (Hz)");
ylabel("Magnitude Response (Linear)");
xlim([0 5000])
set(gcf, 'Position',  [100, 100, 1000, 500])

subplot(1,2,2);
plot(f2,phase);
xline(B,'-r',{'B = 800Hz'});
title("Decimation Filter Phase Response");
xlabel("Frequency (Hz)");
ylabel("Phase (deg)");
xlim([0 5000])
set(gcf, 'Position',  [100, 100, 1000, 500])

%Since the passband has ripples, the filter would have distortion and since
%the passband of the phase plot is linear, it would have no phase distortion

%% Question 4

%b)
e = [x(1), zeros(1, N-1)];
y = zeros(1, N);

for n = 2:1:N
    y(n) = e(n-1) + y(n-1);
    e(n) = x(n) - y(n);
end
%variance of error
e_var = var(e);
fprintf("Varriance of Error = %f \n", e_var);
%% Question 5

%a) 
prob = 0.01;
% alpha value CDF takes y, uniformly distributed on (0, 1), and defines x = Φ−1(y),
% where Φ is the CDF. The function norminv(prob) which can calculate Φ−1(y). 
alpha = norminv(prob);


%b
%System I
Qbits = 5;
scale = 10;
out_filter1 = downsample(filter(b0,a0,x),M);
out_filterq = quantop(out_filter1,Qbits,scale);
outputI = abs(fft(out_filter1,N));
outputIq = abs(fft(out_filterq,N));
%SNR System 1
SNRI = 10*log10(mean(outputI.^2)/mean((outputI - outputIq).^2));
fprintf("\nSNRI = %f \n", SNRI);

%plot System I with Quantization 
temp = fso/N*(0:N-1);
fsys = temp(1:N/2);

figure;
subplot(3, 1, 1);
plot(fsys,10*log10(abs(outputIq(1:N/2))));
xlabel("Frequency (Hz)");
ylabel("Magnitude (dB)");
title("System I: with Quantization");
set(gcf, 'Position',  [100, 100, 1000, 1000])

%--------------------------

%System II
Qbits = 5;
scale = 5;
x2q = quantop(x,Qbits,scale);
out_filter2q = downsample(filter(b0,a0,x2q),M);
outputIIq = abs(fft(out_filter2q,N));
%SNR System II
SNRII = 10*log10(mean(outputI.^2)/mean((outputI - outputIIq).^2));
fprintf("SNRII = %f \n", SNRII);

%plot System I and II: No Quantization 
subplot(3, 1, 2);
plot(fsys,10*log10(abs(outputI(1:N/2))));
xlabel("Frequency (Hz)");
ylabel("Magnitude (dB)");
title("System I and II: No Quantization");
set(gcf, 'Position',  [100, 100, 1000, 1000])

%plot System II with Quantization 
subplot(3, 1, 3);
plot(fsys,10*log10(abs(outputIIq(1:N/2))));
xlabel("Frequency (Hz)");
ylabel("Magnitude (dB)");
title("System II: with Quantization");
set(gcf, 'Position',  [100, 100, 1000, 1000])

%--------------------------

%System III
Qbits = 5;
scale = 0.5;

%No quantization
outfilter3 = downsample(filter(b0,a0,y),M);
outputIII = abs(fft(outfilter3,N));

yout = zeros(1,N);
err = quantop(x(1),Qbits,scale);
eq = [err,zeros(1,N-1)];
e1 = [x(1),zeros(1,N-1)];

for n=2:N
    yout(n) = eq(n-1) + yout(n-1);
    e1(n) = x(n) - yout(n);
    eq(n) = quantop(e1(n),Qbits,scale);
end

outfilterIIIq = downsample(filter(b0,a0,yout),M);
outputIIIq = abs(fft(outfilterIIIq,N));
% SNR III
SNRIII = 10*log10(mean(outputIII.^2)/mean((outputIII - outputIIIq).^2));
fprintf("SNRIII = %f \n", SNRIII);

figure;
subplot(2,1,1)
plot(fsys,10*log10(abs(outputIII(1:N/2))));
xlabel("Frequency (Hz)");
ylabel("Magnitude (dB)");
title("System III: No Quantization)");

subplot(2,1,2)
plot(fsys,10*log10(abs(outputIIIq(1:N/2))));
xlabel("Frequency (Hz)");
ylabel("Magnitude (dB)");
title("System III: with Quantization");



%% Question 6

%a) 

%SNRI = ~26
%SNRII = ~36
%SNRIII = ~53

%SNRIII > SNR II > SMR I

%SNR per Octave
SNR2v1_Oc = abs(SNRII - SNRI)/log2(M);
SNR3v1_Oc = abs(SNRIII - SNRI)/log2(M);

fprintf("\n\nSNR benefit per octave of System II over System I, = %f\n", SNR2v1_Oc);
fprintf("SNR benefit per octave of System III over System I, = %f\n", SNR3v1_Oc);


%b
SNR2v1bits = abs(SNRII - SNRI)/6;
SNR3v1bits = abs(SNRIII - SNRI)/6;

fprintf("\n\n~%.0f bit of extra precision will be required for system 2 over System 1\n",floor(SNR2v1bits));
fprintf("~%.0f bit of extra precision will be required for system 3 over System 1",floor(SNR3v1bits));

%c) yes, the 3db per octave benefit from system II is present.

%d) AWGN



function xq = quantop(x,Qbits,scale)
    x = x/scale;
    xq = round((x+1)*2^(Qbits-1)) / 2^(Qbits-1) - 1;
    xq = max(min(xq,1-2^-(Qbits-1)),-1);
    xq = xq*scale;
end

