%Ahmad Malik
%ECE310 HW2

clc;
close all;
clear;

%% Question 1

%a)
c = 1/8;
Y = [1.0688, .4505, .4594, 1.0688];
X = [1, -1.5055, 1.2630, -.3778];
figure();
zplane(Y*c,X);
title('Pole-Zero Plot of H(z)');

%b)
c2 = 1/2;
y1 = [-0.4954 1];
x1 = [1 -0.4954];
y2 = [0.7626 -1.0101 1];
x2 = [1 -1.0101 0.7626];
%all pass transfer function
Yall = (conv(y1, x2)+conv(y2, x1));
Xall = conv(x1,x2);

%exact poles and zeros of H(z) and H(z) Allpass
[zeros,poles] = tf2zpk(Y*c,X);  
%H(z):   
%zeros [-0.9968 ; 0.2876 + 0.9594i ; 0.2876 - 0.9594i]
%poles [.5050 + 0.7124i ; 0.5050 - .7124i ; 0.4954]
[zeros_all,poles_all] = tf2zpk(Yall*c2,Xall);
%H(z) Allpass:
%zeros [-1.000 ; 0.2895 + 0.9572i ; 0.2895 - 0.9572i]
%poles [.5051 + 0.7124i ; 0.5051 - .7124i ; 0.4954]

figure();
subplot(2,1,1);
zplane(Y*c,X);
title('Pole-Zero Plot of H(z)');
subplot(2,1,2);
zplane(Yall*c2,Xall);
title('Pole-Zero Plot of H(z) using Allpass');

%The graph and the calculated zeros/poles using tf2zpk show that the
%transfer functions are the same.

%c)
h = freqz(Y*c, X);
Hdb = 20*log10(abs(h));
Hph = 180*(unwrap(angle(h))./pi);
w = linspace(0, 2*pi, 512);

figure;
subplot(2, 1, 1)
plot(w, Hdb);
title('Magnitude Response of H(z)');
ylabel('Magnitude (dB)');
xlabel("Frequency (radians)");
xlim([0 2*pi]);
xticks([0 pi/2 pi 3*pi/2 2*pi]);
xticklabels({'0', '\pi/2', '\pi', '3\pi//2', '2\pi',})
subplot(2, 1, 2)
plot(w, Hph);
title('Phase Response of H(z)');
ylabel('Magnitude (dB)');
xlabel("Frequency (radians)");
xlim([0 2*pi]);
xticks([0 pi/2 pi 3*pi/2 2*pi]);
xticklabels({'0', '\pi/2', '\pi', '3\pi//2', '2\pi',})

%d)
hf1 = freqz(y1,x1*2);
hf2 = freqz(y2,x1*2);
Hf1_ph = 180*(unwrap(angle(hf1))./pi);
Hf2_ph = 180*(unwrap(angle(hf2))./pi);
figure();
plot(w,Hf1_ph);
hold on;
plot(w, Hf2_ph);
title("Phase Response of Allpass Factors");
ylabel("Phase (deg)");
xlabel("Frequency (rad)");
xlim([0 2*pi]);
xticks([0 pi/2 pi 3*pi/2 2*pi]);
xticklabels({'0', '\pi/2', '\pi', '3\pi//2', '2\pi',})
legend(["First Factor", "Second Factor"]);
hold off;

%e
Y16 = round16(Y)*c;
X16 = round16(X);
[zeros16, poles16] = tf2zpk(Y16,X16);

Yall_16 = (conv(round16(y1), round16(x2))+ ...
        conv(round16(y2), round16(x1)))*c2;
Xall_16 = conv(round16(x1),round16(x2));
[zerosAll_16, polesAll_16] = tf2zpk(Yall_16,Xall_16);

figure();
subplot(3,1,1);
zplane(Y*c,X);
title("Pole-Zero plot H(z)");
subplot(3,1,2);
zplane(Y16,X16);
title("Pole-Zero plot H(z) Rounded to 16th ");
subplot(3,1,3);
zplane(Yall_16,Xall_16);
title("Pole-Zero plot of H(z)Allpass Rounded to 16th")

%H(z) Original:   
    %zeros [-0.9968 ; 0.2876 + 0.9594i ; 0.2876 - 0.9594i]
    %poles [.5050 + 0.7124i ; 0.5050 - .7124i ; 0.4954]
%H(z) Rounded 1/16th; 
    %zeros [-1.0000; 0.2941 + 0.9558i ; 0.2941 - 0.9558i]
    %poles [.5000 + 0.7071i ; 0.5000 - .7071i ; 0.5000]
%H(z) Allpass Rounded 1/16th; 
    %zeros [-1.0000; 0.2500 + 0.9682i ; 0.2500 - 0.9682i]
    %poles [.5000 + 0.7071i ; 0.5000 - .7071i ; 0.5000]

%NOTE: The poles are located within the unit circle, so the function is
%stable.
    
%f)
h16 = freqz(Y16,X16);
h16All = freqz(Yall_16,Xall_16);
Hdb16 = 20*log10(abs(h16));
Hdb16_all = 20*log10(abs(h16All));

figure();
hold on
plot(w,Hdb);
plot(w, Hdb16);
plot(w, Hdb16_all);
hold off
title("Magnitude Response of H(z)and 1/16th rounded coefficients");
xlabel("Frequency (rad)");
ylabel("Magnitude (dB)");
xlim([0 2*pi]);
xticks([0 pi/2 pi 3*pi/2 2*pi]);
xticklabels({'0', '\pi/2', '\pi', '3\pi//2', '2\pi',})
legend(["H(z)", "H(z) Rounded to 1/16", "H(z) Allpass Rounded to 1/16"]);

%g
Y4 = round4(Y)*c;
X4 = round4(X);
[zeros4, poles4] = tf2zpk(Y4,X4);

Yall_4 = (conv(round4(y1), round4(x2))+ ...
        conv(round4(y2), round4(x1)))*c2;
Xall_4 = conv(round4(x1),round4(x2));
[zerosAll_4, polesAll_4] = tf2zpk(Yall_4,Xall_4);

%H(z) Original:   
    %zeros [-0.9968 ; 0.2876 + 0.9594i ; 0.2876 - 0.9594i]
    %poles [.5050 + 0.7124i ; 0.5050 - .7124i ; 0.4954]
%H(z) Rounded 1/4th; 
    %zeros [-1.0000; 0.2500 + 0.9682i ; 0.2500 - 0.9682i]
    %poles [.3867 + 0.7339i ; 0.3867 - .7339i ; 0.7267]
%H(z) Allpass Rounded 1/4th; 
    %zeros [-1.0000; 0.2500 + 0.9682i ; 0.2500 - 0.9682i]
    %poles [.5000 + 0.7071i ; 0.5000 - .7071i ; 0.5000]

%NOTE: The poles are located within the unit circle, so the function is
%still stable.

h4 = freqz(Y4,X4);
h4All = freqz(Yall_4,Xall_4);
Hdb4 = 20*log10(abs(h4));
Hdb4_all = 20*log10(abs(h4All));

figure();
plot(w,Hdb);
hold on
plot(w, Hdb4);
plot(w, Hdb4_all);
hold off
title("Magnitude Response of H(z)and 1/4th rounded coefficients");
xlabel("Frequency (rad)");
ylabel("Magnitude (dB)");
xlim([0 2*pi]);
xticks([0 pi/2 pi 3*pi/2 2*pi]);
xticklabels({'0', '\pi/2', '\pi', '3\pi//2', '2\pi',})
legend(["H(z)", "H(z) Rounded to 1/4", "H(z) Allpass Rounded to 1/4"]);


%% Question 2

%a)
Yf = [1, -.3];
Xf = [.3, -1];
Wpass = (pi/2);
Wstop = (pi/2 + 0.2);
Rpass = 2;
Rstop = 30;
[N, x] = ellipord(Wpass/pi, Wstop/pi, Rpass, Rstop);
[B, A] = ellip(N, Rpass, Rstop, Wpass/pi);
figure;
freqz(B, A);
title('Magnitude Response of H(F(z))');
%It is a high pass filter with a stopband edge 0.856 rads and
%passband edge of 0.985 rads.

%%
%b)
n = linspace(0, pi, 20e3);
omega = exp(1j*n);
tmp = polyval(-1*Yf, omega)./polyval(Xf, omega);
Hzfneg = polyval(B, tmp)./polyval(A, tmp);

figure;
subplot(2,1,1);
freqz(B, A);
title('Magnitude Response of H(F(z))');
subplot(2,1,2);
plot(n./pi, 20*log10(abs(Hzfneg)));
title('Magnitude Response of H(-F(z))');
ylabel("Magnitude (dB)");
grid('on');
%This is a low pass filer with stopband edge at 1.074rads and passband at
%0.985 rads

%d)
Yf2 = [1, 0.4];
Xf2 = [.4,  1];
tmp = polyval(conv(Yf, Yf2), omega)./polyval(conv(Xf, Xf2), omega);
Hzf_pd = polyval(B, tmp)./polyval(A, tmp);

figure;
subplot(2, 1, 1);
plot(n/pi, 20*log10(abs(Hzf_pd)));
title('Magnitude Response of H(F(z)): part d');
grid('on');
%this is a bandpass filter whose order is 8.

tmp = polyval(-1*conv(Yf, Yf2), omega)./polyval(conv(Xf, Xf2), omega);
Hzfneg_pd = polyval(B, tmp)./polyval(A, tmp);
subplot(2, 1, 2)
plot(n/pi, 20*log10(abs(Hzfneg_pd)));
grid('on');
title('Magnitude Response of H(-F(z)): part d');
%this is a bandstop filter whose order is 8.


function [y] = round16(x)
	y = round(x.*16)./16;
end

function [y] = round4(x)
    y = round(x.*4)./4;
end
