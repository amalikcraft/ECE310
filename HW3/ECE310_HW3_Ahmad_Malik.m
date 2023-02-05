%Ahmad Malik
%ECE 310
%HW 3

clc
clear
close all

%% Question 4 Part A
n = 1e4;
N = 1;
wname = ['db', int2str(N)];
[h0_1,h1_1,f0_1,f1_1] = wfilters(wname);

%Plotting Fiters when N = 1
[h0, w_h0] = freqz(h0_1, 1, n);
[h1, w_h1] = freqz(h1_1, 1, n);
[f0, w_f0] = freqz(f0_1, 1, n);
[f1, w_f1] = freqz(f1_1, 1, n);

figure;
tcl = tiledlayout(2,2);
nexttile
plot(w_h0, abs(h0));
title("H0")
xlabel("Frequency (rad)")
ylabel("Magntiude")
nexttile
plot(w_h1, abs(h1));
title("H1")
xlabel("Frequency (rad)")
ylabel("Magntiude")
nexttile
plot(w_f0, abs(f0));
title("F0")
xlabel("Frequency (rad)")
ylabel("Magntiude")
nexttile
plot(w_f1, abs(f1));
title("F1")
xlabel("Frequency (rad)")
ylabel("Magntiude")
title(tcl,'Magnitude Responses when N = 1')

%Paraunity Check
para_h0 = conv(h0_1,fliplr(h0_1));
para_h1 = conv(h1_1, fliplr(h1_1));
para_f0 = conv(f0_1, fliplr(f0_1));
para_f1 = conv(f1_1, fliplr(f1_1));

%% Part B #1 and 2
N_4 = 4;
wname_4 = ['db', int2str(N_4)];
[h0_4,h1_4,f0_4,f1_4] = wfilters(wname_4);

[m0_4, w_h0_4] = freqz(h0_4, 1, n);
[m1_4, w_h1_4] = freqz(h1_4, 1, n);
figure
hold on
plot(w_h0_4, abs(m0_4))
plot(w_h1_4, abs(m1_4))
title("Superimposed Magnitude Responses when N = 4")
xlabel("Frequency (rad)")
ylabel("Magntiude")
legend("H0(w)","H1(w)")
xlim([0 pi]);
hold off

N_8 = 8;
wname_8 = ['db', int2str(N_8)];
[h0_8,h1_8,f0_8,f1_8] = wfilters(wname_8);

[m0_8, w_h0_8] = freqz(h0_8, 1, n);
[m1_8, w_h1_8] = freqz(h1_8, 1, n);

figure
hold on
plot(w_h0_8, abs(m0_8))
plot(w_h1_8, abs(m1_8))
title("Superimposed Magnitude Responses when N = 8")
xlabel("Frequency (rad)")
ylabel("Magntiude")
legend("H0(w)","H1(w)")
xlim([0 pi]);
hold off

%% Part B #3
maxdiff4 = max((2 - abs(m0_4.^2) - abs(m1_4.^2)));
maxdiff8 = max((2 - abs(m0_8.^2) - abs(m1_8.^2)));

%% Part B #4-6
% N = 4
E004 = h0_4(1:2:length(h0_4));
E104 = h0_4(2:2:length(h0_4));
E014 = h1_4(1:2:length(h1_4));
E114 = h1_4(2:2:length(h1_4));

R004 = f0_4(2:2:length(f0_4));
R104 = f0_4(1:2:length(f0_4));
R014 = f1_4(2:2:length(f1_4));
R114 = f1_4(1:2:length(f1_4));

P004 = conv(E004, R004) + conv(E014, R014);
P014 = conv(E104, R004) + conv(E114, R014);
P104 = conv(E104, R004) + conv(E114, R014);
P114 = conv(E104, R104) + conv(E114, R114);

% N = 8
E008 = h0_8(1:2:length(h0_8));
E108 = h0_8(2:2:length(h0_8));
E018 = h1_8(1:2:length(h1_8));
E118 = h1_8(2:2:length(h1_8));

R008 = f0_8(2:2:length(f0_8));
R108 = f0_8(1:2:length(f0_8));
R018 = f1_8(2:2:length(f1_8));
R118 = f1_8(1:2:length(f1_8));

P008 = conv(E008, R008) + conv(E018, R018);
P018 = conv(E108, R008) + conv(E118, R018);
P108 = conv(E108, R008) + conv(E118, R018);
P118 = conv(E108, R108) + conv(E118, R118);

%% Part B #7
% N = 4
T_4 = (1/2) * (conv(f0_4, h0_4) + conv(f1_4, h1_4));
fliph0_4 = fliplr(h0_4);
fliph1_4 = fliplr(h1_4);
A_4 = (1/2) * (conv(f0_4, fliph0_4) + conv(f1_4, fliph1_4));

% N = 8
T_8 = (1/2) * (conv(f0_8, h0_8) + conv(f1_8, h1_8));
fliph0_8 = fliplr(h0_8);
fliph1_8 = fliplr(h1_8);
A_8 = (1/2) * (conv(f0_8, fliph0_8) + conv(f1_8, fliph1_8));

%% Part B #8
% N = 4
deriv1_4 = diff(abs(m0_4).^2) / (w_h0_4(2) - w_h0_4(1));
deriv2_4 = (diff(diff(abs(m0_4).^2))/(w_h0_4(2) - w_h0_4(1))^2);
figure
hold on
plot(w_h0_4, [deriv1_4; 0])
plot(w_h0_4, [deriv2_4; 0; 0])
title("First and Second Derivatives when N = 4")
xlabel("Frequency (rad)")
ylabel("Magntiude")
legend('First Derivative','Second Derivative')
xlim([0 pi])

% N = 8
deriv1_8 = diff(abs(m0_8).^2) / (w_h0_8(2)- w_h0_8(1));
deriv2_8 = (diff(diff(abs(m0_8).^2))/(w_h0_8(2)-w_h0_8(1))^2);
figure
hold on
plot(w_h0_8, [deriv1_8; 0])
plot(w_h0_8, [deriv2_8; 0; 0])
title("First and Second Derivatives when N = 8")
xlabel("Frequency (rad)")
ylabel("Magntiude")
legend('First Derivative','Second Derivative')
xlim([0 pi])

%% Part B #9 Two Level
G04 = conv(h0_4, upsample(h0_4, 2));
G14 = conv(h0_4, upsample(h1_4, 2));
G24 = conv(h1_4, upsample(h0_4, 2));
G34 = conv(h1_4, upsample(h1_4, 2 ));

[M_G0_4, W_G0_4] = freqz(G04, 1, n);
[M_G1_4, W_G1_4] = freqz(G14, 1, n);
[M_G2_4, W_G2_4] = freqz(G24, 1, n);
[M_G3_4, W_G3_4] = freqz(G34, 1, n);

figure;
hold on
plot(W_G0_4, abs(M_G0_4));
plot(W_G1_4, abs(M_G1_4));
plot(W_G2_4, abs(M_G2_4));
plot(W_G3_4, abs(M_G3_4));
hold off
title('Two Level Tree Structure Filter Bank when N = 4');
xlabel("Frequency (rad)")
ylabel("Magntiude")
legend('G1','G2', 'G3', 'G4');

%%
G08 = conv(h0_4, conv(upsample(h0_4, 2), upsample(h0_4, 4)));
G18 = conv(h0_4, conv(upsample(h0_4, 2), upsample(h1_4, 4)));
G28 = conv(h0_4, conv(upsample(h1_4, 2), upsample(h0_4, 4)));
G38 = conv(h0_4, conv(upsample(h1_4, 2), upsample(h1_4, 4)));
G48 = conv(h1_4, conv(upsample(h0_4, 2), upsample(h0_4, 4)));
G58 = conv(h1_4, conv(upsample(h0_4, 2), upsample(h1_4, 4)));
G68 = conv(h1_4, conv(upsample(h1_4, 2), upsample(h0_4, 4)));
G78 = conv(h1_4, conv(upsample(h1_4, 2), upsample(h1_4, 4)));


[M_G0_8, W_G0_8] = freqz(G08, 1, n);
[M_G1_8, W_G1_8] = freqz(G18, 1, n);
[M_G2_8, W_G2_8] = freqz(G28, 1, n);
[M_G3_8, W_G3_8] = freqz(G38, 1, n);
[M_G4_8, W_G4_8] = freqz(G48, 1, n);
[M_G5_8, W_G5_8] = freqz(G58, 1, n);
[M_G6_8, W_G6_8] = freqz(G68, 1, n);
[M_G7_8, W_G7_8] = freqz(G78, 1, n);

figure;
hold on
plot(W_G0_8, abs(M_G0_8));
plot(W_G1_8, abs(M_G1_8));
plot(W_G2_8, abs(M_G2_8));
plot(W_G3_8, abs(M_G3_8));
plot(W_G4_8, abs(M_G4_8));
plot(W_G5_8, abs(M_G5_8));
plot(W_G6_8, abs(M_G6_8));
plot(W_G7_8, abs(M_G7_8));
hold off
title('Three Level Tree Structure Filter Bank when N =  8');
xlabel("Frequency (rad)")
ylabel("Magntiude")
legend('G1','G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8');
