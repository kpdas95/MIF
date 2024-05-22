%#########################################################
%####                                                 ####
%####      Multivariate Iterative Filtering(MIF)      ####
%####                                                 ####
%#########################################################

% Please cite the following paper if you are using this code or
% part of the code.
%
% [1]  Kritiprasanna Das and Ram Bilas Pachori. "Schizophrenia 
% detection technique using multivariate iterative filtering and
% multichannel EEG signals." Biomedical Signal Processing and 
% Control 67 (2021): 102525.
% [2] Kritiprasanna Das and Ram Bilas Pachori. 
% "Electroencephalogram based motor imagery brain computer 
% interface using multivariate iterative filtering and spatial 
% filtering." IEEE Transactions on Cognitive and Developmental 
% Systems (2022).
% [3] Antonio Cicone, Jingfang Liu, and Haomin Zhou. "Adaptive 
% local iterative filtering for signal decomposition and 
% instantaneous frequency analysis." Applied and Computational
% Harmonic Analysis 41.2 (2016): 384-411.
% 
% For any queries or help please feel free to write an mail to 
% kpdas95@gmail.com. I will be happy to help.
%% Create a random signal of 5 channels and a length of 1000

x = rand(5,1000);
opt=Settings_IF_v1('IF.Xi',2,'IF.alpha','ave','IF.delta',.001,'IF.NIMFs',20);
MIMF = IterFiltMulti(x,opt);

%% Load a multichannel EEG signal
load('Subject00_1_edfm.mat')

%% Decompose using MIF
sig = val(1:5,:); %Select signal from first 10 channels
opt=Settings_IF_v1('IF.Xi',2,'IF.alpha','ave','IF.delta',.001,'IF.NIMFs',20);
MIMF = IterFiltMulti(sig,opt);

%% Plot EEG signal and five MIMFs
t = linspace(0,1,Fs);

subplot(5,1,1);
plot(t,sig(:,1:Fs)');
ylabel('EEG')

% MIMF 1
subplot(5,1,2)
plot(t, MIMF{1,1}(:,1:Fs)')
ylabel('MIMF_1')

% MIMF 2
subplot(5,1,3)
plot(t, MIMF{1,2}(:,1:Fs)')
ylabel('MIMF_2')

% MIMF 3
subplot(5,1,4)
plot(t, MIMF{1,3}(:,1:Fs)')
ylabel('MIMF_3')

% MIMF 4
subplot(5,1,5)
plot(t, MIMF{1,4}(:,1:Fs)')
ylabel('MIMF_4')




