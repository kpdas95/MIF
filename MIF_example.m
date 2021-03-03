%% Example of Multivariate Iterative Filtering (MIF) method


data = load('E:\Matlab Code\MIF\Example_EEG_signal.mat', 'x');
% EEG signal of a healthy subject from schizophrenia database
% Downloded from "https://repod.icm.edu.pl/dataset.xhtml?persi
% stentId=doi:10.18150/repod.0107441".
%
% Chang, S., Fang, K., Zhang, K. and Wang, J., 2015. Network-based
% analysis of schizophrenia genome-wide association data to detect
% the joint functional association signals. PLoS One,
% 10(7), p.e0133404. https://doi.org/10.1371/journal.pone.0133404



%% Plot of EEG Signal of first five channels
x = data.x;
subplot(5,2,1:2)
plot(x(:,1));
xlim([0,1000]);
ylabel('Fp2');
title('EEG signal of first five channels (Healthy)');

subplot(5,2,3:4)
plot(x(:,2));
xlim([0,1000]);
ylabel('F8');

subplot(5,2,5:6)
plot(x(:,3));
xlim([0,1000]);
ylabel('T4');

subplot(5,2,7:8)
plot(x(:,4));
xlim([0,1000]);
ylabel('T6');

subplot(5,2,9:10)
plot(x(:,5));
xlim([0,1000]);
ylabel('O2');
xlabel('Samples');


%% Decompose EEG signal using MIF
%  options=Settings_IF_v1('IF.Xi',2,'IF.alpha','ave','IF.delta',.001,'IF.NIMFs',100
%  ); Default settings

multiIMF = MIF1(x);


%% Plot of MIMF (For visibility MIMF of first two channel are shown) 

subplot(5,2,1:2)
plot(multiIMF{1,1}(1,:));
hold on;
plot(multiIMF{1,1}(2,:));
xlim([0,1000]); xticks('');
ylabel('MIMF_1');
title('EEG signal of first five channels (Healthy)');

subplot(5,2,3:4)
plot(multiIMF{1,2}(1,:));
hold on;
plot(multiIMF{1,2}(2,:));
xlim([0,1000]); xticks('');
ylabel('MIMF_2');

subplot(5,2,5:6)
plot(multiIMF{1,3}(1,:));
hold on;
plot(multiIMF{1,3}(2,:));
xlim([0,1000]); xticks('');
ylabel('MIMF_3');

subplot(5,2,7:8)
plot(multiIMF{1,4}(1,:));
hold on;
plot(multiIMF{1,4}(2,:));
xlim([0,1000]); xticks('');
ylabel('MIMF_4');

subplot(5,2,9:10)
plot(multiIMF{1,5}(1,:));
hold on;
plot(multiIMF{1,5}(2,:));
xlim([0,1000]);
ylabel('O2');
ylabel('MIMF_5');
%%


