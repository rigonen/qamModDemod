clc; clear;

MinSNR = 0; 
MaxSNR = 10;
step = 1;
SNR = MinSNR:step:MaxSNR;
BER1 = [0.079974, 0.056849, 0.0369792, 0.022474, 0.0124219, 0.005625, ...
0.00242187, 0.000833333, 0.000182292, 0, 0];
BER2 = [0.143411, 0.117865, 0.0972135, 0.0764583, 0.0582031,...
    0.0419792, 0.0273958, 0.0178385, 0.0090625, 0.00466146, 0.00161458];
BER3 = [0.197943, 0.176797, 0.156172, 0.13651, 0.116719, ...
    0.101927, 0.0847396, 0.0683594, 0.0530208, 0.0376302, 0.0261979];
BER4 = [0.25474, 0.237057, 0.217005, 0.197005, 0.183932,...
    0.159766, 0.14612, 0.128203, 0.114818, 0.0984635, 0.0813281];




[BER1_ref ~]=berawgn(SNR, 'qam', 4);
[BER2_ref ~]=berawgn(SNR, 'qam', 16);
[BER3_ref ~]=berawgn(SNR, 'qam', 64);
[BER4_ref ~]=berawgn(SNR, 'qam', 256);

figure(1);

semilogy(SNR, BER1, 'k:');
hold on;
semilogy(SNR, BER1_ref, 'k^');
hold on;
semilogy(SNR, BER2, 'r:');
hold on;
semilogy(SNR, BER2_ref, 'r^');
hold on;
semilogy(SNR, BER3, 'b:');
hold on; 
semilogy(SNR, BER3_ref, 'b^');
hold on;
semilogy(SNR, BER4, 'm:');
hold on; 
semilogy(SNR, BER4_ref, 'm^');
hold on; grid on;
% xlim([0 10]);
ylim([10e-5 1])
title('Семейство кривых помехоустойчивости QAM')
legend('QAM4', 'QAM4 теория', 'QAM16', 'QAM16 теория',...
    'QAM64','QAM64 теория', 'QAM256','QAM256 теория')
xlabel('E_b/N_0, дБ');
ylabel('BER');