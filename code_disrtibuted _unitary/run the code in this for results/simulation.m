% Simulation Results
%% AF 2 relays comparision between the partially coherent and non-coherent
%% receivers and varing the number of symbols L=8 and L=257
N=10^(5);
M=3*10^5
x=0:30;
snr_db=[0:30];
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\Af_rayleigh2_partial.mat');
a=error_af2;
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\Af_rayleigh2_non.mat');
b=error_af1;
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\Af_partial_L=8.mat');
c=error_af2;
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\Af_non_L=8.mat');
d=error_af1;
figure;
semilogy(snr_db,(a(1:31))/N,'-rx',snr_db,(b(1:31))/N,'-.b+',snr_db,(c(1:31))/M,'-.bs',snr_db,(d(1:31))/M,'-.bd');
grid on;
h = legend('Partially coherent L=257','Non coherent L=257','Partially coherent L=8','Non coherent L=8',4);
set(h,'Interpreter','none');
axis([0 30 10^-4 0.5]);
xlabel('SNR(dB)');
ylabel('SER');
%% AF 2 Relay system under differnt channel conditions rayleigh, correlated
%% rayleigh with gamma = 0.3 and gamma = 0.9, ricean channel with k=2 and
%% k=3 L=257
N=10^(5);
x=0:30;
snr_db=[0:30];
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\Af_rayleigh2_non.mat');
a=error_af1;
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\Af_corr_rayleigh_three.mat');
c=error_af1;
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\Af_corr_rayleigh.mat');
d=error_af1;
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\Af_rice.mat');
e=error_af1;
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\Af_rice_3.mat');
f=error_af1;
figure;
semilogy(snr_db,(a(1:31))/N,'-rx',snr_db,(c(1:31))/N,'-.bo',snr_db,(d(1:31))/N,'-.b^',snr_db,(e(1:31))/N,'-.bv',snr_db,(f(1:31))/N,'-.bs');
grid on;
h = legend('Rayleigh','Low correlated Rayleigh ? =0.3','High correlated Rayleigh ? =0.9','Rice k=2','Rice k=3',5);
set(h,'Interpreter','none');
axis([0 30 10^-4 0.5]);
xlabel('SNR(dB)');
ylabel('SER');
%% AF 2 relay system under different power allocation schemes   l=267
N=10^(5);
x=0:30;
snr_db=[0:30];
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\Af_rayleigh2p.mat');
a=error_af1;
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\Af_optimal2p10.mat');
b=error_af1;
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\Af_2p_10_suboptimalpower.mat');
c=error_af1;
figure;
semilogy(snr_db,(a(1:31))/N,'-rx',snr_db,(b(1:31))/N,'-.bo',snr_db,(c(1:31))/N,'-.b^');
grid on;
h = legend('Rayleigh ?_F =1','optimal ?_F =10','suboptimal ?_F =10',3);
set(h,'Interpreter','none');
axis([0 30 10^-4 0.5]);
xlabel('P(dB)');
ylabel('SER');
%% DF and AF performance comparision 
N=10^(5);
M=3*10^(5);
x=0:30;
snr_db=[0:30];
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\Af_non_L=8.mat');
a=error_af1;
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\Af_L=8_2p_optimal.mat');
b=error_af1;
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\AF_p2_10_L=8.mat');
c=error_af1;
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\DF_unitary_diff10.mat');
d=error_af1;
figure;
semilogy(snr_db,(a(1:31))/M,'-rx',snr_db,(b(1:31))/M,'-.bo',snr_db,(c(1:31))/M,'-.b^',snr_db,(d(1:31))/M,'-.bd');
grid on;
h = legend('AF-Rayleigh ?_F =1','AF-optimal ?_F =10','AF-suboptimal ?_F =10','DF ?_F =10',4);
set(h,'Interpreter','none');
axis([0 30 10^-4 0.5]);
xlabel('P(dB)');
ylabel('SER');
%% DF comparision between different encoding scheme
M=3*10^(5);
x=0:30;
snr_db=[0:30];
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\DF_unitary_rayleigh_1.mat');
a=error_af1;
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\DF_unitary_diff.mat');
b=error_af1;
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\Af_non_L=8.mat');
c=error_af1;
figure;
semilogy(snr_db,(a(1:31))/M,'-bd',snr_db,(b(1:31))/M,'-.bo',snr_db,(c(1:31))/M,'-.r^');
grid on;
h = legend('DF-Rayleigh scheme1','DF-Rayleigh scheme2','AF',3);
set(h,'Interpreter','none');
axis([0 30 10^-4 0.5]);
xlabel('P(dB)');
ylabel('SER');
%%DF 2-Relay system under differnt channel conditions rayleigh, correlated
%% rayleigh with gamma = 0.3 and gamma = 0.9, ricean channel with k=2 and
%% k=3 L=8 Scheme2
N=3*10^(5);
x=0:30;
snr_db=[0:30];
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\DF_unitary_diff.mat');
a=error_af1;
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\DF_corr_rayleigh_three_diff.mat');
c=error_af1;
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\DF_corr_rayleigh_nine_diff.mat');
d=error_af1;
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\DF_rice_two_diff.mat');
e=error_af1;
load('C:\Users\n.abhilash\Desktop\code_disrtibuted _unitary\DF_rice_three_diff.mat');
f=error_af1;
figure;
semilogy(snr_db,(a(1:31))/N,'-rx',snr_db,(c(1:31))/N,'-.bo',snr_db,(d(1:31))/N,'-.b^',snr_db,(e(1:31))/N,'-.bv',snr_db,(f(1:31))/N,'-.bs');
grid on;
h = legend('Rayleigh','Low correlated Rayleigh ? =0.3','High correlated Rayleigh ? =0.9','Rice k=2','Rice k=3',5);
set(h,'Interpreter','none');
axis([0 30 10^-4 0.5]);
xlabel('SNR(dB)');
ylabel('SER');
