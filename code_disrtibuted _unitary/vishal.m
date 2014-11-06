clc;
clear;
t = 2;
r = 2;
frm_l = 5000;  %no. of bits in each frame
frm_n = 30;   %no. of frames
bits = (rand(1,frm_l*frm_n)>0.5);
bits_est = zeros(1,frm_l*frm_n);
bit_f = zeros(1,frm_l);
snr_db = 0:0.5:20;
snr = 10.^(snr_db./10);
ber = zeros(1,length(snr_db));
%h = zeros(2,2);
g1 = [-1 0;0 -1];
g2 = [1 0;0 1];
d = [1 -1;1 1];
h = randn(2,2)+1j*randn(2,2);
hsd = 1;  %standard deviation of elements of h;
noise = zeros(2,2);
nsd = 1;  %standard deviation of elements of noise;
x = zeros(2,2*frm_l+2);
y = zeros(2,2*frm_l+2);
for i = 1:length(snr)  %for each snr
    
    snr_t = snr(i)/t;  %snr per transmitter
    
    for f = 1:30  %for each of the 30 frames
        %h = complex(normrnd(0,hsd),normrnd(0,hsd));  %h specific to each frame
        %h = randn(2,2)+1j*randn(2,2);
        bit_f = bits((f-1)*frm_l+1:f*frm_l);
        x(:,1:2) = d;
        for b = 1:frm_l    
                if bit_f(b) == 1
                    x(:,2*b+1:2*b+2) = x(:,2*b-1:2*b)*g1;
                elseif bit_f(b) == 0
                    x(:,2*b+1:2*b+2) = x(:,2*b-1:2*b)*g2;
                end
        end
        y = sqrt(snr_t)*h*x+randn(2,2*frm_l+2)+1j*randn(2,2*frm_l+2);
        d = x(:,2*frm_l+1:2*frm_l+2);
        for k = 1:frm_l
            m = [real(trace(g1*ctranspose(y(:,2*k+1:2*k+2))*y(:,2*k-1:2*k))) real(trace(g2*ctranspose(y(:,2*k+1:2*k+2))*y(:,2*k+1:2*k+2)))];
            if m(1)>m(2)
                bits_est(frm_l*(f-1)+k) = 1;
            elseif m(1) <= m(2)
                bits_est(frm_l*(f-1)+k) = 0;
            end
        end
    end
    ber(i) = nnz(bits-bits_est)/length(bits);
end
figure(2)
semilogy(snr_db,ber,'r-s');
grid on;
title('dstm 2x2 (ber vs snr)');