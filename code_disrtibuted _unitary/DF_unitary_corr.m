 clc;
clear; 

N=3*10^5; %number of symbols

%Pt=1;% Total Tranmitted Power
%Pb=Pt/2;% Transmitted power of source
%Pr1=Pt/4; % Transmitted power of Relay
%Pr2=Pt/4;% Transmitted power of the Relay
T=8; %Coherence time
snr_db = [0:40];% SNR in dB
u1=[1 3 7 6 5 0 4 2];
%u1=[1 3 0 7 2 5 6 7];
u2=[1 3 0 7 2 5 6 7];
L=8;
M=N/log2(L);
p1 = [exp(1i*2*(pi/L)*u1)];
theta1=diag(p1);
p2= [exp(1i*2*(pi/L)*u2)];
theta2=diag(p2);
data_BS10 = (rand(1,N)>0.5);
data1=zeros(T,M);
data_RS1_tx=zeros(T,M);
data_RS2_tx=zeros(T,M);
 Rx_data=zeros(T,M);

A1= eye(L);
a=[0 1 2 3 4 5 6 7];

A2=diag(exp(1i*2*pi*a/L));
unitary_array1 =zeros(T,L);
unitary_array2 =zeros(T,L);

for a= 1:L
   unitary1=theta1^(a-1);
   for b=1:T;
       unitary_array1(b,a) = unitary1(b,b).*(1/sqrt(T));
   end
   
end
for a= 1:L
   unitary2=theta2^(a-1);
   for b=1:T;
       unitary_array2(b,a) = unitary2(b,b).*(1/sqrt(T));
   end
   
end
data_BS1=zeros(T,M);

%%Modulation at Source
for c = 1:3:N-2
   % for j=1:T
      data_BS1(:,(c-1)/3+1) = map_3bits_to_unitary(data_BS10(c:c+2),unitary_array1);
    %end
end


    

gamma=0.3;
Rt=sqrtm([1 gamma;gamma 1]);
Rr=sqrtm([1 gamma;gamma 1])
    

%% Channel gains
        Bs_Rs=ctranspose(Rr)*((1/sqrt(2)).*[randn(1,M)+1i*randn(1,M);randn(1,M)+1i*randn(1,M)]);% Channel coefficients for Source-Relay1 Link
        Rs_Ms=Rt*((1/sqrt(2)).*[randn(1,M)+1i*randn(1,M);randn(1,M)+1i*randn(1,M)]);
         channel_BS1_RS1 = Bs_Rs(1,:); 
       
                
        channel_RS1_MS1 = Rs_Ms(1,:);
                
        channel_BS1_RS2 =  Bs_Rs(2,:); 
        
        channel_RS2_MS1 = Rs_Ms(2,:);       
        
 %% AWGN     
        nbr1=(1/sqrt(2))*(randn(T,M)+1i*randn(T,M));%AWGN for at Relay for Soure-Relay1 Link
        nr1m=(1/sqrt(2))*(randn(T,M)+1i*randn(T,M));%AWGN for at Desination  for Relay1-Desination Link
        nbr2=(1/sqrt(2))*(randn(T,M)+1i*randn(T,M));%AWGN for at Desination  for source-relay2 Link
     %   nr2m=(1/sqrt(2))*(randn(T,M)+1i*randn(T,M));%AWGN for at Desination  for relay2-Desination Link
        
        error_af1=zeros(1,length(snr_db));
        error_af2=zeros(1,length(snr_db));
 %%MAJOR LOOP
        
for       i=1:length(snr_db)
                snr_linear(i)=10.^(snr_db(i)/10);   % Linear value of SNR                           
               % No(i)=Pt./snr_linear(i);            % Noise power
                Pb = snr_linear(i)/4 ;
                Pr1 = Pb/2; %optimal power allocation
                Pr2 = Pr1;
   %% AMPLIFICATION FACTOR
               Beta1=sqrt(Pr1*T);
               Beta2=sqrt(Pr2*T);
               data_BS1_RS1= zeros(T,M);
               data_RS1_MS1= zeros(T,M);
               data_BS1_RS2= zeros(T,M);
               %data_MS1_RS2= zeros(8,N/3);
               data_RS2_MS1= zeros(T,M);
               
    for j = 1:T
   %% RECEIVED DATA AT RELAY1 FROM BS         
               data_BS1_RS1(j,:)= sqrt(Pb*T)*channel_BS1_RS1.*data_BS1(j,:)+ nbr1(j,:);
              % data_BS1_RS1(2,:)=(channel_BS1_RS1*sqrt(Pb*T)).*data_BS1(2,:)+nbr1(2,:);
   %% RECEIVED DATA AT RELAY2 FROM BS         
               data_BS1_RS2(j,:)=sqrt(Pb*T).*(data_BS1(j,:)).*channel_BS1_RS2+nbr2(j,:);
               %data_BS1_RS2(2,:)=sqrt(Pb*T).*(-data_BS1(2,:)).*channel_BS1_RS2+nbr2(2,:);              
   
              
             
    end
     for k=1:M
                data_RS1_tx(:,k) = demod_mod_DF(data_BS1_RS1(:,k),unitary_array1,unitary_array2,T);
                data_RS2_tx(:,k) = demod_mod_DF(data_BS1_RS2(:,k),unitary_array1,unitary_array2,T);
     end
    for j = 1:T
    %% RECEIVED DATA AT MS FROM RS1     
               data_RS1_MS1(j,:)=Beta1.*( data_RS1_tx(j,:)).*channel_RS1_MS1+nr1m(j,:);
              % data_RS1_MS1(2,:)=Beta1.*data_BS1_RS1(2,:).*channel_RS1_MS1+nr1m(2,:);
              
 
    %% RECEIVED DATA AT MS FROM RS2     
               data_RS2_MS1(j,:) = Beta2.*(exp(1i*2*pi*(j-1)/8)).*data_RS2_tx(j,:).*channel_RS2_MS1;
                %data_RS2_MS1(2,:) = Beta2.*data_BS1_RS2(2,:).*channel_RS2_MS1+nr2m(2,:);
   
                
    end
     %% TOTAL RECEIVED DATA AT MS
    data_MS = data_RS2_MS1 + data_RS1_MS1 ; 
    %% DEMODULATION NON COHERENT SUB OPTIMAL DETECTOR
            for k=1:M
                [data1(:,k),Symbol_index(k)] = demod_unitary_noncoherent_DF(data_MS(:,k),unitary_array2,T);
                Rx_data(:,k)=unitary_array1(:,Symbol_index(k));
                %% CALCULATING ERRORS
              if abs(Rx_data(:,k)- data_BS1(:,k)) == 0
                  error_af1(i)= error_af1(i);
              else
                  error_af1(i)= error_af1(i)+ 1;
             end
                                  
            end
           
    
       
   
        
        
    
end 

save('DF_corr_rayleigh_three_diff.mat','error_af1');
%% PLOTTING THE SIMULATION AND THEORATICAL RESULTS
    figure 
    %% SIMULATION RESULTS
    simber_af1=error_af1/N;
    simber_af2=error_af2/N;
    semilogy(snr_db,simber_af1,'g.-')
     hold on   
    semilogy(snr_db,simber_af2,'b.-')
    grid on
    axis([0 40 10^-5 0.5]);
    xlabel('SNR(dB)');
    ylabel('SER');