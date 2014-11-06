clc;
clear; 

N=10^5; %number of symbols
M=N/8;
%Pt=1;% Total Tranmitted Power
%Pb=Pt/2;% Transmitted power of source
%Pr1=Pt/4; % Transmitted power of Relay
%Pr2=Pt/4;% Transmitted power of the Relay
T=8; %Coherence time
p_db = [0:40];% SNR in dB
u=[1 7 60 79 187 125 198 154];
L=257;
p = [exp(1i*2*(pi/L)*u)];
theta=diag(p);
data_BS10 = (rand(1,N)>0.5);%% Data generation so that s*s  =1
data1=zeros(T,N/8);
data2=zeros(T,N/8);% Received data matrix
A1= eye(L);
a=[0 1 2 3 4 5 6 7];%size = L

A2=diag(exp(1i*2*pi*a/L));
%unitary_array = zeros(8,8);
%u=[7 60 79 187 125 198 154];
%theta=diag(exp(1i*2*pi/L*u))/sqrt(8);
% unitary_array=[s1 s2 s3 s4 s5 s6 s7 s8];
unitary_array =zeros(T,L);
%L=8;%number of symbols
for a= 1:L
   unitary=theta^(a-1);
   for b=1:T;
       unitary_array(b,a) = unitary(b,b).*(1/sqrt(T));
   end
   
end
data_BS1=zeros(T,M);
for c = 1:8:N-7
   % for j=1:T
      data_BS1(:,(c-1)/8+1) = map_3bits_to_unitary(data_BS10(c:c+7),unitary_array);
    %end
end


    

%% Channel gains
        channel_BS1_RS1=(sqrt(10)/sqrt(2)).*[randn(1,M)+1i*randn(1,M)];% Channel coefficients for Source-Relay1 Link
        a11=sum((abs(channel_BS1_RS1)))/M; %Average channel gain for Soure-Relay Link Channel
                
        channel_RS1_MS1=(1/sqrt(2)).*[randn(1,M)+1i*randn(1,M)];% Channel coefficients for Relay1-Desination Link
        b11=sum((abs(channel_RS1_MS1)))/M; %Average channel gain for Relay-Desination Link Channel
                
        channel_BS1_RS2=(sqrt(10)/sqrt(2))*[randn(1,M)+1i*randn(1,M)];% Channel coefficients for Source-Relay2 Link
        a12=sum((abs(channel_BS1_RS2)))/M; %Average channel gain for Soure-Desination Link Channel
        
        channel_RS2_MS1=(1/sqrt(2)).*[randn(1,M)+1i*randn(1,M)];% Channel coefficients for Relay2-Desination Link
        b21=sum((abs(channel_RS2_MS1)))/M; %Average channel gain for Relay-Desination Link Channel
        
 %% AWGN     
        nbr1=(1/sqrt(2))*(randn(T,M)+1i*randn(T,M));%AWGN for at Relay for Soure-Relay1 Link
        nr1m=(1/sqrt(2))*(randn(T,M)+1i*randn(T,M));%AWGN for at Desination  for Relay1-Desination Link
        nbr2=(1/sqrt(2))*(randn(T,M)+1i*randn(T,M));%AWGN for at Desination  for source-relay2 Link
        nr2m=(1/sqrt(2))*(randn(T,M)+1i*randn(T,M));%AWGN for at Desination  for relay2-Desination Link
        
        error_af1=zeros(1,length(p_db));
        
 %%MAJOR LOOP
        
for       i=1:length(p_db)
                p_linear(i)=10.^(p_db(i)/10);   % Linear value of SNR                           
               % No(i)=Pt./snr_linear(i);            % Noise power
                Pb = p_linear(i)/2;
                Pr1 = Pb/2; %optimal power allocation
                Pr2 = Pr1;
   %% AMPLIFICATION FACTOR
               Beta1=sqrt(Pr1./(10*Pb+1));
               Beta2=sqrt(Pr2./(10*Pb+1));
               data_BS1_RS1= zeros(8,M);
               data_RS1_MS1= zeros(8,M);
               data_BS1_RS2= zeros(8,M);
               %data_MS1_RS2= zeros(8,N/3);
               data_RS2_MS1= zeros(8,M);
               
    for j = 1:T
   %% RECEIVED DATA AT RELAY1 FROM BS         
               data_BS1_RS1(j,:)= sqrt(Pb*T)*channel_BS1_RS1.*data_BS1(j,:)+ nbr1(j,:);
              % data_BS1_RS1(2,:)=(channel_BS1_RS1*sqrt(Pb*T)).*data_BS1(2,:)+nbr1(2,:);
                            
   %% RECEIVED DATA AT MS FROM RS1     
               data_RS1_MS1(j,:)=Beta1.*(data_BS1_RS1(j,:)).*channel_RS1_MS1+nr1m(j,:);
              % data_RS1_MS1(2,:)=Beta1.*data_BS1_RS1(2,:).*channel_RS1_MS1+nr1m(2,:);
        
  %% RECEIVED DATA AT RELAY2 FROM BS         
               data_BS1_RS2(j,:)=sqrt(Pb*T).*(data_BS1(j,:)).*channel_BS1_RS2+nbr2(j,:);
               %data_BS1_RS2(2,:)=sqrt(Pb*T).*(-data_BS1(2,:)).*channel_BS1_RS2+nbr2(2,:);
    %% RECEIVED DATA AT MS FROM RS2     
               data_RS2_MS1(j,:) = Beta2.*(exp(1i*2*pi*(j-1)/8)).*data_BS1_RS2(j,:).*channel_RS2_MS1+ nr2m(j,:);
                %data_RS2_MS1(2,:) = Beta2.*data_BS1_RS2(2,:).*channel_RS2_MS1+nr2m(2,:);
   
                
    end
     %% TOTAL RECEIVED DATA AT MS
    data_MS = data_RS2_MS1 + data_RS1_MS1 ; 
    %% DEMODULATION NON COHERENT SUB OPTIMAL DETECTOR
            for k=1:M
                data1(:,k) = demod_unitary_noncoherent_2(data_MS(:,k),unitary_array,T);
                %% CALCULATING ERRORS
              if abs(data_BS1(:,k)- data1(:,k)) == 0
                  error_af1(i)= error_af1(i);
              else
                  error_af1(i)= error_af1(i)+ 1;
             end
                                  
            end
           
         %   error_af1(i)=size(col1,2);
   
   
        
        
    
end 
save('Af_2p_10_suboptimalpower.mat','error_af1');

%% PLOTTING THE SIMULATION AND THEORATICAL RESULTS
    figure 
    %% SIMULATION RESULTS
    simber_af1=error_af1/N;
    
    semilogy(p_db,simber_af1,'g.-')
     hold on   
    
    grid on
    axis([0 40 10^-5 0.5]);
    xlabel('p(dB)');
    ylabel('SER');