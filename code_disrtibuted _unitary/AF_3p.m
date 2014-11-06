 clc;
clear;

N= 10^4; %number of symbols
%Pt=1;% Total Tranmitted Power
%Pb=Pt/2;% Transmitted power of source
%Pr1=Pt/4; % Transmitted power of Relay
%Pr2=Pt/4;% Transmitted power of the Relay
T=8; %Coherence time
p_db = [0:40];% SNR in dB
u=[220 191 6 87 219 236 173 170];
L=257;
M=N/log2(L-1);
p = [exp(1i*2*(pi/L)*u)];
theta=diag(p);
data_BS10 = (rand(1,N)>0.5);%% Data generation so that s*s  =1
data1=zeros(T,M);
data2=zeros(T,M);% Received data matrix
A1= eye(T);
a1=[0 5 2 7 4 1 6 3];%size = T
a2=[0 6 4 2 0 6 4 2];
A2=diag(exp(1i*2*pi*a1/L));
A3=diag(exp(1i*2*pi*a2/L));
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
for c = 1:log2(L-1):N-log2(L-1)+1
   % for j=1:T
      data_BS1(:,((c-1)/8)+1) = map_3bits_to_unitary(data_BS10(c:c+7),unitary_array);
    %end
end


    

%% Channel gains
        channel_BS1_RS1=(1/sqrt(2)).*[randn(1,M)+1i*randn(1,M)];% Channel coefficients for Source-Relay1 Link
        a11=sum((abs(channel_BS1_RS1)))/M; %Average channel gain for Soure-Relay Link Channel
                
        channel_RS1_MS1=(1/sqrt(2)).*[randn(1,M)+1i*randn(1,M)];% Channel coefficients for Relay1-Desination Link
        b11=sum((abs(channel_RS1_MS1)))/M; %Average channel gain for Relay-Desination Link Channel
                
        channel_BS1_RS2=(1/sqrt(2))*[randn(1,M)+1i*randn(1,M)];% Channel coefficients for Source-Relay2 Link
        a12=sum((abs(channel_BS1_RS2)))/M; %Average channel gain for Soure-Desination Link Channel
        
        channel_RS2_MS1=(1/sqrt(2)).*[randn(1,M)+1i*randn(1,M)];% Channel coefficients for Relay2-Desination Link
        b21=sum((abs(channel_RS2_MS1)))/M; %Average channel gain for Relay-Desination Link Channel
       
        channel_RS3_MS1=(1/sqrt(2)).*[randn(1,M)+1i*randn(1,M)];
        
        channel_BS1_RS3=(1/sqrt(2))*[randn(1,M)+1i*randn(1,M)];
 %% AWGN     
        nbr1=(1/sqrt(2))*(randn(T,M)+1i*randn(T,M));%AWGN for at Relay for Soure-Relay1 Link
        nr1m=(1/sqrt(2))*(randn(T,M)+1i*randn(T,M));%AWGN for at Desination  for Relay1-Desination Link
        nbr2=(1/sqrt(2))*(randn(T,M)+1i*randn(T,M));%AWGN for at Desination  for source-relay2 Link
        nr2m=(1/sqrt(2))*(randn(T,M)+1i*randn(T,M));%AWGN for at Desination  for relay2-Desination Link
        nbr3=(1/sqrt(2))*(randn(T,M)+1i*randn(T,M));
        error_af1=zeros(1,length(p_db));
        %error_af2=zeros(1,length(snr_db));
 %%MAJOR LOOP
        
for       i=1:length(p_db)
                p_linear(i)=10.^(p_db(i)/10);   % Linear value of SNR                           
               % No(i)=Pt./snr_linear(i);            % Noise power
                Pb =  p_linear(i);
                Pr1 = Pb/3; %optimal power allocation
                Pr2 = Pr1;
                Pr3 = Pr1;
   %% AMPLIFICATION FACTOR
               Beta1=sqrt(Pr1./(Pb+1));
               Beta2=sqrt(Pr2./(Pb+1));
               Beta3=sqrt(Pr2./(Pb+1));
               data_BS1_RS1= zeros(8,M);
               data_RS1_MS1= zeros(8,M);
               data_BS1_RS2= zeros(8,M);
               %data_MS1_RS2= zeros(8,N/3);
               data_RS2_MS1= zeros(8,M);
               data_RS3_MS1= zeros(8,M);
               data_BS1_RS3= zeros(8,M);
    for j = 1:T
   %% RECEIVED DATA AT RELAY1 FROM BS         
               data_BS1_RS1(j,:)= sqrt(Pb*T)*channel_BS1_RS1.*data_BS1(j,:)+ nbr1(j,:);
              % data_BS1_RS1(2,:)=(channel_BS1_RS1*sqrt(Pb*T)).*data_BS1(2,:)+nbr1(2,:);
                            
   %% RECEIVED DATA AT MS FROM RS1     
               data_RS1_MS1(j,:)=Beta1.*(data_BS1_RS1(j,:)).*channel_RS1_MS1;
              % data_RS1_MS1(2,:)=Beta1.*data_BS1_RS1(2,:).*channel_RS1_MS1+nr1m(2,:);
        
  %% RECEIVED DATA AT RELAY2 FROM BS         
               data_BS1_RS2(j,:)=sqrt(Pb*T).*(data_BS1(j,:)).*channel_BS1_RS2+nbr2(j,:);
               %data_BS1_RS2(2,:)=sqrt(Pb*T).*(-data_BS1(2,:)).*channel_BS1_RS2+nbr2(2,:);
    %% RECEIVED DATA AT MS FROM RS2     
               data_RS2_MS1(j,:) = Beta2.*(A2(j,j)).*data_BS1_RS2(j,:).*channel_RS2_MS1+ nr2m(j,:);
                %data_RS2_MS1(2,:) = Beta2.*data_BS1_RS2(2,:).*channel_RS2_MS1+nr2m(2,:);
     %% RECEIVED DATA AT RELAY3 FROM BS         
               data_BS1_RS3(j,:)=sqrt(Pb*T).*(data_BS1(j,:)).*channel_BS1_RS3+nbr3(j,:);
               %data_BS1_RS2(2,:)=sqrt(Pb*T).*(-data_BS1(2,:)).*channel_BS1_RS2+nbr2(2,:);
    %% RECEIVED DATA AT MS FROM RS3     
               data_RS3_MS1(j,:) = Beta3.*(A3(j,j)).*data_BS1_RS3(j,:).*channel_RS3_MS1;
   
                
    end
     %% TOTAL RECEIVED DATA AT MS
    data_MS = data_RS2_MS1 + data_RS1_MS1+data_RS3_MS1 ; 
    %% DEMODULATION NON COHERENT SUB OPTIMAL DETECTOR
            for k=1:M
                data1(:,k) = demod_unitary_noncoherent_3(data_MS(:,k),unitary_array,T);
                %% CALCULATING ERRORS
              if abs(data_BS1(:,k)- data1(:,k)) == 0
                  error_af1(i)= error_af1(i);
              else
                  error_af1(i)= error_af1(i)+ 1;
             end
                                  
            end
           
    
   
        
        
    
end 
save('Af_rayleigh3_power4.mat','error_af1');

%% PLOTTING THE SIMULATION AND THEORATICAL RESULTS
    figure 
    %% SIMULATION RESULTS
    simber_af1=error_af1/N;
    
    semilogy(p_db,simber_af1,'g.-')
     hold on   
   
    grid on
    axis([0 40 10^-5 0.5]);
    xlabel('P(dB)');
    ylabel('SER');