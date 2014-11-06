clc;
clear; 

N=2*10^2; %number of symbols
L=4;
s=log2(L);
M=N/log2(L);

T=2; %Coherence time
snr_db = [0:40];% SNR in dB

data_BS10 = (rand(1,N)>0.5);
A1= [1 0;0 1];


A2=[0 -1;1 0];
u1= (1/sqrt(T))*[1 -1;1 1];
u2= (1/sqrt(T))*[-1 -1;1 -1];
u3= (1/sqrt(T))*[1 1;-1 1];
u4= (1/sqrt(T))*[-1 1;-1 -1];

unitary_array =[u1 u2 u3 u4 ];



       data_BS1= zeros(2,M);
        error_af1=zeros(1,length(snr_db));
        error_af2=zeros(1,length(snr_db));
        
        for f=1:s:N-s+1
          data_BS1(:,(f-1)/s+1) =modulate_diff_1(data_BS10(:,f:f+s-1),unitary_array,4); 
        end
        
 %%MAJOR LOOP
        
for       i=1:length(snr_db)
    
               
                    
                    
                 
                   
                   channel_BS1_RS1=(1/sqrt(2)).*[randn(1,1)+1i*randn(1,1)];% Channel coefficients for Source-Relay1 Link
        a11=sum((abs(channel_BS1_RS1)))/M; %Average channel gain for Soure-Relay Link Channel
                
        channel_RS1_MS1=(1/sqrt(2)).*[randn(1,1)+1i*randn(1,1)];% Channel coefficients for Relay1-Desination Link
        b11=sum((abs(channel_RS1_MS1)))/M; %Average channel gain for Relay-Desination Link Channel
                
        channel_BS1_RS2=(1/sqrt(2))*[randn(1,1)+1i*randn(1,1)];% Channel coefficients for Source-Relay2 Link
        a12=sum((abs(channel_BS1_RS2)))/M; %Average channel gain for Soure-Desination Link Channel
        
        channel_RS2_MS1=(1/sqrt(2)).*[randn(1,1)+1i*randn(1,1)];% Channel coefficients for Relay2-Desination Link
        b21=sum((abs(channel_RS2_MS1)))/M; %Average channel gain for Relay-Desination Link Channel
        
        nbr1=(1/sqrt(2))*(randn(T,M)+1i*randn(T,M));%AWGN for at Relay for Soure-Relay1 Link
        nr1m=(1/sqrt(2))*(randn(T,M)+1i*randn(T,M));%AWGN for at Desination  for Relay1-Desination Link
        nbr2=(1/sqrt(2))*(randn(T,M)+1i*randn(T,M));%AWGN for at Desination  for source-relay2 Link
        nr2m=(1/sqrt(2))*(randn(T,M)+1i*randn(T,M));%AWGN for at Desination  for relay2-Desination Link
        
        
                snr_linear(i)=10.^(snr_db(i)/10);   % Linear value of SNR                           
               % No(i)=Pt./snr_linear(i);            % Noise power
                Pb = snr_linear(i)/4 + sqrt(snr_linear(i)*snr_linear(i) + 4*snr_linear(i))/4;
                Pr1 = Pb/2; %optimal power allocation
                Pr2 = Pr1;
   %% AMPLIFICATION FACTOR
               Beta1=sqrt(Pr1./(Pb+1));
               Beta2=sqrt(Pr2./(Pb+1));
               data_BS1_RS1= zeros(2,M);
               data_RS1_MS1= zeros(2,M);
               data_BS1_RS2= zeros(2,M);
               %data_MS1_RS2= zeros(8,N/3);
               data_RS2_MS1= zeros(2,M);
               
    for j = 1:T
   %% RECEIVED DATA AT RELAY1 FROM BS         
               data_BS1_RS1(j,:)= sqrt(Pb*T)*channel_BS1_RS1.*data_BS1(j,:)+ nbr1(j,:);
              % data_BS1_RS1(2,:)=(channel_BS1_RS1*sqrt(Pb*T)).*data_BS1(2,:)+nbr1(2,:);
                            
   %% RECEIVED DATA AT MS FROM RS1     
               data_RS1_MS1(j,:)=Beta1.*data_BS1_RS1(j,:).*channel_RS1_MS1+nr1m(j,:);
              % data_RS1_MS1(2,:)=Beta1.*data_BS1_RS1(2,:).*channel_RS1_MS1+nr1m(2,:);
        
  %% RECEIVED DATA AT RELAY2 FROM BS         
               data_BS1_RS2(j,:)=sqrt(Pb*T).*(data_BS1(j,:)).*channel_BS1_RS2+nbr2(j,:);
                       
    end
    %% RECEIVED DATA AT MS FROM RS2     
     data_RS2_MS1(1,:) = -Beta2.*data_BS1_RS2(2,:).*channel_RS2_MS1 ;
      data_RS2_MS1(2,:) = Beta2.*data_BS1_RS2(1,:).*channel_RS2_MS1 ;
     %% TOTAL RECEIVED DATA AT MS
    data_MS = data_RS2_MS1 + data_RS1_MS1 ; 
  
    %dummy1 = sqrt(Pb*Pr1*T/(Pb+1))*[1 -1;1 1]*[channel_BS1_RS1*channel_RS1_MS1;0];
    %% DEMODULATION NON COHERENT SUB OPTIMAL DETECTOR
            for k=1:M-1
                %data1(:,1)   = ddstm_demod2_1(dummy1,data_MS(:,1),unitary_array,L);
                data1(:,k+1) = ddstm_demod2_1(data_MS(:,k),data_MS(:,k+1),unitary_array,L);
                
            end
            %data10 = zeros(1,N);
            
             

             for v=1:length(data1)-1  
                 if data1(v+1)== 0
                 data10(s*v+1:1:s*v+s) = [0 0];
                 elseif data1(v+1)== 1
                 data10(s*v+1:1:s*v+s) = [1 0];
                 elseif data1(v+1)== 2
                 data10(s*v+1:1:s*v+s) = [0 1];
                 elseif data1(v+1)== 3
                 data10(s*v+1:1:s*v+s) = [1 1];
                            
                 end
             end
             
           
             error_af1(i)= size(find(data10(s+1:length(data10))-data_BS10(s+1:length(data10))),2);       
    
    
end
%% PLOTTING THE SIMULATION AND THEORATICAL RESULTS
    figure 
    %% SIMULATION RESULTS
    simber_af1=error_af1/N;
   % simber_af2=error_af2/N;
    semilogy(snr_db,simber_af1,'g.-')
   %  hold on   
   % semilogy(snr_db,simber_af2,'b.-')
    grid on
    axis([0 40 10^-5 0.5]);
    xlabel('SNR(dB)');
    ylabel('SER');