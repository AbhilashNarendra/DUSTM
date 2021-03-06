 clc;
clear; 

N=3*10^2; %number of symbols
L=8;
M=N/log2(L);
%Pt=1;% Total Tranmitted Power
%Pb=Pt/2;% Transmitted power of source
%Pr1=Pt/4; % Transmitted power of Relay
%Pr2=Pt/4;% Transmitted power of the Relay
T=2; %Coherence time
snr_db = [0:40];% SNR in dB
u=[1 5];

p = [exp(1i*2*(pi/L)*u)];
theta=diag(p);
data_BS10 = (rand(1,N)>0.5);
A1= [1 0;0 1];


A2=[-1 0;0 1];

unitary_array =zeros(T,L);
for a= 1:L
   unitary = theta^(a-1);
   for b=1:T;
       unitary_array(b,a) = unitary(b,b)*(1/sqrt(2));
   end
   
end


       
        error_af1=zeros(1,length(snr_db));
        error_af2=zeros(1,length(snr_db));
 %%MAJOR LOOP
        frame_size = 10;
for       i=1:length(snr_db)
    
                for frame=1:frame_size:M-frame_size+1
                    
                    
                   data_BS1(:,frame:frame+frame_size-1) = modulate_diff(data_BS10(3*frame-2:(3*(frame+frame_size-1))),unitary_array); 
                   
                   channel_BS1_RS1=(1/sqrt(2)).*[randn(1,1)+1i*randn(1,1)];% Channel coefficients for Source-Relay1 Link
        a11=sum((abs(channel_BS1_RS1)))/M; %Average channel gain for Soure-Relay Link Channel
                
        channel_RS1_MS1=(1/sqrt(2)).*[randn(1,1)+1i*randn(1,1)];% Channel coefficients for Relay1-Desination Link
        b11=sum((abs(channel_RS1_MS1)))/M; %Average channel gain for Relay-Desination Link Channel
                
        channel_BS1_RS2=(1/sqrt(2))*[randn(1,1)+1i*randn(1,1)];% Channel coefficients for Source-Relay2 Link
        a12=sum((abs(channel_BS1_RS2)))/M; %Average channel gain for Soure-Desination Link Channel
        
        channel_RS2_MS1=(1/sqrt(2)).*[randn(1,1)+1i*randn(1,1)];% Channel coefficients for Relay2-Desination Link
        b21=sum((abs(channel_RS2_MS1)))/M; %Average channel gain for Relay-Desination Link Channel
        
        nbr1=(1/sqrt(2))*(randn(T,frame_size)+1i*randn(T,frame_size));%AWGN for at Relay for Soure-Relay1 Link
        nr1m=(1/sqrt(2))*(randn(T,frame_size)+1i*randn(T,frame_size));%AWGN for at Desination  for Relay1-Desination Link
        nbr2=(1/sqrt(2))*(randn(T,frame_size)+1i*randn(T,frame_size));%AWGN for at Desination  for source-relay2 Link
        nr2m=(1/sqrt(2))*(randn(T,frame_size)+1i*randn(T,frame_size));%AWGN for at Desination  for relay2-Desination Link
        
        
                snr_linear(i)=10.^(snr_db(i)/10);   % Linear value of SNR                           
               % No(i)=Pt./snr_linear(i);            % Noise power
                Pb = snr_linear(i)/4 + sqrt(snr_linear(i)*snr_linear(i) + 4*snr_linear(i))/4;
                Pr1 = Pb/2; %optimal power allocation
                Pr2 = Pr1;
   %% AMPLIFICATION FACTOR
               Beta1=sqrt(Pr1./(Pb+1));
               Beta2=sqrt(Pr2./(Pb+1));
               data_BS1_RS1= zeros(2,frame_size);
               data_RS1_MS1= zeros(2,frame_size);
               data_BS1_RS2= zeros(2,frame_size);
               %data_MS1_RS2= zeros(8,N/3);
               data_RS2_MS1= zeros(2,frame_size);
               
    for j = 1:T
   %% RECEIVED DATA AT RELAY1 FROM BS         
               data_BS1_RS1(j,frame:frame+frame_size-1)= sqrt(Pb*T)*channel_BS1_RS1.*data_BS1(j,frame:frame+frame_size-1)+ nbr1(j,:);
              % data_BS1_RS1(2,:)=(channel_BS1_RS1*sqrt(Pb*T)).*data_BS1(2,:)+nbr1(2,:);
                            
   %% RECEIVED DATA AT MS FROM RS1     
               data_RS1_MS1(j,frame:frame+frame_size-1)=Beta1.*(data_BS1_RS1(j,frame:frame+frame_size-1)).*channel_RS1_MS1+nr1m(j,:);
              % data_RS1_MS1(2,:)=Beta1.*data_BS1_RS1(2,:).*channel_RS1_MS1+nr1m(2,:);
        
  %% RECEIVED DATA AT RELAY2 FROM BS         
               data_BS1_RS2(j,frame:frame+frame_size-1)=sqrt(Pb*T).*(data_BS1(j,frame:frame+frame_size-1)).*channel_BS1_RS2+nbr2(j,:);
               %data_BS1_RS2(2,:)=sqrt(Pb*T).*(-data_BS1(2,:)).*channel_BS1_RS2+nbr2(2,:);
    %% RECEIVED DATA AT MS FROM RS2     
               data_RS2_MS1(j,frame:frame+frame_size-1) = Beta2.*A2(j,j).*data_BS1_RS2(j,frame:frame+frame_size-1).*channel_RS2_MS1 ;
                %data_RS2_MS1(2,:) = Beta2.*data_BS1_RS2(2,:).*channel_RS2_MS1+nr2m(2,:);
   
                
    end
     %% TOTAL RECEIVED DATA AT MS
    data_MS = data_RS2_MS1 + data_RS1_MS1 ; 
  
    dummy1 = sqrt(Pb*Pr1*T/(Pb+1))*[1 -1;1 1]*[1;0];%*[channel_BS1_RS1*channel_RS1_MS1;channel_BS1_RS2*channel_RS2_MS1];
    %% DEMODULATION NON COHERENT SUB OPTIMAL DETECTOR
            for k=1:frame_size-1
                data1(:,frame)   = ddstm_demod2(dummy1,data_MS(:,frame),unitary_array,L);
                data1(:,frame+k) = ddstm_demod2(data_MS(:,frame+k-1),data_MS(:,frame+k),unitary_array,L);
                
            end
            %data10 = zeros(1,N);
            
             
               end
             for v=0:length(data1)-1  
                 if de2bi(data1(v+1))== 0
                 data10(2*v+1:1:2*v+3) = [0 0 0];
                 elseif de2bi(data1(v+1))== 1
                 data10(2*v+1:1:2*v+3) = [0 0 1];
                 elseif de2bi(data1(v+1))== 2
                 data10(2*v+1:1:2*v+3) = [0 1 0];
                 elseif de2bi(data1(v+1))== 3
                 data10(2*v+1:1:2*v+3) = [0 1 1];
                 elseif de2bi(data1(v+1))== 4
                 data10(2*v+1:1:2*v+3) = [1 0 0];
                  elseif de2bi(data1(v+1))== 5
                 data10(2*v+1:1:2*v+3) = [1 0 1];
                  elseif de2bi(data1(v+1))== 6
                 data10(2*v+1:1:2*v+3) = [1 1 0];
                  elseif de2bi(data1(v+1))== 7
                 data10(2*v+1:1:2*v+3) = [1 1 1];
                 
                 end
             end
             
             
             error_af1(i)= size(find(data10-data_BS10(1:length(data10))),2);       
    
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