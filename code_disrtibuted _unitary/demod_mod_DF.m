function out = demod_mod_DF(rx,unitary_array_demod,unitary_array_mod,T)
%% demodulation part
T=8;
L=size(unitary_array_demod,2);
calc= zeros(1,L);
for i=1:size(unitary_array_demod,2)
   
    X=unitary_array_demod(:,i);
    calc(i) = ctranspose(rx)*X*ctranspose(X)*rx;
end

[~,index]=max(calc);

%%modulation
out= unitary_array_mod(:,index);


