function [out1,out2]=demod_unitary_noncoherent_DF(y,unitary_array,L)
A1= eye(L);
a=[0 1 2 3 4 5 6 7];%size = L

A2=diag((exp(1i*2*pi*a/L)));


y1= ctranspose(y);
calc= zeros(1,L);

for i=1:size(unitary_array,2)
   
    X=[A1*unitary_array(:,i) A2*unitary_array(:,i)];
    calc(i) = y1*X*ctranspose(X)*y;
end

[~,index]=max(calc);

out1 = unitary_array(:,index);
out2=index;

end