function out=demod_unitary_noncoherent(y,unitary_array,L)
A1= eye(L);
a=[0 1 2 3 4 5 6 7];%size = L

A2=diag((exp(1i*2*pi*a/L)));
%A = [A1 A2];

y1= ctranspose(y);
calc= zeros(1,L);
 %b=zeros(8,8);
% X = [b b b b b b b b]; 
%GLOBAL unitary_array;
for i=1:8
   
    X=[A1*unitary_array(:,i) A2*unitary_array(:,i)];
    calc(i) = y1*X*ctranspose(X)*y;
end

[~,index]=max(calc);
%r=rem(index,4);
%if rem(index,4)~=0
 %   out=[a((index-r)/4);a(r)];
%else
 %   out=[a((index/4));a(4)];
    
%end
out = unitary_array(:,index);

end


    
    