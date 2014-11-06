function out= ddstm_demod1_1(y0,y1,unitary_array,L)
calc = zeros(1,L);
%A1=eye(2);
%A2=[-1 0;0 1];
%X1 = [A1*y1 A2*y1];
%X0  = [A1*y0 A2*y1];

for i=1:L
    
    calc(i) = norm( (y1 - unitary_array(:,2*i-1:2*i)*y0));
end
[~,index]=max(calc);

out = index-1;
end
