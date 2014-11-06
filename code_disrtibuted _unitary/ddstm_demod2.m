function out= ddstm_demod2(y0,y1,unitary_array,L)
calc = zeros(1,L);


for i=1:L
    
    calc(i) =  real(trace(y0*ctranspose(y1)*diag(unitary_array(:,i))));
end
[~,index]=max(calc);

out = index-1;
end