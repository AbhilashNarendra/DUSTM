function out= ddstm_demod2_1(y0,y1,unitary_array,L)
calc = zeros(1,L);


for i=1:L
    
    calc(i) =  real(trace(y0*ctranspose(y1)*unitary_array(:,2*i-1:2*i)));
end
[~,index]=max(calc);

out = index-1;
end