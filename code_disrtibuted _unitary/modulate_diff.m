function tx = modulate_diff(bits,unitary_array)
L=8;
N=size(bits,2);
tx = zeros(2,N/log2(L));
dummy = [1;0];
for i=1:3:N-2
    if i==1
        binary_2_dec = bi2de(bits(i:i+1)); 
        tx(:,i) = diag(unitary_array(:,binary_2_dec+1))*dummy;
    else
        binary_2_dec = bi2de(bits(i:i+1));
        tx(:,(i-1)/3 +1)=diag(unitary_array(:,binary_2_dec+1))*tx(:,(i-1)/3);
    end
end
end
