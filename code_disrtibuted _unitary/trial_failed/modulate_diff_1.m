function tx = modulate_diff_1(bits,unitary_array,L)

s=log2(L);
N=size(bits,2);
tx = zeros(2,N/log2(L));
dummy = [1;0];
for i=1:s:N-s+1
    if i==1
        binary_2_dec = bi2de(bits(i:i+s-1)); 
        tx(:,i) = unitary_array(:,2*binary_2_dec+1:2*binary_2_dec+2)*dummy;
    else
        binary_2_dec = bi2de(bits(i:i+s-1));
        tx(:,(i-1)/s +1)=unitary_array(:,2*binary_2_dec+1:2*binary_2_dec+2)*tx(:,(i-1)/s);
    end
end
end
