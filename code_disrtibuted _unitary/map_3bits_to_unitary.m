function out1= map_3bits_to_unitary(binary_symbol_vector,unitary_array)
% 1+j1 = s1
% 1-j1 = s2
% -1+j1 = s3
% -1-j1 = s4

%u=[7 60 79 187 125 198 154];
%theta=diag(exp(1i*2*pi/L*u))/sqrt(8);
% unitary_array=[s1 s2 s3 s4 s5 s6 s7 s8];
%unitary_array =zeros(L,L);
%for i= 1:L
 %  unitary=theta^(i-1);
  % for j=1:L;
   %    unitary_array(j,i) = unitary(j,j);
   %end
   
%end
binary_2_dec = bi2de(binary_symbol_vector); 
    

out1=unitary_array(:,binary_2_dec+1);
%out2=unitary_array;
end

