
function y = quantize_v(x,B)  
 
% The same as quantize.m, except that x is a vector or a matrix.

[M N] = size(x);

for i=1:M
 for j=1:N
  y(i,j) = quantize(x(i,j),B);
 end
end
