
function y = quantize(x,B)  
 
% The function rounds x into a binary fixed point 
% representation with B bits, including the sign bit. 

% If overflow or underflow occurs, the largest possible 
% fixed-point representation is returned. 
 
% The largest numbers that can be represented are
% +- 2^(B-1)-1

% If x is less than these extreme values, the routine finds
% the number of bits that are needed to represent the integer
% part of x and uses the remaining bits for its fractional part.

y = 0;
frac = x - fix(x);
if abs(x) >= (2^(B-1)-1) 
   if x > 0
      y = 2^(B-1)-1;    % overflow
   else
      y = -(2^(B-1)-1); % underflow
   end;   
else   
   for i=0:B-1,
      if abs(x) < 2^i; % i bits are needed to represent the integer part of x
         M = B-i-1;  % M bits are used for its fractional part
         y = (2^(-M)*round(frac/2^(-M))) + fix(x); 
         break
      end
  end
end
