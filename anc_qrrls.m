close all
clear all
load('channel_32.mat');
Mp=32;
M=32;
Ms=32;
I = 25;
N=5000;

SNR = 30;
var_v = 10^(-SNR/10);
mu=0.01;
lambda = 0.995;
epsilon = 10;
learn_rls = zeros(1,N);

for l=1:I
    l
    input = rand(1,N);
    mod_input = wgn(1,N,0);
    
    w_lms = zeros(Mp,1);
    
    g = zeros(1,N);
    f = zeros(1,N);
    u  = zeros(1,Mp);
    v  = zeros(1,Ms);
    yw = zeros(1,Ms);
    f_u = zeros(1,Mp);
    
   %for rls
   Phi = epsilon^(1/2)*eye(M,M); 
   q = zeros(M,1);
    for i=1:N
        u = [input(i) u(1:Mp-1)];
        v = [mod_input(i) v(1:Ms-1)];
        y = u*w_lms;
        yw = [y yw(1:Ms-1)];
        
        filt_u = u*channels;
        f_u = [filt_u f_u(1:Mp-1)];
        g(i)= u*channelp - yw*channels;
        f(i) = g(i) - v*channels;
        d1 = g(i)+f_u*w_lms;

A_qrrls = [lambda^(1/2)*Phi  f_u'; lambda^(1/2)*q' d1'; zeros(1,M) 1];
    [c_qrrls,s_qrrls] = givens(A_qrrls(1,1),A_qrrls(1,M+1));
    
    B_qrrls = A_qrrls*[c_qrrls zeros(1,M-1) s_qrrls; zeros(M-1,1) eye(M-1,M-1) zeros(M-1,1);conj(s_qrrls) zeros(1,M-1) -c_qrrls]; 

    for kk=2:M
        [c_qrrls,s_qrrls] = givens(B_qrrls(kk,kk),B_qrrls(kk,M+1));  
        B_qrrls = B_qrrls*[eye(kk-1,kk-1) zeros(kk-1,M-kk+2); zeros(1,kk-1) c_qrrls zeros(1,M-kk) s_qrrls; zeros(M-kk,kk) eye(M-kk,M-kk) zeros(M-kk,1);zeros(1,kk-1) conj(s_qrrls) zeros(1,M-kk) -c_qrrls];
    end
    
     Phi = B_qrrls(1:M,1:M);
     q = (B_qrrls(M+1,1:M))';
     w_lms = (inv(Phi))'*q;
    end
      learn_rls = learn_rls + abs(g).^2;
end

learn_rls = 10*log10(learn_rls/I);
plot(1:N,learn_rls,'r');
title('Mean Square Error (M=32)')
xlabel('Iterations')
ylabel('MSE (dB)')
grid
axis tight