close all
clear all
load('channel_4.mat');
Mp=4;
Ms=4;
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
    s_rls = zeros(Ms,1);
    
    g = zeros(1,N);
    f = zeros(1,N);
    u  = zeros(1,Mp);
    v  = zeros(1,Ms);
    yw = zeros(1,Ms);
    f_u = zeros(1,Mp);
    
   %for rls
   P = (1/epsilon)*eye(Ms,Ms); 
    P1 = (1/epsilon)*eye(Mp,Mp); 
    for i=1:N
        u = [input(i) u(1:Mp-1)];
        v = [mod_input(i) v(1:Ms-1)];
        y = u*w_lms;
        yw = [y yw(1:Ms-1)];
        
        filt_u = u*channels;
        f_u = [filt_u f_u(1:Mp-1)];
        g(i)= u*channelp - yw*channels;

        gamma1 = 1/(1 +(1/lambda)*f_u*P1*f_u'); 
        gg1 = (1/lambda)*P1*f_u'*gamma1;
        w_lms = w_lms + gg1*g(i);
         P1 = (1/lambda)*P1 - (gg1*gg1'/gamma1);
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