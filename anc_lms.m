close all
clear all
load('channel_4.mat');
Mp=32;
Ms=32;
I = 25;
N=5000;
SNR = 30;
var_v = 10^(-SNR/10);
mu=0.01;

learn_lms = zeros(1,N);

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
    
    for i=1:N
        u = [input(i) u(1:Mp-1)];
        v = [mod_input(i) v(1:Ms-1)];
        y = u*w_lms;
        yw = [y yw(1:Ms-1)];
        
        filt_u = u*channels;
        f_u = [filt_u f_u(1:Mp-1)];
        g(i)= u*channelp - yw*channels;
        
        w_lms = w_lms + mu*f_u'*g(i);
    end
      learn_lms = learn_lms + abs(g).^2;
end

learn_lms = 10*log10(learn_lms/I);
plot(1:N,learn_lms,'r');
title('Mean Square Error (M=32)')
xlabel('Iterations')
ylabel('MSE (dB)')
grid
axis tight