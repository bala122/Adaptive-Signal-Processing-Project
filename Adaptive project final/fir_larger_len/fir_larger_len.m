time = 20;
N_samp=8000*time; % 20s duration @ 8kHz sampling rate
power =1;
%training input data
noise = wgn(N_samp,1,power);
Fs = 8000;

%Hann window lowpass filter
Fc = 0.5;
N = 20;
lowpass_hann = designfilt('lowpassfir','FilterOrder',N,'CutoffFrequency',Fc*Fs/2, 'Window',hann(N+1),'SampleRate',Fs);
coeffs = lowpass_hann.Coefficients;
%parameters of comparing fir filter is same, order = 18;
%two different orders of adaptive and comparing filters
N_comp=18;
N_adap = 18;
lowpass_hann_comp = designfilt('lowpassfir','FilterOrder',N_comp,'CutoffFrequency',Fc*Fs/2, 'Window',hann(N_comp+1),'SampleRate',Fs);
x1 = transpose(noise); % assuming a white noise input
%d is the desired output
d1 = filter(coeffs,1,x1);
zer1 = zeros(1,N_adap);
x1 = [zer1, x1];
%NLMS algorithm
mew1= 0.01;
[w1 A] = nlms(x1, d1, N_adap,mew1);
Hd_adap = dfilt.df1(w1,1);
%Plotting
hfvt = fvtool(lowpass_hann,lowpass_hann_comp,Hd_adap,'Color','white')
legend('Target FIR(order 20)','Comparing FIR(order 18)','Adaptive FIR(order 18)')
title('Hanning window design Adaptive FIR')

%An observation noted is that with increasing order, the affect on
% stopband performance with slight changes in filter order is high.
%So this is experimented comparing for two cases
% original: 20,18 and new: 20 17
%original: 70 68 and new: 70 67

%order 20 target fir
lowpass_hann1 = lowpass_hann;
%order 18 adaptive fir
Hd_adap1 = Hd_adap; 

%for order 17 adaptive fir
N_adap2 = 17;
N_adap5 = 16;
x2 = transpose(noise); % assuming a white noise input
x5 = transpose(noise); % assuming a white noise input
%coeffs are target fir (order 20) coeffs
d2 = filter(coeffs,1,x2);
d5 = filter(coeffs,1,x5);
x2 = [zeros(1,N_adap2) , x2];
x5 = [zeros(1,N_adap5) , x5];
%NLMS algorithm
mew1= 0.01;
[w2 A2] = nlms(x2, d2, N_adap2,mew1);
Hd_adap2 = dfilt.df1(w2,1);
[w5 A5] = nlms(x5, d5, N_adap5,mew1);
Hd_adap5 = dfilt.df1(w5,1);
%plotting
hfvt2 = fvtool(lowpass_hann1,Hd_adap1,Hd_adap2,Hd_adap5,'Color','white')
title('Hanning window design adaptive fir(low order)')
legend('Target FIR(order 20)','Adaptive FIR(order 18)','Adaptive FIR(order 17)','Adaptive FIR(order 16)')

%order 70 target fir
N2 = 70;
lowpass_hann2 = designfilt('lowpassfir','FilterOrder',N2,'CutoffFrequency',Fc*Fs/2, 'Window',hann(N2+1),'SampleRate',Fs);
coeffs2 = lowpass_hann2.Coefficients;
%order 68 and 67 adaptive fir
N_adap3 = 68;
N_adap4 = 67;
x3 = transpose(noise); % assuming a white noise input
%coeffs are target fir (order 70) coeffs
d_70 = filter(coeffs2,1,x3);
x3_1 = [zeros(1,N_adap3) , x3];
x3_2 = [zeros(1,N_adap4) , x3];
%NLMS algorithm
mew1= 0.01;
[w3_1 A3_1] = nlms(x3_1, d_70, N_adap3,mew1);
Hd_adap3 = dfilt.df1(w3_1,1); %order 68

[w3_2 A3_2] = nlms(x3_2, d_70, N_adap4,mew1);
Hd_adap4 =dfilt.df1(w3_2,1); %order 67
%plotting
hfvt3 = fvtool(lowpass_hann2,Hd_adap3,Hd_adap4,'Color','white')
title('Hanning window design adaptive fir(high order)')
legend('Target FIR(order 70)','Adaptive FIR(order 68)','Adaptive FIR(order 67)')


function [w_out,A_out] = nlms(x,d,order,mew)

wi = (zeros(order+1,1));  %weight vector initially zero
eps = 0.0001;        % epsilon chosen as a small positive parameter
ei = 1;
A=[];
for i= 1 :length(d)
    di = d(i); %at time i
    c = i+order; %index for x
    ui = flip(x(c-order:c)); % extracting inputs of size = filter order +1
    
    ei = di - ui*wi; %error 
    wi = wi + (mew/(eps + ui*ui'))* ui'* ei; %estimating weights
    ei;
    A= [A wi'*wi];
end
w_out = wi;
A_out = A;
end