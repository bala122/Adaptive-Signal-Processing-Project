time = 20;
N=8000*time; % 20s duration @ 8kHz sampling rate
power =1;
%training input data
noise = wgn(N,1,power);
Fs = 8000;
Ndft = fft(noise);
Ndft = Ndft(1:N/2+1);
psdN = (1/(Fs*N)) * abs(Ndft).^2;
psdN(2:end-1) = 2*psdN(2:end-1);
freq = 0:Fs/N:Fs/2;
plot(freq,10*log10(psdN))
grid on
title('Periodogram Using FFT of input signal')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')

%Lower order IIR Filter target.
%Bandpass using butterworth iir method.
%Filter order 50
%Fc1 = 0.3 rad/s (normalized)
%Fc2 = 0.7 rad/s (normalized)

%making the filter
[b,a] = butter(25,[0.3 0.7],'bandpass');
%lower order bandapss butterworth iir
[b2,a2] = butter(7,[0.3 0.7],'bandpass');

% Bandpass Chebyshev Type I
%same freq specs.
%Filter order 10
[b3,a3] = cheby1(5,1,[0.3 0.7]);


'coeffs orig.'
'num'
b
'den'
a


%We use an adaptive filter of order 99
order =  99;
%wi_f;
%wi_f2;
%wi_f3;

%step size
mew = 0.01;

%white input data
x1 = transpose(noise); 

%d is the desired output
d1 = filter(b,a,x1);
d2 = filter(b2,a2,x1);
d3 = filter(b3,a3,x1);

x1 = [zeros(1,order), x1];

%esimtating
[wi_f  A] = nlms(x1, d1, order,mew);

%lower order target
[wi_f2  A2] = nlms(x1, d2, order,mew);

%chebyshev filter
[wi_f3  A3] = nlms(x1, d3, order,mew);


%estimate
b_est = wi_f;
b_est2 = wi_f2;
b_est3 = wi_f3;


%plotting
h = fvtool(b,a,b_est,1);
title('Mag. Resp. for butterworth target filter order 50 using adaptive fir filter')
legend (h,'original(iir)', 'estimate(fir order 99)')

%As we can see for a higher order iir, it becomes hard to implement via an fir
%filter. For the given frequency specifications, all we can do is to reduce
%the order of the target filter to get it to converge. After hit and trial,
%for an adaptive filter of order 99 it's seen it slowly starts diverging at
%target filter order 14. We saw earlier that it converged for order = 10
h2 = fvtool(b2,a2,b_est2,1);
title('Magnitude Resp. for butterworth target filter order 14')
legend (h2,'original(iir)', 'estimate(fir order 99)')

%Moreover, it depends on the kind of iir filter, for sharp filters like
%chebyshev-type 1 it starts diverging earlier itself( around order 10 whereas 
% for butterworth iir order 10 we saw previously that the mag resp. converged )
h3 = fvtool(b3,a3,b_est3,1);
title('Magnitude Resp. for chebyshev target filter order 10')
legend (h3,'original(iir)', 'estimate(fir order 99)')




function [w_out,A_out] = nlms(x,d,order,mew)

%size(x)
%size(d)
wi = (zeros(order+1,1));  %weight vector initially zero
%wi = (10^-10)*randn(order+1,1)
%size(wi);
%mew=0.01;            %step size
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
