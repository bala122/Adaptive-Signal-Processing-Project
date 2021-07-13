time = 5;
N=8000*time; %  @ 8kHz sampling rate
power =1;
noise = wgn(N,1,power);
Fs = 8000;
periodogram(noise)
grid on
title('Periodogram of input signal')
xlabel('Frequency(normalized *pi)')
ylabel('Power/Frequency (dB/Hz)')

%FIR target with not a larger length
%Multiband filter designed by equiripple FIR.

order = 98;

%freq vector
f = [0, 0.28, 0.3, 0.48, 0.5, 0.69, 0.7, 0.8 ,0.81, 1];

%mag vector of mag. values for the respective frequencies.
a = [0, 0, 1, 1, 0, 0, 1, 1, 0, 0];

%weight vector for the weight of corresponding ripples in the frequency
%bands. Assume its the same ripple everywhere here
w = [1,1,1,1,1];

b =firpm(order,f,a,w);
'coeffs'
transpose(b)


%fir filter


%feeding input signal to Target filter and getting output d(n)

x1 = transpose(noise); % assuming a white noise input
%d is the desired output
d1 = filter(b,1,x1);

zer1 = zeros(1,order);
x1 = [zer1, x1];

%'A check'

%for k = 1:length(d) 
    
   % c= k+order; %index for x
   % d(k)-(flip(x(c-order:c)))*transpose(b)
%end

% found to be approximately same.

%NLMS algorithm

[wi A] = nlms(x1, d1, order);

%now we have the estimated weights



%estimate
b_est = wi
h = fvtool(b_est,1,b,1);
h1 = fvtool(b_est,1,b,1);
legend (h,'estimate', 'original');
legend (h1,'estimate', 'original');
h1.Analysis='phase'

%plotting wi progress
t1 = linspace(0,time,length(A)) ;
figure;
plot(t1,A)
grid
title('progress of weight matrix')
xlabel('time')
ylabel('inner product of weight matrix')

%similarly we can test it out for other kinds of design methods, like the
%hamming window based multiband fir filter, or a kaiser window based low pass 


%hamming window designed multiband filter of same order

order2 = 98;
low = 0.4;
bnd = [0.6, 0.9];
b2 = fir1(order2,[low,bnd]);

x2 = transpose(noise); % assuming a white noise input

d2 = filter(b2,1,x2);
%d is the desired output

zer2 = zeros(1,order2);
x2 = [zer2, x2];

[wi2 A2] = nlms(x2, d2, order2);

%estimate

b_est2 = wi2
h2 = fvtool(b_est2,1,b2,1);
h3 = fvtool(b_est2,1,b2,1);

legend (h2,'estimate', 'original');
legend (h3,'estimate', 'original');
h2.Analysis='Magnitude'

h3.Analysis='phase'
title('Phase via Hamming window')


%kaiser
%making a low pass filter of order 108 via kaiser window
fcuts = [1000 1500];
mags = [1 0];
devs = [0.05 0.01];
%3*Fs is done to get a filter order of 108, we anyway deal with normalized
%frequencies in the digital domain.
[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,3*Fs);
hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');
order3 = n;
b3 = hh;
%feeding input signal to Target filter and getting output d(n)

x3= transpose(noise); % assuming a white noise input

d3 = filter(b3,1,x3);
%d is the desired output

zer3 = zeros(1,order3);
x3 = [zer3, x3];

[wi3 A3] = nlms(x3, d3, order3);

%now we have the estimated weights

%estimate
b_est3 = wi3
h4 = fvtool(b_est3,1,b3,1);
h5 = fvtool(b_est3,1,b3,1);
legend (h4,'estimate', 'original');
legend (h5,'estimate', 'original');
h4.Analysis='Magnitude'
h5.Analysis='phase'
title('Phase via Kaiser window method')



%NLMS algorithm
function [w_out,A_out] = nlms(x,d,order)

size(x)
size(d)
wi = (zeros(order+1,1));  %weight vector initially zero
size(wi)
mew=0.01;            %step size
eps = 0.0001;        % epsilon chosen as a small positive parameter
%ei = 1;
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
