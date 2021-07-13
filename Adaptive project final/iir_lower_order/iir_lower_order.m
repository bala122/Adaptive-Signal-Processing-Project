time = 20;
N=8000*time; % 20s duration @ 8kHz sampling rate
power =1;
%training input data
noise = wgn(N,1,power);
Fs = 8000;

%Lower order IIR Filter target.
%Bandpass using butterworth iir method.
%Filter order 10
%Fc1 = 0.3 rad/s (normalized)
%Fc2 = 0.7 rad/s (normalized)


%NLMS Fir adaptive filter length 99
order =  98;
order2 =  100;
%making the filter
[b,a] = butter(5,[0.3 0.7],'bandpass');

wi_f = zeros(order+1,1);
%for changed mew
wi_fnew = zeros(order+1,1); 
wi_f2 = zeros(order2+1,1);

%step size
mew = 0.01;

'coeffs orig.'
'num'
b
'den'
a

% 1 experiment

% assuming a white noise input
noise = wgn(N,1,power);

x1 = transpose(noise); 


%d is the desired output
d1 = filter(b,a,x1);

zer1 = zeros(1,order);
zer2 = zeros(1,order2);

%zero padding both inputs
x2 = [zer2, x1];
x1 = [zer1, x1];

%esimtating
[wi A] = nlms(x1, d1, order,mew);
wi_f= wi+wi_f;

[wi2 A2] = nlms(x2, d1, order2,mew);
wi_f2= wi2+wi_f2;

%changed mew
mew2 = 0.001;
%1 experiment
exps=1;
for L=1:exps
    %Using 10 times the training data size.
noise2 = wgn(10*N,1,1*power);
x_new = transpose(noise2);
d_new = filter(b,a,x_new);
x_new = [zer1, x_new];
[wi_new A_new ]= nlms(x_new,d_new,order,mew2);
wi_fnew = wi_new+wi_fnew;
end



%estimate
b_est = wi_f;
b_est2 = wi_f2;
b_est_new = wi_fnew/exps;

'est coeffs fir'
wi_f
%plotting
h = fvtool(b,a,b_est,1);
title('Magnitude Resp. for Fir filter length 99')
h1 = fvtool(b,a,b_est,1);
legend (h,'original', 'estimate')
legend (h1,'original', 'estimate')
h1.Analysis='phase'
title('Phase Resp. for Fir filter length 99')

%phase





%The group delay seems to be the same, hence the phase is just off by some
%constant (approximately 2*pi) which means the phase is having a very low
%error actually.
%The graphs are different because of phase errors at dc and at normalized
%freq of 1 as the phase is more sensitive at these frequencies due to error
%in coeffs.

%Now we can see how to get the phase closer to the desired phase response
%by increasing the order 
 
h2 = fvtool(b,a,b_est2,1);
title('Magnitude Resp. for Fir filter length 101')

h3 = fvtool(b,a,b_est2,1);

legend (h2,'original', 'estimate')
legend (h3,'original', 'estimate')
h3.Analysis='phase'
title('Phase Resp. for Fir filter length 101')

%for order = 100 phase fits almost perfectly for most experiments. Hence,
%changing order is a good way to improvise the performance at a low cost.
%Here, we just tweaked the order by 1.

%Changed mew
%plotting
h6 = fvtool(b,a,b_est_new,1);
title('Magnitude Resp. for Fir filter length 99(changed mew)')
h7 = fvtool(b,a,b_est_new,1);
legend (h6,'original', 'estimate')
legend (h7,'original', 'estimate')
h7.Analysis='phase'
title('Phase Resp. for Fir filter length 99(changed mew)')
%we can see for this, even with 10 times the training data and decreasing step
%size, assuming it has reached steady state the mag. and phase is diverging
%in the ends of the spectrum. Hence, tweaking step size isn't a good way to
%improvise the response.



%An alternative way of predicting lower order iir response is
%to model the adaptive filter as iir itself. This would require length of
%numerator and denominator weight vectors( b and a).
%we compare this and the fir case.

noise = wgn(N,1,power);
x3 = transpose(noise); 

%d is the desired output
%d is same
d3 = filter(b,a,x3);

%order to be sent: order = length(b) + length(a)-1
filt_len3 = length(b) + length(a) -1;

%zero padding for initial samples
x3 = [ zeros(1,length(b)) x3];
d3 = [ zeros(1,length(a)) d3];

%Here the input regressors would be different, so we define the nlms
%algorithm separately here.
wi_3 = (zeros(filt_len3,1));  %weight vector initially zero
%size(wi);
mew3=0.01;            %step size
eps3 = 0.0001;        % epsilon chosen as a small positive parameter
A3=[];

for i= 1 :length(d3)-length(a)
    
    
    %extracting ui of length b
    c = i+length(b)-1; %index for x
    ui = flip(x3(c-(length(b)-1):c)); 
    
    %extracting di_block of length 'a' not including latest sample 
    c2 = i+length(a)-2;
    di_block = flip(d3(c2 - (length(a)-2):c2));

    %latest sample extracted separately or d(n)
    d_latest = d3(c2+1); 
    
    %concatenating di_block ( length = length(den.) -1) and ui till len(num.)
    ui_net = [ui di_block];
    
    
    %a check
    %'check'
    %d_latest - ui_net* transpose([ (b) (-a(2:end))]) %error 
    
    
    % so the resulting weight vector would be
    % [b1 b2 .....bm -a2 -a3....-an]
    %where m is num. length , n is den. length
    %a1 is 1 and we include d(n) corresponding to a1 in the nlms recursion.
    
    ei3 = d_latest - ui_net*wi_3; %error 
    wi_3 = wi_3 + (mew3/(eps3 + ui_net*ui_net'))* ui_net'* ei3; %estimating weights
    %ei3
    A3= [A3 wi_3'*wi_3];
end
w_out = wi_3;
A_out = A3;

'final est coeffs est.'

'numerator'
numf= transpose(wi_3(1:length(b)))
'actual num'
b

'denominator'
denf = [ 1 -transpose(wi_3(length(b)+1:end))]

'actual den'
a

h4 = fvtool(b,a,numf,denf);
title('Magnitude Resp. for adaptive iir model')

h5 = fvtool(b,a,numf,denf);

legend (h4,'original', 'estimate')
legend (h5,'original', 'estimate')
h5.Analysis='phase'

figure
t1 = linspace(0,time,length(A3)) ;
plot(t1,A3)
grid on
title('Progress of weight vector in Adaptive IIR Model')
xlabel('Time')
ylabel('Squared norm of weight vector')
%Hence, we see that even with modelling the iir filter as exactly with
%separate 'a', 'b' coeffs. We dont get the performance of the fir case with
% order 98 or 100
%intuitively this is because we involve di_block in the input regressors which
%isn't in our hands, while earlier we could generate ui and use a sliding
%window to give an estimate. Hence, the instantaneous approximation may not
%be so accurate in this case. Although it may converge better in the RLS
%algorithm.






function [w_out,A_out] = nlms(x,d,order,mew)

%size(x)
%size(d)
wi = (zeros(order+1,1));  %weight vector initially zero
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

