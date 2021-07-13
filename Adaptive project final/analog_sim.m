
simOut = sim('analog_adapt2','ReturnWorkspaceOutputs','on')
'predicted weights of the filter'
wts_est = simOut.weights

%Comparing bode responses of original and est.
%original 
%order of analog lowpass butterworth= 8
order = 8;
fc = 2000;
fs = 8000;
[zb,pb,kb]= butter(order,fc*2*pi,'s');
[bb,ab] = zp2tf(zb,pb,kb);

%original
[hb,wb] = freqs(bb,ab,4096);
%estimate
[hb2,wb2] = freqz(wts_est,1,'whole',4096);

%plotting
%estimate
semilogy(wb2*fs/(2*pi),abs(hb2))
hold on
%original
semilogy(wb/(2*pi),abs(hb))
grid
hold off
title('Magnitude response of analog butterworth filter and estimate')
ylabel('Magnitude')
xlabel('Frequency in Hz')
legend('estimate(fir length 99)','original (analog)')
xlim([0,4000])
