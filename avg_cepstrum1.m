%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Function description:                                   %
%--------------------------------------------------------------------------
%This function computes real cepstrum, complex
%cepstrum , complex cepstrum using fft(phase unwrapping  method) , minimum phase 
%impulse response and mixed-phase impulse response. Also for these calculations 
%implies signal vector partitioning according to the time length dt. For a
%specific channel and a selected sleep stage, it partitions epochs into K
%segments with 256 x dt samples. Then, it finds the average segment of the whole set of segments
%that are derived from partitioning all epochs records. In the end , it computes
%the products using the values of the average segment.
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%
%Function Arguments:
%
%==========================================================================
%Z : the data table of a patient which contains the polysomnographic
%recordings
%channel:the selected channel(EEG_F4_M1 , EEG_C4_M1 , EEG_O2_M1 , EEG_C3_M2 , EMG_chin , EOG_E1_M2 , EOG_E2_M2 , ECG)
%sleep stage: the selected sleep stage()
%dt: dt (in seconds) is an integer < 30 seconds 

%--------------------------------------------------------------------------
%
%Return Products:
%
%==========================================================================
%RC: the computed real cepstrum index
%CC1: the computed complex cepstrum index
%CC2: the computed complex cepstrum index( FFT method without using
%unwrapping phase)
%Hm: the estimated minimum phase impulse response
%H: the estimated mixed-phase impulse response


function [RC, CC1, CC2 ,Hm, H]=avg_cepstrum1(Z,channel, sleep_stage,dt)






%frequency of sampling
fs=256;

% number of segments per epoch (epoch partitioning)
K = floor(30/dt) ; 

%samples per segment (epoch partitioning into K segments)
M = floor(dt * fs) ;


% sum is used for finding the average sleep stage values
sum = zeros(M,1);


% channel matching with  number of the column in Z matrix
if channel == "EEG_F4_M1"
    
    j = 1;
    
elseif   channel == "EEG_C4_M1"
    
     j = 2;
    
elseif   channel == "EEG_O2_M1"
    
     j = 3;
    
elseif   channel == "EEG_C3_M2"
    
     j = 4;
    
elseif   channel == "EMG_chin"
    
     j = 5;
    
elseif   channel == "EOG_E1_M2"
    
     j = 6;
    
elseif   channel == "EOG_E2_M2"
    
     j = 7;
    
else
    
     j = 8;
    
end

a=0;




for i=1:1:height(Z)
    
    
    
 if string(Z.Annotations{i}) == sleep_stage 
     
    
     B  = Z{i,j}{1,1};
     L = length(B);
     
     
     for k=1:1:K 
         
         if k*M<=L
             
             A = B(((k-1)*M+1):k*M);
             %A = A(~isnan(A)) ;
             sum = sum + A;
             a=a+1;
             
         end
         
     end
  
 end
    
end



%A1 is the average signal from 30s epochs for that particular sleep stage and
%channel selection
 A1 = sum./a ;

 A1 = A1(~isnan(A1)) ;
 
 L=length(A1);
  
 
 %Time vector 
 t = (0:1:L-1)*(1/fs);
 
 
  
 %We subtract the mean value of the signal A1 from A1 (detrend) and then calculate
 %the Fourier Transform of A1 using fft. This is an effort to find
 %periodicities in the signal A1.
 
%  A1=detrend(A1,1);
 

 
 
%  Y=fft(A1);
%  P2 = abs(Y/L)';
%  P1 = P2(1:round(L/2+1));
%  P1(2:end-1) = 2*P1(2:end-1);
%  
%  f = fs*(0:(L/2))/L;
 
%  
%  figure(1)
%  clf;
%  hold on;
%  plot(f,P1) 
%  title("Single-Sided Amplitude Spectrum of A1(t)")
%  xlabel("f (Hz)")
%  ylabel("|P1(f)|")
%  hold off;
 
%  figure(2)
%  clf;
%  plot(1:L,A1)
%  
 
%  figure(5)
%  clf;
% findpeaks(P1(1:length(P1)/2),1:length(P1)/2)
% [p,peaks]=findpeaks(P1(1:length(P1)/2),1:length(P1)/2,'MinPeakDistance',2,'MinPeakHeight',max(P1)/4);
% m=mean(peaks);
 

%Applying hanning window on the signal with the aim of smoothing

w = hanning(L, 'periodic') ;
 
A1 = A1 .* w ;

%Find the average cepstrum coefficients for every case (real cepstrum , complex cepstrum , complex cepstrum existing from bicepstrum )

[RC , Hm] = rceps1(A1);
 
CC1 = cceps(A1) ; 
 
[H,CC2,~,~] = bicepsf(A1,round(sqrt(L)),L, 0,'unbiased', 512, 0);
 
 
%Plot the average cepstrum coefficients for every case (real cepstrum , complex cepstrum , complex cepstrum existing from bicepstrum )

% figure(4)
% clf;
% subplot(4,1,1)
% plot(t(1:round(L/2)) , RC(1:round(L/2)) , 'r'); title(sprintf('Average Complex Cepstrum of channel "%s" for the sleep stage %s ', channel , sleep_stage)); xlabel('Time(s)'); ylabel('Cepstrum');
% 
% 
% subplot(4,1,2)
% plot(t(1:round(L/2)) , CC1(1:round(L/2)) , 'r'); title(sprintf('Average Real Cepstrum of channel "%s" for the sleep stage %s ', channel , sleep_stage)); xlabel('Time(s)'); ylabel('Cepstrum');
%  
% 
% subplot(4,1,3)
% plot(-256/2:256/2-1, CC2),title(sprintf(' Complex Cebtrum using bicepstrum(FFT method) of channel "%s" for the sleep stage %s ', channel , sleep_stage)); xlabel('samples'); ylabel('Cepstrum'); grid on
%  
% subplot(4,1,4)
% plot(-256/2:256/2-1, H),title(sprintf(' Impulse response of channel "%s" for the sleep stage %s ', channel , sleep_stage)); xlabel('samples'); ylabel('Impulse response'); grid on