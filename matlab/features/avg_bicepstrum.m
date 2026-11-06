%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 Function description:
%                 
%--------------------------------------------------------------------------
%This function computes the Bicepstrum for a specific sleep stage and a
%selected channel of a patient. To compute the Bicepstrum we implement
%partitioning in epochs data vectors (256 x 30 samples) according to the
%selected length time dt (partition an epoch into K segments with 256 x dt
%samples).Then, we find the average segment from all derived segments from partitioning
%and we calculate the bicepstral coefficients to this segment. 
%
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
%nfft: The samples for FFT computation
%--------------------------------------------------------------------------
%
%Return Products:
%
%==========================================================================
%Bceps:the bicepstral coefficients
%m:the 2nd moment measure of the calculated bicepstrum



function [Bceps,m]=avg_bicepstrum(Z,channel, sleep_stage,nfft,dt)

%axis determination
om = [-nfft/2:nfft/2-1] / nfft;


idx = Z.Annotations == sleep_stage;


m=zeros(sum(idx),1);

%frequency of sampling
fs=256;

% number of segments per epoch (epoch partitioning)
K = floor(30/dt) ; 

%samples per segment (epoch partitioning into K segments)
M = floor(dt * fs) ;


% sum1 is used for finding the average sleep stage values
sum1 = zeros(M,1);


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
p=0;



for i=1:1:height(Z)
    
    
    
 if string(Z.Annotations{i}) == sleep_stage 
     
    
     B  = Z{i,j}{1,1};
     L = length(B);
     p=p+1;
     [~,~,bic,~]=bicepsf(B,20,L, 0,'unbiased', nfft, 1);
     m(p)= mom2(bic,om);
     
     for k=1:1:K 
         
         if k*M<=L
             
             A = B(((k-1)*M+1):k*M);
             sum1 = sum1 + A;
             a=a+1;
             
         end
         
     end
  
 end
    
end



%A1 is the average signal from dt for all epochs for that specific sleep stage and
%channel selection
 A1 = sum1./a ;

 
 L=length(A1);
  
 
%  %Time vector 
%  t = (0:1:L-1)*(1/fs);
 
 
  


%Applying hanning window on the signal with the aim of smoothing

w = hanning(L, 'periodic') ;
 
A1 = A1 .* w ;

%Find the average bicepstrum coefficients 


 
[~,~,Bceps,~] = bicepsf(A1,20,L, 0,'unbiased', nfft, 1);
 
 
%Plot the average bicepstrum coefficients using FFT method

% figure(5)
% clf;
% hold on;
% xlabel( 'Lag m ');
% ylabel('Lag n');
% title('Bicepstrum')
% s = mesh( -nfft/2:nfft/2-1, -nfft/2:nfft/2-1 , abs(Bceps(:,:)));
% 
% s.FaceColor = 'interp' ;
% 
% colorbar;
% hold off;
% 
% figure(6)
% clf;
% om = [-nfft/2:nfft/2-1] / nfft;
% clf, contour(om,om,abs(Bceps),8), grid,
% title('bicepstrum')