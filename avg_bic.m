%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 Function description:
%                 
%--------------------------------------------------------------------------
%This function computes the Bicoherecy index for a specific sleep stage and a
%selected channel of a patient. To compute the Bicoherency we implement
%partitioning in epochs data vectors (256 x 30 samples) according to the
%selected length time dt (partition an epoch into K segments with 256 x dt
%samples).Then, we find the average segment from all derived segments from partitioning
%and we calculate the bicoherency index to this segment. 
% 
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
%bic:the bicoherency index
%waxis:the axis values for bicoherency graphical representation 




function [bic,waxis] = avg_bic(Z,channel, sleep_stage,dt,nfft)

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
    
    
    
    
 %Partitioning implementation along the ensemble epochs   
    
 if string(Z.Annotations{i}) == sleep_stage 
     
     %data vector
     B  = Z{i,j}{1,1};
     L = length(B);
     
     
     for k=1:1:K 
         
         if k*M<=L
             
             %partitioning of each epoch
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
 
 %Calculating the bicoherence 
 [bic,waxis] = bicoher(A1,  nfft, 0, L/10, 0);
 
 
%   hold off, clf
%   contour(waxis,waxis,bic,4), grid on 
%   title('Bicoherence estimated via the direct (FFT) method')
%   xlabel('f1'), ylabel('f2')
%   set(gcf,'Name','Hosa BICOHER')
