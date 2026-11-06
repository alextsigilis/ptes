%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 Function description:
%                 
%--------------------------------------------------------------------------
%This function computes ar and ma coefficients of an arma-fitted model to
%patient data. First, it selects the channel and the entire set of epochs
%records for a particular sleep stage. Then it implies partitioning in
%epochs data vectors and it finds the average segment from all of
%them.Following this, it estimates the order of the ar and ma parts. In the end, it
%computes the ar and ma coefficients.
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
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%
%Return Products:
%
%==========================================================================
%ar:the computed coefficients from the ar part
%ma:the computed coefficients from the ma part
%p:the order of the ar part
%q:the order of the ma part







function [ar , ma ,p , q] = arma_estimate(Z,channel, sleep_stage,dt)


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
 
 p = arorder(A1,3, 10,10, 1);
 
 q = maorder(A1,0,9,0.05, 0);
 
 [ar, ma] = armarts(A1,p,q, 3, p+q ,L/10,0,'biased');
 
 