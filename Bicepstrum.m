
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 Function description:
%                 
%--------------------------------------------------------------------------
%This function computes Bicepstrum for a data vector. For this
%implementation , it requires the max lag for the cumulants calculation.
%The samples per segment for the data vector A in order to apply
%segmentation. Segmentation is used for computing the cumulants for each
%distinct segment and then it calculates the bicepstral coefficients using
%the FFT method.
%
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%
%Function Arguments:
%
%==========================================================================
%Z : the data table of a patient which contains the polysomnographic
%recordings
%nLag: number of lags to compute
%nsamp: samples per segment 
%nfft: The samples for FFT computation

%--------------------------------------------------------------------------
%
%Return Products:
%
%==========================================================================
%Bceps : the Bicepstrum index 



function [Bceps] = Bicepstrum(A , nLag , nsamp , nfft)



    nLag = min(nLag, nsamp-1);

    C1 = zeros(2*nLag+1, 2*nLag+1);
    Bceps = zeros(2*nLag+1, 2*nLag+1);

    
    
    if (nfft  < 2*nLag+1)  
        
        nfft = 2^nextpow2(nsamp); 
        
    end

   
   [~,~,C] = bispeci(A,nLag,nsamp,0,'b', nfft, 1);
   
  
   
   for i = -nLag:1:nLag
       
       for j = -nLag:1:nLag
       
         
                        
             C1(i+nLag+1,j+nLag+1) = i*C(i+nLag+1,j+nLag+1)  ;
       
           
       
       end
       
   end
   
%    C = C/max(abs(C),[],"all") ;
%    C1 = C1/max(abs(C1),[],"all") ;
   
   
   
   B = fftshift(ifft2(fft2(C1)./fft2(C)));
   
   %[rows , ~] = size(B);
   
   for i = -nLag:1:nLag
       
       for j = -nLag:1:nLag
           
           
          if i ~= 0
           
            Bceps(i+nLag+1,j+nLag+1) = (1/i) * B(i+nLag+1,j+nLag+1) ;
       
          else
           
            Bceps(i+nLag+1,j+nLag+1) =  B(i+nLag+1,j+nLag+1) ;
       
          end
       
       end
       
   end
   
   
    figure(1)
    clf;
    hold on;
    fontsize='\fontsize{10}';
    seconds='(\its \rm)';
    title([fontsize 'Bicepstrum (\itV \rm^{3} )']);
    xlabel([fontsize 'Lag m ' seconds]);
    ylabel([fontsize 'Lag n'  seconds]);
    imagesc( -nLag:1:nLag , -nLag:1:nLag , abs(B(:,:)));
    hold off;

    figure(2)
    clf;
    hold on;
    xlabel([fontsize 'Lag m ' seconds]);
    ylabel([fontsize 'Lag n'  seconds]);
    s = mesh( -nLag:1:nLag , -nLag:1:nLag , abs(B(:,:)));
    s.FaceColor = 'flat'
    hold off;
   
