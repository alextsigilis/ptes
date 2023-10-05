function [Bcep]=bicepstrum_stats(Z,patient)

p = size(Z);


channel = {'EEG_F4_M1' , 'EEG_C4_M1' , 'EEG_O2_M1' ,'EEG_C3_M2' , 'EMG_chin' ,'EOG_E1_M2' , 'EOG_E2_M2' , 'ECG' }; 

names = [ "var" "skw" "krt" "entropia" "squared_entropy" "Patient"  "channel"  "Annotations"];
types = [  "double" "double" "double" "double" "double"  "double" "string" "string"];
sz = [p(1) numel(types)];
Bcep = table('Size',sz,'VariableTypes',types,'VariableNames',names);

%Bcep=[];


  
       
        for i=1:p(1)
    
    
    
          B  = Z{i,1}{1,1};
          [~,~,bic,~]=bicepsf(B,20,length(B), 0,'unbiased', 256, 1);
          Bcep(i,:) =  { mean(std(bic)') mean(skewness(bic)') mean(kurtosis(bic)')  mean(abs(mean(bic .* log2(bic),'omitnan'))') mean(abs(mean(sqrt(bic) .* log2(sqrt(bic)),'omitnan'))') patient channel{1} Z.Annotations{i}};
    

        end
       
        
  

