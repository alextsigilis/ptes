


% names = ["Bicepstrum" "var" "skw" "krt" "entropia" "squared_entropy" "Patient"  "channel"  "Annotations"];
% types = [ "double" "double" "double" "double" "double" "double"  "double" "string" "string"];
% sz = [p(1) numel(types)];
% Bicepstrum_stat = table('Size',sz,'VariableTypes',types,'VariableNames',names);
names = [ "var" "skw" "krt" "entropia" "squared_entropy" "Patient"  "channel"  "Annotations"];
types = [  "double" "double" "double" "double" "double"  "double" "string" "string"];
sz = [140000 numel(types)];
Bicepstrum_stat = table('Size',sz,'VariableTypes',types,'VariableNames',names);
%Bicepstrum_stat=[];

a=0;
for i=[1:13,15:17,19:35,37:63,65:76,78:97,99:119]  %i=[50:63,65:76,78:97,99:119]
    
     S=read_data_of_patient(i);
     pp=size(S);
    
     Bcep=bicepstrum_stats(S,i);
     
     Bicepstrum_stat((a+1):(a+pp),:)= Bcep;
      a=a+pp;
end

% %first patient
% p1 = 1 ;
% 
% %last patient
% p2=119;
% %table for storing statistics
% names = ["Bicepstrum" "var" "skw" "krt" "entropia" "squared_entropy" "moment2nd" "Patient"  "channel"  "Annotations"];
% types = [ "double" "double" "double" "double" "double" "double" "double" "double" "string" "string"];
% sz = [1200*120 numel(types)];
% Bicepstrum_statistics = table('Size',sz,'VariableTypes',types,'VariableNames',names);
% 
% 
% 
% for k=1:8 
%     
%     
%     
%   for j = 1:5 
%     
%       Bcepsmat = zeros(nfft) ;
%       
%       a=0;
%       
% 
%          for i=[1:13,15:17,19:31,37:63,65:76,78:97,99:111,113:119] 
% 
% 
% 
% 
% 
%              S=read_data_of_patient(i);
% 
%              %A = S{3,1}{1,1} ;
%              
%                %calculating bicepstral coefficients
%             [Bceps,~]=avg_bicepstrum(S,channel{k}, sleepstages{j},nfft,30);
%             
%             Bicepstrum_stat(5*(p2-p1+1)*(k-1)+(i-p1)*5+j,:) = {Bceps std(Bceps)' skewness(Bceps)' kurtosis(Bceps)'  abs(mean(Bceps .* log2(Bceps),'omitnan'))' abs(mean(sqrt(Bceps) .* log2(sqrt(Bceps)),'omitnan'))' i channel{k} sleepstages{j}};
%              
%             
%          end
%    end
%   
% end
