%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 Function description:
%                 
%--------------------------------------------------------------------------
%This function finds the average vector data of a type cell matrix, which
%cells store data vectors (epoch records).
%
%--------------------------------------------------------------------------

function [A]=find_avg(C)

sum = zeros(length(C{1,1}),1);



for i=1:1:length(C)
    
  B  = C{i,1};
  
  if isnan(B)==0
      
  sum = sum + B;
  
  end
    
end


A = sum./length(C) ;


