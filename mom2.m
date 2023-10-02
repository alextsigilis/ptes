function s = mom2(B, waxis)
   [N,~] = size(B);
   s = 0;
   for i = 1 : N
      for j = 1 : N
        f1 = waxis(i);
        f2 = waxis(j);
        s = s + f1*f1 * f2*f2 * abs(B(i,j));
      end
    end
end