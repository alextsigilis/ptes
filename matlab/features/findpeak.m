function [x,y] = findpeak(b,f)
n = size(b,1);
x = [];
y = [];
bw = imregionalmax(b);
for i = 1 : n
    for j = 1 : n
        if bw(i,j) == 1
            x = [x f(j)];
            y = [y f(i)];
        end
    end
end
end