
z=[42 52 48 58; 4 5 4 3];
a= mean(z);

function y = mean(x)
    if ~ismatrix(x)
        error('Input must be a matrix')
    end
    [m,n] = size(x);
    y=[];
    for i = 1:1:m
        sum = 0;
        for j = 1:1:n
           sum = sum + x(i,j);
        end
        y=[y sum/n]
    end
end

function y = variance(x)
    if ~ismatrix(x)
        error('Input must be a matrix')
    end
    [m,n] = size(x);
    y=[];
    for i = 1:1:m
        sum = 0;
        for j = 1:1:n
           sum = sum + x(i,j);
        end
        y=[y; sum/n]
    end
end