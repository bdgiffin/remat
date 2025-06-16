clear all
clc

N = 10;
A = zeros(N);
for i=1:10
    A(i,i) = 2;
    if (i > 1)
        A(i,i-1) = -1;
    end
    if (i < N)
        A(i,i+1) = -1;
    end
end

[L,D] = ldl(A)