%% Test for dcInCore
clc

n1 = 3;
n2 = 3;
n3 = 3;
n4 = 3;

x = zeros(n1,n2,n3,n4);
for l = 1:n4
    for k = 1:n3
        for j = 1:n2
            for i = 1:n1
                x(i,j,k,l) = i + j/10. + k/100. + l/1000.;
            end
        end
    end
end

disp(x)
disp('First reshape: ')
x = reshape(x,n1*n2,n3*n4)
disp('Transpose')
x = x'
disp('Reshape again: ')
x = reshape(x,n3,n4,n1,n2)