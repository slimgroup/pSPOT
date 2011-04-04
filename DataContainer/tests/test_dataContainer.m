function test_suite = test_dataContainer
initTestSuite;
end

function test_dataContainer_dcInCore
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
% disp('First reshape: ')
x = reshape(x,n1*n2,n3*n4);
% disp('Transpose')
x = x';
% disp('Reshape again: ')
x = reshape(x,n3,n4,n1,n2);
end

function test_dataContainer_plus
%% plus
clc
x = dcInCore.randn(5,4);
A = opGaussian(5,4);
y = A + x;
end

function test_dataContainer_bsxfun
%% bsxfun
clear, clc
x = dcInCore.randn(5,4);
y = randn(5,1);
z = bsxfun(@minus,x,y);
a = bsxfun(@minus,y,x);
end

function test_dataContainer_dcInCore_reshape
%% Testing intelligent reshape function
n1 = randi(10);
n2 = randi(10);
n3 = randi(10);
n4 = randi(10);

x = dcInCore.randn(n1,n2,n3,n4);
x = reshape(x,n1*n2,n3*n4);
cell2mat(isize(x))
y = vec(x);
cell2mat(isize(y))
z = reshape(y,n1,n2,n3,n4);

end

function test_dataContainer_dcInCore_permute
%% Testing permute and unpermute
n1 = randi([2,10]);
n2 = randi([2,10]);
n3 = randi([2,10]);
x = dcInCore.randn(n1,n2,n3);
x = permute(x,3,2,1);
x = unpermute(x);
end