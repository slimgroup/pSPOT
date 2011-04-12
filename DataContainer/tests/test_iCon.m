function test_suite = test_iCon
initTestSuite;
end

function test_iCon_conj
%% conj
n1 = randi(10);
n2 = randi(10);
A  = randn(n1,n2) + 1i*randn(n1,n2);
assertEqual(conj(A),double(conj(iCon(A))));
end % conj

function test_iCon_ctranspose
%% ctranspose
n1 = randi(10);
n2 = randi(10);
A  = randn(n1,n2) + 1i*randn(n1,n2);
B  = iCon(A);
assertEqual(A',double(B'));
end % ctranspose

function test_iCon_horzcat
%% horzcat
n1 = randi(10);
n2 = randi(10);
A  = randn(n1,n2) + 1i*randn(n1,n2);
B  = iCon(A);
assertEqual( [A A], double([B B]) );
end % horzcat

function test_iCon_imag
%% horzcat
n1 = randi(10);
n2 = randi(10);
A  = randn(n1,n2) + 1i*randn(n1,n2);
B  = iCon(A);
assertEqual( imag(A), double( imag(B) ) );
end % horzcat

function test_iCon_ldivide
%% ldivide
n1 = randi(10);
n2 = randi(10);
A  = randn(n1,n2) + 1i*randn(n1,n2);
B  = randn(n1,n2) + 1i*randn(n1,n2);
C  = A .\ B;
assertEqual( double( iCon(A) .\ B ), C);
assertEqual( double( A .\ iCon(B) ), C);
assertEqual( double( iCon(A) .\ iCon(B) ), C);
end % ldivide

function test_iCon_minus
%% minus
n1 = 6;
n2 = 1;
A  = randn(n1,n2) + 1i*randn(n1,n2);
B  = randn(n1,n2) + 1i*randn(n1,n2);
C  = A - B;
assertEqual( double( iCon(A) - B ), C);
assertEqual( double( A - iCon(B) ), C);
assertEqual( double( iCon(A) - iCon(B) ), C);
end % ldivide

function test_iCon_plus
%% plus
x = iCon.randn(5,4);
A = opGaussian(5,4);
y = A + x;
end

function test_iCon_bsxfun
%% bsxfun
n1 = randi([2 10]);
n2 = randi([2 10]);
A  = randn(n1,n2) + 1i*randn(n1,n2);
B  = iCon(A);
y  = randn(n1,1);
assertEqual( bsxfun(@minus,A,y), double( bsxfun(@minus,B,y) ) );
end

function test_iCon_reshape
%% Testing reshape function
n1 = randi(10);
n2 = randi(10);
n3 = randi(10);
n4 = randi(10);
x = iCon.randn(n1,n2,n3,n4);
x = reshape(x,n1*n2,n3*n4);
y = vec(x);
z = reshape(y,n1,n2,n3,n4);

end

function test_iCon_permute
%% Testing permute and unpermute
n1 = randi([2,10]);
n2 = randi([2,10]);
n3 = randi([2,10]);
x = iCon.randn(n1,n2,n3);
x = permute(x,[3 2 1]);
x = invpermute(x);
end