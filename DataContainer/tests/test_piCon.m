function test_suite = test_piCon
initTestSuite;
end

function test_piCon_conj
%% conj
n1 = randi(10);
n2 = randi(10);
A  = randn(n1,n2) + 1i*randn(n1,n2);
A  = distributed(A);
assertEqual(conj(A),double(conj(piCon(A))));
end % conj

function test_piCon_ctranspose
%% ctranspose
n1 = randi(10);
n2 = randi(10);
A  = randn(n1,n2) + 1i*randn(n1,n2);
A  = distributed(A);
B  = piCon(A);
assertEqual(A',double(B'));
end % ctranspose

function test_piCon_horzcat
%% horzcat
n1 = randi(10);
n2 = randi(10);
A  = randn(n1,n2) + 1i*randn(n1,n2);
A  = distributed(A);
B  = piCon(A);
assertEqual( [A A], double([B B]) );
end % horzcat

function test_piCon_imag
%% horzcat
n1 = randi(10);
n2 = randi(10);
A  = randn(n1,n2) + 1i*randn(n1,n2);
A  = distributed(A);
B  = piCon(A);
assertEqual( imag(A), double( imag(B) ) );
end % horzcat

function dtest_piCon_ldivide
%% ldivide
n1 = randi(10);
n2 = randi(10);
A  = randn(n1,n2) + 1i*randn(n1,n2);
B  = randn(n1,n2) + 1i*randn(n1,n2);
A  = distributed(A);
B  = distributed(B);
C  = A .\ B;
assertEqual( double( piCon(A) .\ B ), C);
assertEqual( double( A .\ piCon(B) ), C); % EXTREME DANGER cuz
% distributed class does not know how to handle piCon. NOT ANYMORE CUZ I
% SPECIFIED DISTRIBUTED AS INFERIOR CLASS NYAHAHAHAHAHAHAHAHAHA~~~~~~
assertEqual( double( piCon(A) .\ piCon(B) ), C);
end % ldivide

function test_piCon_minus
%% minus
n1 = randi(10);
n2 = randi(10);
A  = distributed( randn(n1,n2) + 1i*randn(n1,n2) );
B  = distributed( randn(n1,n2) + 1i*randn(n1,n2) );
C  = A - B;
assertEqual( double( piCon(A) - B ), C);
assertEqual( double( A - piCon(B) ), C);
assertEqual( double( piCon(A) - piCon(B) ), C);
end % minus

function test_piCon_mldivide
%% mldivide
n1 = randi(10);
n2 = n1;
A  = distributed( randn(n1,n2) + 1i*randn(n1,n2) );
x  = distributed( randn(n2,2)  + 1i*randn(n2,2)  );
y  = A*x;
x  = gather(x);
assertElementsAlmostEqual( gather( double( piCon(A) \ y ) ), x);
assertElementsAlmostEqual( gather( double( A \ piCon(y) ) ), x);
assertElementsAlmostEqual( gather( double( piCon(A) \ piCon(y) ) ), x);
end % mldivide

function test_piCon_mrdivide
%% mrdivide
n1 = randi(10);
n2 = n1;
A  = distributed( randn(n1,n2) + 1i*randn(n1,n2) );
x  = distributed( randn(n2,2)  + 1i*randn(n2,2)  );
y  = A*x;
x  = gather(x);
assertElementsAlmostEqual( gather(double( piCon(y') / A' )'), x);
assertElementsAlmostEqual( gather(double( y' / piCon(A') )'), x);
assertElementsAlmostEqual( gather(double( piCon(y') / piCon(A') )'), x);
end % mrdivide

function test_piCon_mtimes
%% mtimes
n1 = randi(10);
n2 = randi(10);
A  = distributed( randn(n1,n2) + 1i*randn(n1,n2) );
B  = distributed( randn(n2,2)  + 1i*randn(n2,2)  );
C  = A*B;
C  = gather(C);
assertEqual( gather(double( piCon(A) * B )), C);
assertEqual( gather(double( A * piCon(B) )), C);
assertEqual( gather(double( piCon(A) * piCon(B) )), C);
end % mtimes

function test_piCon_permute
%% Testing permute and unpermute
n1 = randi([2,10]);
n2 = randi([2,10]);
n3 = randi([2,10]);
x  = piCon.randn(n1,n2,n3);
assertEqual(invpermute(permute(x,[3 2 1])),x);
end % permute

function test_piCon_plus
%% minus
n1 = randi(10);
n2 = randi(10);
A  = distributed( randn(n1,n2) + 1i*randn(n1,n2) );
B  = distributed( randn(n1,n2) + 1i*randn(n1,n2) );
C  = A + B;
assertEqual( double( piCon(A) + B ), C);
assertEqual( double( A + piCon(B) ), C);
assertEqual( double( piCon(A) + piCon(B) ), C);
end % plus

function test_piCon_power
%% plus
n1 = randi(10);
n2 = randi(10);
A  = distributed( randn(n1,n2) + 1i*randn(n1,n2) );
n  = randi(10);
C  = A .^ n;
assertEqual( double( piCon(A) .^ n ), C);
assertEqual( double( A .^ piCon(n) ), C);
assertEqual( double( piCon(A) .^ iCon(n) ), C);
end % power

function test_piCon_rdivide
%% rdivide
n1 = randi(10);
n2 = randi(10);
A  = distributed( randn(n1,n2) + 1i*randn(n1,n2) );
B  = distributed( randn(n1,n2) + 1i*randn(n1,n2) );
C  = A ./ B;
assertEqual( double( piCon(A) ./ B ), C);
assertEqual( double( A ./ piCon(B) ), C);
assertEqual( double( iCon(A) ./ piCon(B) ), C);
end % rdivide

function test_piCon_real
%% real
n1 = randi(10);
n2 = randi(10);
A  = distributed( randn(n1,n2) + 1i*randn(n1,n2) );
assertEqual( double(real(piCon(A))), real(A) );
end % real

function test_piCon_reshape
%% Testing piCon reshape
n1 = randi(5);
n2 = randi(5);
n3 = randi(5);
x  = piCon.randn(n1,n2,n3);
x  = reshape(x,n1,n2*n3);
x  = invec(x);
end % reshape

function test_piCon_redistribute
%% Testing piCon reshape
n1 = randi(5);
n2 = randi(5);
n3 = randi(5);
x  = piCon.randn(n1,n2,n3);
% x.codistInfo;
x  = redistribute(x,2);
% x.codistInfo;
end % redistribute