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

function test_piCon_empty
%% empty
n1 = randi(10);
n2 = randi(10);
n3 = randi(10);
A  = zeros(n1,n2,0);
A  = distributed(A);  
B  = piCon.empty(n1,n2,n3);
assertEqual(A,double( B ));
end % empty

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
% assertEqual( double( A .\ piCon(B) ), C); % EXTREME DANGER cuz
% distributed class does not know how to handle piCon
assertEqual( double( piCon(A) .\ piCon(B) ), C);
end % ldivide

function test_piCon_reshape
%% Testing piCon reshape
n1 = randi(5);
n2 = randi(5);
n3 = randi(5);
x  = piCon.randn(n1,n2,n3);
x  = reshape(x,n1,n2*n3);
x  = reform(x);
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

function test_piCon_permute
%% Testing piCon permute
n1 = randi([2 5]);
n2 = randi([2 5]);
n3 = randi([2 5]);
x  = piCon.randn(n1,n2,n3);
% x.codistInfo;
x  = permute(x,2,3,1);
% x.codistInfo;
end % permute