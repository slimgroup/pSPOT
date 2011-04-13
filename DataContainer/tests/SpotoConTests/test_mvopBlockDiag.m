function test_suite = test_mvopBlockDiag
%test_opBlockDiag  Unit tests for the opBlockDiag operator
initTestSuite;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seed = setup
   iCon.randn('state',0);
   rand('state',0);
   seed = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function test_mvopBlockDiag_prod(~)
   n = randi(100); m = randi(100);
   A = opMatrix(iCon.randn(m,m));
   B = opMatrix(iCon.randn(n,n));
   D = opBlockDiag(A,B);
   x = iCon.randn(m+n,2);
   y = [A*x(1:m,:); B*x(m+1:end,:)];
   assertElementsAlmostEqual( y, D*x )
   % assertElementsAlmostEqual is used instead of assertEqual because of a
   % inherent error (of the magnitude e-15) of the matlab built-in A'*x
   
   y2 = [A'*x(1:m,:); B'*x(m+1:end,:)];
   assertElementsAlmostEqual( y2, D'*x )
end

function test_mvopBlockDiag_overlap(~)
   m1 = randi([5,100]); n1 = randi([5,100]);
   m2 = randi([5,100]); n2 = randi([5,100]);
   A  = opMatrix(iCon.randn(m1,n1));
   B  = opMatrix(iCon.randn(m2,n2));
   
   % row overlap
   ov = 5;
   D  = opBlockDiag(A,B,ov);
   assertFalse(dottest(D));
   
   % column overlap
   ov = -5;
   D  = opBlockDiag(A,B,ov);
   assertFalse(dottest(D));

   % row-anti-diag overlap
   ov = m1+1;
   D  = opBlockDiag(A,B,ov);
   assertFalse(dottest(D));

   % col-anti-diag overlap
   ov = -(m1+1);
   D  = opBlockDiag(A,B,ov);
   assertFalse(dottest(D));
end