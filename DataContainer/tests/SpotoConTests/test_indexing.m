function test_suite = test_indexing
%test_indexing  Unit tests for indexing
initTestSuite;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setup
   iCon.randn('state',0); rand('state',0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_indexing_slices
   A = iCon.randn(4,5);
   B = opMatrix(A);
   idx1 = [1,3];
   idx2 = [2,4,5];
   assertElementsAlmostEqual(...
      A(:,:),...
      double(B(:,:)));
   assertElementsAlmostEqual(...
      A(:),...
      B(:));
   assertElementsAlmostEqual(...
      A(idx1,idx2),...
      double(B(idx1,idx2)));
   assertElementsAlmostEqual(...
      A(idx1,:),...
      double(B(idx1,:)));
   assertElementsAlmostEqual(...
      A(:,idx2),...
      double(B(:,idx2)));
   assertElementsAlmostEqual(...
      A(end-2:end,1:end),...
      double(B(end-2:end,1:end)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_indexing_absolute_slices
   A = iCon.randn(4,5);
   B = opMatrix(A);
   idx1 = [1,2,3,5; 6,10,11,18]';
   idx2 = [1,2,3,9; 11,12,17,18]';
   assertElementsAlmostEqual(...
      A(idx1),...
      B(idx1));
   assertElementsAlmostEqual(...
      A(idx2),...
      B(idx2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_indexing_logical
   A = iCon.randn(4,5);
   B = opMatrix(A);
   idx1 = [1,2,3,5; 6,10,11,18]';
   idx2 = [1,2,3,9; 11,12,17,18]';
   mask1 = zeros(size(A)); mask1(idx1) = 1; mask1 = (mask1 > 0.5);
   mask2 = zeros(size(A)); mask2(idx2) = 1; mask2 = (mask2 > 0.5);
   mask3 = [true, false, true];
   assertElementsAlmostEqual(...
      A(mask1),...
      B(mask1));
   assertElementsAlmostEqual(...
      A(mask2),...
      B(mask2));
   assertElementsAlmostEqual(...
      A(mask3),...
      B(mask3));
end