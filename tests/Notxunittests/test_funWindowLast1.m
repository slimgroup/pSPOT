function test_funWindowLast1()
    p=matlabpool('size');
    for h=1:1:5
        for m=100:100:500
            for t=1:5
                % prepare input vector
                l=(3*h+randi(p*100,1))*p;
                odims=[m,l];
                x=distributed.randn(odims);
                x=x(:);
                n=prod(size(x));
                X=gather(x);
        
                % run function form with MPI
                tic;
                [ N xs ys ] = pSPOT.pWindow.funWindowLast1HaloShape(l,p,h);
                y=pSPOT.pWindow.funWindowLast1HaloMake(x,l,p,h);
                z=pSPOT.pWindow.funWindowLast1HaloAverage(y,l,p,h);
                v=pSPOT.pWindow.funWindowLast1HaloDrop(z,l,p,h);
                tm=toc;
                y=gather(y);
                z=gather(z);
                v=gather(v);
        
                % run operator form with MPI
                H=opdWindowLast1Halo(n,l,p,h);
                E=opdWindowLast1HaloAverage(n,l,p,h);
                tic;
                Y=H*x;
                Z=E*Y;
                V=H'*Z;
                TM=toc;
                Y=gather(Y);
                Z=gather(Z);
                V=gather(V);
        
                % run operator form with Kronecker
                A=oplWindow1Davg(l,p,h);
                D=opDirac(m);
                W=oppKron2Lo(A,D);
                tic;
                YY=W*x;
                ZZ=W*W'*YY;
                VV=W'*ZZ;
                KT=toc;
        
                fprintf('%3d %3d %3d %3d | %d |',h,m,l,N,t);
                fprintf('%1d %1d | %1d %1d |',isequal(z,y),isequal(X,v),isequal(Z,Y),isequal(X,V));
                fprintf('%6.3f %6.3f %5.2f %6.3f %5.2f\n',TM,tm,tm/TM,KT,KT/TM)
            end
        end
    end
end
