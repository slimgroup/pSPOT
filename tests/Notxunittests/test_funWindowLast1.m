function test_funWindowLast1()
p=matlabpool('size');
for h=1:1:5
    for m=100:100:500
    for t=1:5
        l=(3*h+randi(p*100,1))*p;
        odims=[m,l];
        x=distributed.randn(odims);
        x=x(:);
        n=prod(size(x));
        X=gather(x);

        tic;
        [y N PART]=pSPOT.pWindow.funWindowLast1HaloMake(x,l,p,h);
        [z N PART]=pSPOT.pWindow.funWindowLast1HaloAverage(y,l,p,h);
        [v N part]=pSPOT.pWindow.funWindowLast1HaloDrop(z,l,p,h);
        tm=toc;
        y=gather(y);
        z=gather(z);
        v=gather(v);

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

        A=oplWindow1Davg(l,p,h);
        D=opDirac(m);
        W=oppKron2Lo(A,D);
        tic;
        YY=W*x;
        ZZ=W*W'*YY;
        VV=W'*ZZ;
        KT=toc;

        fprintf('%3d %3d %3d | %d |',h,m,l,t);
        fprintf('%1d %1d | %1d %1d |',isequal(z,y),isequal(x,v),isequal(Z,Y),isequal(x,V));
        fprintf('%6.3f %6.3f %5.2f %6.3f %5.2f\n',TM,tm,tm/TM,KT,KT/TM)
    end
    end
end
end
