M=[12 20];
P=4;
H=3;

disp('start test');

for m=M(1):M(2)
    for p=2:P
        for h=0:H
    
            [A a]=pWindow.fd_window_1Dfor(m,p,h);
            [B b]=pWindow.fd_window_1Dbck(m,p,h);

            try
                n=m+(p-1)*2*h;
                x0=rand(m,1);
                y=A*x0;
                x1=B*y;
                check=sum(abs((x1-x0)'));
                fprintf('\tm=%d n=%d p=%d h=%d: %d\n',m,n,p,h,check);
            catch me
                fprintf('m=%d n=%d p=%d h=%d\n',m,n,p,h);
                whos a A b B
            end


        end
    end
end

disp('end test');
