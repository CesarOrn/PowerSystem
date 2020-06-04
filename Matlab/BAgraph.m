function A=BAgraph(dimF)
    dim=dimF;%number of nodes
    A=[0 1 1; 1 0 1; 1 1 0];%adjacency matrix
    A=sparse(A);%information stored as a sparse matrix

    kmedio=3;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=(length(A)+1):dim
        z=full(sum(A));
        zc=cumsum(z/sum(z));
        for jj=1:kmedio
            a1=find((zc-rand)>0);
            A(j,a1(1))=1;
            A(a1(1),j)=1;
        end
    end
    z=full(sum(A));
    edges=dim*kmedio;
    grmin=min(z);
    grmax=max(z);

    for i=grmin:grmax
        pf(i)=sum(z==i);
        if pf(i)>0
            grad(i)=i;
        end
    end


end
%set(figure)
%plot(log(grad(pf>0)),log(pf(pf>0)),'k*')

