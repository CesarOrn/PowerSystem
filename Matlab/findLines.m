%coupling: is the coupling matrix of the network we are looking at.
%eigV: is the eigenvector of the laplacian matrix.

%data are the Theta_j - Theta_i coefficients
%line is a cell that contian the Theta_j- Theta_i.

function [line,data]= findLines(coupling,eigV)
    k=0;
    dataSize=size(eigV,1);
    % looks a upper triangule of the matrix. If A(i,j) is not equal to 0. We have
    % a connection to j,i.
    for i = 1:dataSize
        upper=dataSize-i;
        for j = dataSize-upper:dataSize
            if(coupling(i,j)~=0)
                k=k+1;
                linePower=eigV(j,:)-eigV(i,:);
                lineStr=strcat('Line:',num2str(j),',',num2str(i));
                line{k,:}=[lineStr];
                data(k,:)=linePower;

            end
        end
    end

end
