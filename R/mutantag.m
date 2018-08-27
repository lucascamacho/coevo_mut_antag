dim=200;
matriz=ones(dim);
matriz=matriz-eye(dim);
tempo=200;
z=10.*rand(dim,1);
theta=10.*rand(dim,1);
intforca=0.9;
alfa=0.2;
resultados=zeros(tempo,dim);
antprob=0.1;
vmatriz=zeros(dim);
phi=0.2;
barreira=3;
% % of antagonistic links
for i=1:dim
    for j=1:dim
        if matriz(i,j)==1
            lixo=rand(1,1);
            if lixo<=antprob
                matriz(i,j)=0;
                vmatriz(i,j)=1;
            end
        end
    end
end
Q=matriz+vmatriz;
F=Q;

for tempoi=1:tempo
    
    %record the output
    for i=1:dim
        resultados(tempoi,i)=z(i);
    end
    
    %define interaction strengths
    for i=1:dim
        for j=1:dim
            if Q(i,j)>0
                Q(i,j)=F(i,j).*exp(-alfa.*((z(j)-z(i))^2));
            end
        end
    end
    
    %standardize the matriz
    for i=1:dim
        Q(i,:)=Q(i,:)./sum(Q(i,:));
    end
    
    Q=intforca*Q;
    S=zeros(dim,1);
    %compute the selection differentials
    for i=1:dim
        S(i)=S(i)+phi.*(1-intforca).*(theta(i)-z(i));
        for j=1:dim
            if matriz(i,j)>0
                S(i)=S(i)+phi.*Q(i,j).*(z(j)-z(i));
            else
                if vmatriz(i,j)>0
                    if abs(z(j)-z(i))<barreira
                        if z(i)>z(j)
                            S(i)=S(i)+phi.*Q(i,j).*(z(j)+barreira-z(i));
                        else
                            S(i)=S(i)+phi.*Q(i,j).*(z(j)-barreira-z(i));
                        end
                    end
                end
            end
        end
    end
    
    
            z=z+S;
    
end
x=linspace(1,tempo,tempo);
plot(x, resultados)

B=zeros(dim,1);
for i=1:dim
    for j=1:dim
        if z(j)>z(i)
        B(i)=B(i)-Q(i,j).*vmatriz(i,j).*barreira;
        else
            B(i)=B(i)+Q(i,j).*vmatriz(i,j).*barreira;
        end
    end
end

zest=inv(eye(dim)-Q)*(B+(1-intforca)*theta);
zest2=inv(eye(dim)-Q)*((1-intforca)*theta);
R=[z zest zest2]
mean(abs(zest-zest2)./zest)
scatter(z,zest)