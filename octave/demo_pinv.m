clear all
P=65536;
x=rand(P,1)-0.5;
signal=rand(P,1)-0.5;
xavant=x;
for p=[10 20]
  x(p:P)=x(p:P)+signal(1:end-p+1);
end
Nlag=44;
matrice=zeros(Nlag+1,length(x));
for N=0:Nlag
    matrice(N+1,N+1:P)=signal(1:end-N);
end
poidspinv=x'*pinv(matrice);
poids=x'*(matrice'*inv(matrice*matrice'));
size(matrice*matrice')
subplot(211)
plot(poids)
subplot(212)
plot(x')
hold on
plot(xavant')
plot(xavant'-(x'-poidspinv*matrice))
std(x')
std(xavant')
std(xavant'-(x'-poidspinv*matrice))
