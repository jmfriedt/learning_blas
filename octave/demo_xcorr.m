% motif pseudo aleatoire `a trouver
code=rand(1,1024);code=code-mean(code);

% signal recu : code decal'e de 4 indices
signal=rand(1,16*1024);signal=signal-mean(signal);
for p=[4 11]
  signal(p:p+length(code)-1)=signal(p:p+length(code)-1)+code;
end

% matrice contenant les copies du code retard'e de 1 `a Nlag
for Nlag=[20 2000]
  matrice=zeros(2*Nlag+1,length(signal)); % longueur signal >> longueur code
  for N=-Nlag:Nlag
      matrice(N+Nlag+1,N+Nlag+1:Nlag+N+length(code))=code;
  end
  
  % correlation comme produit matriciel
  tic
  s1=signal*matrice.';
  toc
  figure
  plot(abs(s1))
  
  % correlation dans le domaine de Fourier
  hold on
  pkg load signal
  tic
  s2=xcorr(signal,code,2*Nlag)(2*Nlag+1:end);
  toc
  plot(abs(s2)+100)
end

% Elapsed time is 0.000374079 seconds. : matrix product, Nlag=20
% Elapsed time is 0.00329089 seconds.  : Fourier product, Nlag=20
% Elapsed time is 0.0215409 seconds.   : matrix product, Nlag=2000
% Elapsed time is 0.00253201 seconds.  : Fourier product, Nlag=2000
