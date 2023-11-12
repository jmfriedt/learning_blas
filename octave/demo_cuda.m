N=8192
if (exist("state")!=0) % 2 for file
  load("state")
  rand ("state", state)
  printf("state loaded")
else
  state=rand ("state");
  save state state
  printf("state saved")
end

% a=rand(N,N)-0.5+j*(rand(N,N)-0.5);
% b=rand(N,N)-0.5+j*(rand(N,N)-0.5);
a=single(rand(N,N)-0.5);
b=single(rand(N,N)-0.5);
tic
start=clock();
% c=a*b;
a*b;
toc
% c(1,1)

% $ octave ./demo_cuda.m 
% N = 8192
% Elapsed time is 5.22813 seconds.
% ans =  3.8054 - 1.3838i

