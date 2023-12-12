function fy = myfunc(A,B,u)
% A is the eigen-value matrix
% B is the disorder-matrix
fy = -1i*(A+B)*u;
end