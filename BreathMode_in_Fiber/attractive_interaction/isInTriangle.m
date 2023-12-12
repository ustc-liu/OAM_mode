function z = isInTriangle(a,x,y)

% A = [0 sqrt(3)*a/3];
% B = [-a/2 -sqrt(3)*a/6];
% C = [a/2 -sqrt(3)*a/6];

k_AB = sqrt(3);
k_BC = 0;
k_AC = -sqrt(3);

y1 = k_AB*x + sqrt(3)/3*a;
y2 = -sqrt(3)/6*a;
y3 = -k_AC*x + sqrt(3)/3*a;

if (y <= y1) && ((y >= y2) && (y <= y3))
    z = true;
else
    z = false;
end