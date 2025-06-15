% main.m
% Author1:   Susana Lopez Moreno <susanalopezmoreno@pusan.ac.kr>
% Author2:   Jeong-Hoon Ju       <jjh793012@naver.com>
% Date:     2025-06-15
% Affiliation: Department of Mathematics, Pusan National University
% 
% Description:
%   Experiment checking the inequality of the determinantal identity
%   of the tensor geometric mean

syms x y
v = [ x, y ];  

%%%%%%%%%%%%%%%%%%%%%%%
% Resultant of tensor A
%%%%%%%%%%%%%%%%%%%%%%%

% Construction of A
A = zeros(2,2,2);
A(:,:,1) = [ 6  5
             10  11 ];
A(:,:,2) = [ 7  15
             10  22 ];

% Computation of the resultant of A
f_A=expand(v*squeeze(A(:,1,:))*v.');
g_A=expand(v*squeeze(A(:,2,:))*v.');

[C1, M1]= coeffs(f_A,[x,y]);
[C2, M2]= coeffs(g_A,[x,y]);

% Find the position of x^2, x*y, and y^2 in M1 and M2
[~, ix2] = ismember(x^2, M1);    
[~, ixy] = ismember(x*y, M1);    
[~, iy2] = ismember(y^2, M1);    

a1 = C1(ix2);  
a2 = C1(ixy);  
a3 = C1(iy2);

[~, ix2] = ismember(x^2, M2);    
[~, ixy] = ismember(x*y, M2);    
[~, iy2] = ismember(y^2, M2);    

b1 = C2(ix2);  
b2 = C2(ixy);  
b3 = C2(iy2);

% Calculation of the resultant
res_A = resultant_coeffs(a1,a2,a3,b1,b2,b3);

%%%%%%%%%%%%%%%%%%%%%%%
% Resultant of tensor B
%%%%%%%%%%%%%%%%%%%%%%%

% Construction of tensor B
B = zeros(2,2,2);
B(:,:,1) = [ 14  13
             26  31 ];
B(:,:,2) = [ 20  42
             28  60 ];

% Computation of the resultant of B
f_B=expand(v*squeeze(B(:,1,:))*v.');
g_B=expand(v*squeeze(B(:,2,:))*v.');

[C1_B, M1_B]= coeffs(f_B,[x,y]);
[C2_B, M2_B]= coeffs(g_B,[x,y]);

% Find the position of x^2, x*y, and y^2 in M1_B and M2_B
[~, ix2] = ismember(x^2, M1_B);    
[~, ixy] = ismember(x*y, M1_B);    
[~, iy2] = ismember(y^2, M1_B);    

c1 = C1_B(ix2);  
c2 = C1_B(ixy);  
c3 = C1_B(iy2);

[~, ix2] = ismember(x^2, M2_B);    
[~, ixy] = ismember(x*y, M2_B);    
[~, iy2] = ismember(y^2, M2_B);    

d1 = C2_B(ix2);  
d2 = C2_B(ixy);  
d3 = C2_B(iy2);

% Calculation of the resultant of B
res_B = resultant_coeffs(c1,c2,c3,d1,d2,d3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resultant of the tensor geometric mean G=A#_{*_M}B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computation of the tensor geometric mean
G = tensor_gm(A,B);

f_G=expand(v*squeeze(G(:,1,:))*v.');
g_G=expand(v*squeeze(G(:,2,:))*v.');

[C1_G M1_G]= coeffs(f_G,[x,y]);
[C2_G M2_G]= coeffs(g_G,[x,y]);

% find the position of x^2, x*y, and y^2 in M1
[~, ix2] = ismember(x^2, M1_G);    
[~, ixy] = ismember(x*y, M1_G);    
[~, iy2] = ismember(y^2, M1_G);    

e1 = C1_G(ix2);  
e2 = C1_G(ixy);  
e3 = C1_G(iy2);

[~, ix2] = ismember(x^2, M2_G);    
[~, ixy] = ismember(x*y, M2_G);    
[~, iy2] = ismember(y^2, M2_G);    

f1 = C2_G(ix2);  
f2 = C2_G(ixy);  
f3 = C2_G(iy2);

res_G = resultant_coeffs(e1,e2,e3,f1,f2,f3);


% Display the results
disp(sqrt(res_A*res_B));
disp(res_G);

% Compare results
if root~=res_G
    disp('not equal');
end