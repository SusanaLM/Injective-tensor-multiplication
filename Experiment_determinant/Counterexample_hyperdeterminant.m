% Counterexample_hyperdeterminant.m
% Author1:   Jeong-Hoon Ju       <jjh793012@naver.com>
% Author2:   Susana Lopez Moreno <susanalopezmoreno@pusan.ac.kr>
% Date:     2025-06-20
% Affiliation: Department of Mathematics, Pusan National University
% 
% Description:
%   Experiment checking the inequality of the determinantal identity
%   of the tensor geometric mean

%%%%%%%%%%%%%%%%%%%%%%%
% Hyperdeterminant of A
%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(2,2,2);
A(:,:,1) = [ 2  1
             1  2 ];
A(:,:,2) = [ 4  1
             1  4 ];
hyper_A=hyperdeterminant(A);

%%%%%%%%%%%%%%%%%%%%%%%
% Hyperdeterminant of B
%%%%%%%%%%%%%%%%%%%%%%%
B = zeros(2,2,2);
B(:,:,1) = [ 2  -1
             -1  2 ];
B(:,:,2) = [ 2  -1
             -1  2 ];
hyper_B=hyperdeterminant(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hyperdeterminant of G=A#_{*_M}B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G=zeros(2,2,2);
G(:,:,1) = matrix_gm(A(:,:,1),B(:,:,1));
G(:,:,2) = matrix_gm(A(:,:,2),B(:,:,2));
hyper_G=hyperdeterminant(G);

% Display results
disp(hyper_G);
disp(root);
if hyper_G~=root
    disp('not equal');
end