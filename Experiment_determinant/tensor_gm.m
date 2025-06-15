function G = tensor_gm(A,B)
G=zeros(2,2,2);
G(:,:,1) = matrix_gm(A(:,:,1),B(:,:,1));
G(:,:,2) = matrix_gm(A(:,:,2),B(:,:,2));
end