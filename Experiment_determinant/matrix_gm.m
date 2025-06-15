function[M]=matrix_gm(M1,M2)
M=sqrtm(M1)*sqrtm(inv(sqrtm(M1))*M2*inv(sqrtm(M1)))*sqrtm(M1);