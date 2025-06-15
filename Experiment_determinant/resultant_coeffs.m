function res = resultant_coeffs(a_1,a_2,a_3,b_1,b_2,b_3)
res = a_3^2*b_1^2-a_2*a_3*b_1*b_2+a_1*a_3*b_2^2+a_2^2*b_1*b_3-2*a_1*a_3*b_1*b_3-a_1*a_2*b_2*b_3+a_1^2*b_3^2;
end