%find the mode shape and resonant frequencies of the system
%INPUTS
%M_mat: the n x n mass (inertia) matrix
%K_mat: the n x n stiffness matrix
%OUTPUTS
%U_r: the mode shape vector of the system
%K_mat: the resonant frequency of the system
function [U_r, omega_r] = find_mode_shape_and_resonant_frequencies(M_mat, K_mat)
    %Use MATLAB to solve the generalized eigenvalue problem
    [U_r,lambda] = eig(K_mat,M_mat);
    omega_r = sqrt(diag(lambda));
end