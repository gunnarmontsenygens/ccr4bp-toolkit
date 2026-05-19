function X_H_vec = L2H(X_L_vec)

    X_H_vec = [X_L_vec(1);
               X_L_vec(2);
               X_L_vec(3) - X_L_vec(2);
               X_L_vec(4) + X_L_vec(1)];
end