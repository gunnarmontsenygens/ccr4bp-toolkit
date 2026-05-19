function X_L_vec = H2L(X_H_vec)

    X_L_vec = [X_H_vec(1);
               X_H_vec(2);
               X_H_vec(3) + X_H_vec(2);
               X_H_vec(4) - X_H_vec(1)];
end