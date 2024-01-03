function X = dummy(Mx,r)

% 조금더 고민해봐야함
[X,~,~] = svd(Mx);
X = X(:,1:r);

% Compute the Moore-Penrose pseudo-inverse of M
Mx_pseudo_inv = pinv(Mx);

% Compute Y using the formula
X = Mx_pseudo_inv * X * sqrtm(inv(X' * Mx_pseudo_inv * X));




% Check if Y'*M*Y is approximately equal to the identity matrix I
% identity_matrix = eye(n);
% error = norm(Y'*M*Y - identity_matrix);

% if error < 1e-6
%     disp('Y''*M*Y is approximately equal to the identity matrix I.');
% else
%     disp('Y''*M*Y is not approximately equal to the identity matrix I.');
% end

% Display the computed Y
% disp('Matrix Y:');
% disp(Y);


end







