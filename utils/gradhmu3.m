function out = gradhmu3(x, reg_val)
mu = 1;
out = (x - proxmu2(x, reg_val)) * (1 / mu);
end


function out = proxmu2(x, reg_val)

N = length(x)/15;
lambda = reg_val;
w = 1;


X = reshape(x, 15, N);
% X_s=squeeze(sqrt(sum(X.^2,2)))+eps;
X_s = reshape(sqrt(sum(X.^2, 1)), 1,N) + eps;
w = repmat(w, 1,N);
q = max(X_s-w.*lambda, 0);
out = reshape(bsxfun(@times, reshape(q ./ X_s,1,N), X), 15,N);
indx = out(1:3,:) <= 0;
temp = out(1:3,:);
temp(indx > 0)=0;
out(1:3,:) = temp;
out = reshape(out,[],1);

end