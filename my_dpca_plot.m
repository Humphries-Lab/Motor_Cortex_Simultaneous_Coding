function my_dpca_plot(Xfull, W, V)
% produces a plot of the dPCA results. X is the data matrix, W and V
% are decoder and encoder matrices

X = Xfull(:,:)';
Xcen = bsxfun(@minus, X, mean(X));
Z = Xcen * W;
b = V(:,1:numCompToShow)'*V(:,1:numCompToShow);
end