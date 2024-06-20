function MEE = KEP2MEE(KEP)
a = KEP(:,1);
e = KEP(:,2);
i = KEP(:,3);
W = KEP(:,4);
w = KEP(:,5);
TA = KEP(:,6);

p = a .* (1 - e.^2);
f = e .* cos(w + W);
g = e .* sin(w + W);
h = tan(i/2) .* sin(W);
k = tan(i/2) .* cos(W);
L = W + w + TA;
MEE = [p f g h k L];
end

