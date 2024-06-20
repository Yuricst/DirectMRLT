%% Convert Keplerian elements to Modified Equinoctial Elements
%
% Yuri Shimane, yuri.shimane@gatech.edu
%   Created    : 2024/06/20
%   Last edits : 2024/06/20
%
%
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
h = tan(i/2) .* cos(W);
k = tan(i/2) .* sin(W);
L = W + w + TA;
MEE = [p f g h k L];
end

