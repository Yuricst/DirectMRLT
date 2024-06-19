%% Convert MEE to one-revolution R,V history
%
% Yuri Shimane, yuri.shimane@gatech.edu
%   Created    : 2024/06/19
%   Last edits : 2024/06/19
%
function RV = MEE2RVorbit(mu,MEE)
    MEEs = MEE .* ones(1000,6);
    MEEs(:,6) = linspace(0,2*pi,1000);
    RV = MEE2RV(mu,MEEs);
end

