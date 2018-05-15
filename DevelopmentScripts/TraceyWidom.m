function [t,f2] = TraceyWidom(tLeft,tRight,dx)

% TRACEYWIDOM compute Tracey-Widom distribution for complex numbers (f2)
% [T,F] = TRACEYWIDOM(tL,tR,DX) computes the Tracey-Widom distribution
% between the points [tL,tR] - for the left and right-most parts, given
% discretisation DX.
%
% Returns: T: evalution points; F: the density
%
% From: Edelman & Wang 2013

opts = odeset('reltol',1e-12,'abstol',1e-15);

deq = @(t,y) [y(2); t*y(1)+2*y(1)^3; y(4); y(1)^2];

y0 = [airy(tRight); airy(1,tRight); 0; airy(tRight)^2];
[t,y] = ode45(deq,tRight:-dx:tLeft,y0,opts);
F2 = exp(-y(:,3));
f2 = gradient(F2,t);