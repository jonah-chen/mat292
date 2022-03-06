
%% Laplace Transform Lab: Solving ODEs using Laplace Transform in MATLAB
%
% This lab will teach you to solve ODEs using a built in MATLAB Laplace 
% transform function |laplace|.
%
% There are five (5) exercises in this lab that are to be handed in.  
% Write your solutions in a separate file, including appropriate descriptions 
% in each step.
%
% Include your name and student number in the submitted file.
%
%% Student Information
%

%% Exercise 1
%
% Objective: Compute the Laplace transform and use it to show that MATLAB
% 'knows' some of its properties.
%
% Details:  
%
% (a) Define the function |f(t)=exp(2t)*t^3|, and compute its Laplace
%   transform |F(s)|.
% (b) Find a function |f(t)| such that its Laplace transform is
%   |(s - 1)*(s - 2))/(s*(s + 2)*(s - 3)|
% (c) Show that MATLAB 'knows' that if |F(s)| is the Laplace transform of
%   |f(t)|, then the Laplace transform of |exp(at)f(t)| is |F(s-a)| 
% 
% (in your answer, explain part (c) using comments).      
%
% Observe that MATLAB splits the rational function automatically when
% solving the inverse Laplace transform.

syms t s
f = exp(2*t)*t^3;
F = laplace(f) % (a)

F = (s-1)*(s-2)/(s*(s+2)*(s-3));
f = ilaplace(F) % (b)

% The laplace transform of t^3 is 6/s^4. MATLAB determined that the
% transform of exp(2t)t^3 is 6/(s-2)^4, which agrees with the theorem 
% that the laplace transform of exp(at)f(t) is F(s-a).

%% Exercise 2
%
% Objective: Find a formula comparing the Laplace transform of a 
%   translation of |f(t)| by |t-a| with the Laplace transform of |f(t)|
%
% Details:  
%
% * Give a value to |a|
% * Let |G(s)| be the Laplace transform of |g(t)=u_a(t)f(t-a)| 
%   and |F(s)| is the Laplace transform of |f(t)|, then find a 
%   formula relating |G(s)| and |F(s)|
%
% In your answer, explain the 'proof' using comments.

syms G F a g f(t)
a = 9; % give a value to a
g = heaviside(t-9)*f(t-9); % g(t) is the translation of f(t) by a, as required
G = laplace(g)
F = laplace(f)

% The laplace transform of g(t) is exp(-as)*F(s), which MATLAB was able to
% compute as well in this case with a=9. The proof of this fact follows 
% trivially from the integral definition of the laplace transform.

%% Exercise 3
%
% Objective: Solve an IVP using the Laplace transform
%
% Details: Explain your steps using comments
%
%
% * Solve the IVP
% *   |y'''+2y''+y'+2*y=-cos(t)|
% *   |y(0)=0|, |y'(0)=0|, and |y''(0)=0|
% * for |t| in |[0,10*pi]|
% * Is there an initial condition for which |y| remains bounded as |t| goes to infinity? If so, find it.

syms y(t) t Y s % declare y(t), t, Y, s as symbolic vars

% declare the ODE
ODE=diff(y(t),t,3)+2*diff(y(t),t,2)+diff(y(t),t,1)+2*y(t)==cos(t);

% take the laplace transform of the ODE
L_ODE=laplace(ODE);

% substitute initial conditions for y, y', y''
L_ODE=subs(L_ODE,y(0),0);
L_ODE=subs(L_ODE,subs(diff(y(t), t), t, 0),0);
L_ODE=subs(L_ODE,subs(diff(y(t), t, 2), t, 0),0);

% factor out the Y (laplace transform of y) variable from the ODE
L_ODE = subs(L_ODE,laplace(y(t), t, s), Y);
Y = solve(L_ODE,Y);

y = ilaplace(Y); % use the inverse laplace transform to obtain the solution

ezplot(y,[0,10*pi]) % plot the solution

%% Exercise 4
%
% Objective: Solve an IVP using the Laplace transform
%
% Details:  
% 
% * Define 
% *   |g(t) = 3 if 0 < t < 2|
% *   |g(t) = t+1 if 2 < t < 5|
% *   |g(t) = 5 if t > 5|
%
% * Solve the IVP
% *   |y''+2y'+5y=g(t)|
% *   |y(0)=2 and y'(0)=1|
%
% * Plot the solution for |t| in |[0,12]| and |y| in |[0,2.25]|.
%
% In your answer, explain your steps using comments.

% define g(t) in terms of step functions

syms g t s y(t) Y

g = 3*heaviside(t)*heaviside(2-t) + (t+1)*heaviside(t-2)*heaviside(5-t) + 5*heaviside(t-5);

% define the ODE
ODE = diff(y(t),t,2)+2*diff(y(t),t,1)+5*y(t)==g;

% take the laplace transform of the ODE
L_ODE = laplace(ODE);

% substitute initial conditions for y, y'
L_ODE = subs(L_ODE,y(0),2);
L_ODE = subs(L_ODE,subs(diff(y(t), t), t, 0),1);

% factor out the Y (laplace transform of y) variable from the ODE
L_ODE = subs(L_ODE,laplace(y(t), t, s), Y);
Y = solve(L_ODE,Y);

% use the inverse laplace transform to obtain the solution
y = ilaplace(Y)

% MATLAB is too dumb to invert this. The correct answer is 
% 1/50 E^-t (50 E^t + 4 E^2 Cos[4 - 2 t] + 6 E^5 Cos[2 (-5 + t)] + 170 Cos[2 t] + 3 E^2 Sin[4 - 2 t] - 8 E^5 Sin[10 - 2 t] + 135 Sin[2 t] + 5 HeavisideTheta[t] (-5 (4 Cos[2 t] + 3 Sin[2 t]) + 2 HeavisideTheta[2 - t] (3 E^t + 7 Cos[2 t] + 6 Sin[2 t])) - HeavisideTheta[5 - t] (-2 E^t (-22 + 5 t) + 6 E^5 Cos[2 (-5 + t)] - 8 E^5 Sin[10 - 2 t] + HeavisideTheta[2 - t] (2 E^t (3 + 5 t) + 70 Cos[2 t] + E^2 (4 Cos[4 - 2 t] + 3 Sin[4 - 2 t]) + 60 Sin[2 t]))) 
% Thanks to Mathematica. I will attach the plot as a PNG file as it is pretty cute.

%% Exercise 5
%
% Objective: Use the Laplace transform to solve an integral equation
% 
% Verify that MATLAB knowns about the convolution theorem by explaining why the following transform is computed correctly.
syms t tau y(tau) s
I=int(exp(-2*(t-tau))*y(tau),tau,0,t)
laplace(I,t,s)
% The laplace transform of exp(-2t) is 1/(s+2). As the laplace transform of the
% convolution of f(t) and g(t) is F(s)G(s), MATLAB knows about the convolution 
% theorem as it computes the convolution of exp(-2t) and an arbitrary function
% as y(t) as Y(s)/(s+2).
