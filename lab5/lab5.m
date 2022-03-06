%% Second-Order Lab: Second-Order Linear DEs in MATLAB
%
%
% In this lab, you will learn how to use |iode| to plot solutions of 
% second-order ODEs. You will also learn to classify the behaviour of 
% different types of solutions.
%
% Moreover, you will write your own Second-Order ODE system solver, and 
% compare its results to those of |iode|.
%
% Opening the m-file lab5.m in the MATLAB editor, step through each
% part using cell mode to see the results.
%
% There are seven (7) exercises in this lab that are to be handed in on the
% due date of the lab.
%
%% Student Information
%

%% Exercise 1
%
% 100% of the solutions are decay while oscillating.
% The exact solution is y=(C1*cos(2t)+C2*sin(2t))*exp(-t/2), which is the product of 
% a decaying exponential term oscillating sinusoidal function. Hence, non-trivial solutions would
% decay while oscillating.
%

%% Exercise 2
%
% 100% of the solutions are grow. 
% Let g=sqrt(3)/2 which is less than 1, The exact solution is c1*exp((1+g)*t)+c2*exp((1-g)*t), which is
% two growing exponential terms. Hence, non-trivial solutions would always grow.

%% Exercise 3
%
% 100% of the solutions are decay. 
% Let h=sqrt(2)/2 which is less than g. The exact solution is c1*exp((g+h)*t)+c2*exp((g-h)*t), which is
% two decaying exponential terms. Hence, non-trivial solutions would always decay.

%% Exercise 4
%
roots([1 2 6 2 5])

%% 
% The general solution is y(t)= (C1*cos(2t)+C2*sin(2t))*exp(-t) + C3*cos(t) + C4*sin(t)
%
% 100% of the solutions are decay while oscillating, but never decays to zero.
% 0% of the solutions are decay while oscillating and decays to zero.
% 0% of the solutions just oscillate.

%% Exercise 5
%
% Objective: Classify equations given the roots of the characteristic
% equation.
%
% Details: Your answer can consist of just a short sentence, as |grows| or 
% |decays while oscillating|.
%
% Consider a second-order linear constant coefficient homogeneous DE with
% |r1| and |r2| as roots of the characteristic equation.
%
% Summarize your conclusions about the behaviour of solutions for randomly
% chosen initial data when.
%
% The general solution is Ae^r1t + Be^r2t.
% (a) |0 < r1 < r2|
% grows 100% of the time.
%
% (b) |r1 < 0 < r2|
% grows if B non-zero, which occurs 100% of the time for random initial data.
% decay if B zero, which occurs 0% of the time for random initial data.
%
% (c) |r1 < r2 < 0|
% decay 100% of the time.
% 
% (d) |r1 = alpha + beta i| and  |r2 = alpha - beta i| and |alpha < 0|
% shrinks while oscilating 100% of the time.
%
% (e) |r1 = alpha + beta i| and  |r2 = alpha - beta i| and |alpha = 0|
% oscillates forever 100% of the time.
%
% (f) |r1 = alpha + beta i| and  |r2 = alpha - beta i| and |alpha > 0|
% grows while oscilating 100% of the time.
%

%% Exercise 7
%

P=@(t) exp(-t/5);
Q=@(t) (1-exp(-t/5));
G=@(t) sin(2*t);

out = DE2(P,Q,G,0,20,1,0,0.1);
plot(out.x,out.y)
