%% ODE Lab: Creating your own ODE solver in MATLAB
% In this lab, you will write your own ODE solver for the Improved Euler method 
% (also known as the Heun method), and compare its results to those of |ode45|.
% 
% You will also learn how to write a function in a separate m-file and execute 
% it.
% 
% Opening the m-file lab3.m in the MATLAB editor, step through each part using 
% cell mode to see the results. Compare the output with the PDF, which was generated 
% from this m-file.
% 
% There are six (6) exercises in this lab that are to be handed in on the due 
% date. Write your solutions in the template, including appropriate descriptions 
% in each step. Save the .m files and submit them online on Quercus.
%% Student Information

%% Creating new functions using m-files.
% Create a new function in a separate m-file:
% 
% Specifics: Create a text file with the file name f.m with the following lines 
% of code (text):
%%
% 
%  function y = f(a,b,c) 
%  y = a+b+c;
%
%% 
% Now MATLAB can call the new function f (which simply accepts 3 numbers and 
% adds them together). To see how this works, type the following in the matlab 
% command window: sum = f(1,2,3)
%% Exercise 1
% Objective: Write your own ODE solver (using the Heun/Improved Euler Method).
% 
% Details: This m-file should be a function which accepts as variables (t0,tN,y0,h), 
% where t0 and tN are the start and end points of the interval on which to solve 
% the ODE, y0 is the initial condition of the ODE, and h is the stepsize. You 
% may also want to pass the function into the ODE the way |ode45| does (check 
% lab 2).
% 
% Note: you will need to use a loop to do this exercise. You will also need 
% to recall the Heun/Improved Euler algorithm learned in lectures. 
%% Exercise 2
% Objective: Compare Heun with |ode45|.
% 
% Specifics: For the following initial-value problems (from lab 2, exercises 
% 1, 4-6), approximate the solutions with your function from exercise 1 (Improved 
% Euler Method). Plot the graphs of your Improved Euler Approximation with the 
% |ode45| approximation.
% 
% (a) |y' = y tan t + sin t, y(0) = -1/2| from |t = 0| to |t = pi|
% 
% (b) |y' = 1 / y^2 , y(1) = 1| from |t=1| to |t=10|
% 
% (c) |y' = 1 - t y / 2, y(0) = -1| from |t=0| to |t=10|
% 
% (d) |y' = y^3 - t^2, y(0) = 1| from |t=0| to |t=1|
% 
% Comment on any major differences, or the lack thereof. You do not need to 
% reproduce all the code here. Simply make note of any differences for each of 
% the four IVPs.

a0 = heun(@(t,y) y.*tan(t) + sin(t), 0, pi, -1/2, 0.01);
a1 = ode45(@(t,y) y.*tan(t) + sin(t), [0, pi], -1/2);
b0 = heun(@(t,y) y.^-2, 1, 10, 1, 0.01);
b1 = ode45(@(t,y) y.^-2, [1,10], 1);
c0 = heun(@(t,y) 1 - t.*y./2, 0, 10, -1, 0.01);
c1 = ode45(@(t,y) 1 - t.*y./2, [0,10], -1);
d0 = heun(@(t,y) y.^3 - t.^2, 0, 1, 1, 0.01);
d1 = ode45(@(t,y) y.^3 - t.^2, [0, 1], 1);
%%
% (a)
plot(a0.t, a0.y, a1.x, a1.y);
xlabel('t');
ylabel('y');
title('Equation a');
legend('heun','ode45')
%%
% (b)
plot(b0.t, b0.y, b1.x, b1.y);
xlabel('t');
ylabel('y');
title('Equation b');
legend('heun','ode45')
%%
% (c)
plot(c0.t, c0.y, c1.x, c1.y);
xlabel('t');
ylabel('y');
title('Equation c');
legend('heun','ode45')
%%
% (d)
plot(d0.t, d0.y, d1.x, d1.y);
xlabel('t');
ylabel('y');
title('Equation d');
legend('heun','ode45')
%%
% There are no major differences between the two methods for 
% (a) and (b). There is a small difference in (c), where the 
% heun solution is smoother than the ode45 solution, likely
% due to the fact it is using a smaller step size. For (d),
% The numerical solution produced by heun has infinity for
% the last entry, whereas ode45 stopped before reaching 
% infinity and produced a warning as it cannot meet the 
% tolerance.
%% Exercise 3
% Objective: Use Euler's method and verify an estimate for the global error.
% 
% Details: 
% 
% (a) Use Euler's method (you can use euler.m from iode) to solve the IVP
% 
% |y' = 2 t sqrt( 1 - y^2 ) , y(0) = 0|
% 
% from |t=0| to |t=0.5|.

t = linspace(0,0.5,101);
f = @(t,y) 2*t.*sqrt(1-y.^2); %f=y'
y = euler(f, 0, t);
%% 
% (b) Calculate the solution of the IVP and evaluate it at |t=0.5|.

% The exact solution is y=sin(t^2)
sin(0.5^2)
%% 
% (c) Read the attached derivation of an estimate of the global error for Euler's 
% method. Type out the resulting bound for En here in a comment. Define each variable.

% E_n <= \Delta t (exp(M\Delta t n)-1)(1+M)/2
% where:
%   M         :     any positive real that is greater |f| and the norm of both
%                   partial derivatives of f (w.r.t y and t)
%   \Delta t  :     the step size
%   n         :     the step

% For this solution, the error at 0.5 is E_101. So,
%% 
% (d) Compute the error estimate for |t=0.5| and compare with the actual error.

deltaT = 0.5 / 101;
n = 101;

% Compute the derivatives
Dyf=-2*t.*y./sqrt(1-y.^2);
Dtf=2*sqrt(1-y.^2);

M = max([f(t,y), Dyf, Dtf])
err_est = (1+M)*deltaT/2 * (exp(M*deltaT*n)-1)
realerr = abs(y(101)-sin(0.5^2))
%% 
% (e) Change the time step and compare the new error estimate with the actual 
% error. Comment on how it confirms the order of Euler's method.

est_errs = zeros(30,1);
realerrs = zeros(30,1);
delta_Ts = linspace(0.00001, 0.01, 30);
for i = 1:30
    deltaT = delta_Ts(i);
    n = 0.5/deltaT;
    t = 0:deltaT:0.5;
    y = euler(f,0,t);
    Dyf=-2*t.*y./sqrt(1-y.^2);
    Dtf=2*sqrt(1-y.^2);
    M = max([f(t,y), Dyf, Dtf]);
    est_errs(i) = (1+M)*deltaT/2 * (exp(M*deltaT*n)-1);
    realerrs(i) = abs(y(end)-sin(0.5^2));
end

plot(delta_Ts, realerrs, 'o', delta_Ts, est_errs, 'o')

% From the plot, the estimated error is a straight line. This can be seen in 
% the formula as everything is constant except for \Delta t (note that 
% \Delta t n=0.5 always).
% The real errors seem to fit to a straight line, but gets noisy for 
% larger values of \Delta t although the plot still suggests that the
% order of euler's method is O(\Delta t).
%% Adaptive Step Size
% As mentioned in lab 2, the step size in |ode45| is adapted to a specific error 
% tolerance.
% 
% The idea of adaptive step size is to change the step size |h| to a smaller 
% number whenever the derivative of the solution changes quickly. This is done 
% by evaluating f(t,y) and checking how it changes from one iteration to the next.
%% Exercise 4
% Objective: Create an Adaptive Euler method, with an adaptive step size |h|.
% 
% Details: Create an m-file which accepts the variables |(t0,tN,y0,h)|, as in 
% exercise 1, where |h| is an initial step size. You may also want to pass the 
% function into the ODE the way |ode45| does.
% 
% Create an implementation of Euler's method by modifying your solution to exercise 
% 1. Change it to include the following:
% 
% (a) On each timestep, make two estimates of the value of the solution at the 
% end of the timestep: |Y| from one Euler step of size |h| and |Z| from two successive 
% Euler steps of size |h/2|. The difference in these two values is an estimate 
% for the error.
% 
% (b) Let |tol=1e-8| and |D=Z-Y|. If |abs(D)<tol|, declare the step to be successful 
% and set the new solution value to be |Z+D|. This value has local error |O(h^3)|. 
% If |abs(D)>=tol|, reject this step and repeat it with a new step size, from 
% (c).
% 
% (c) Update the step size as |h = 0.9*h*min(max(tol/abs(D),0.3),2)|.
% 
% Comment on what the formula for updating the step size is attempting to achieve.

% abs(D) is the estimated error of that step
% tol is the minimum acceptable tolerance.
% tol/abs(D) is the ratio between the minimum acceptable tolerance and
% the estimated error of the step.
% 
% The formula first reduces the step-size h to 90% of what it was.
% Then, depending on the ratio between the minimum acceptable tolerance and 
% the estimated error, multiply the step size by the ratio. This makes sense  
% because if the estimated error is too small, then the step size should be 
% increased to improve efficiency. But, if the estimated error is too large, 
% the step size should be decreased to to improve accuracy.
%
% However, changing the step size by too much at once may not be good as 
% firstly, the local error is O(h^3) not O(h), and secondly, the error abs(D)
% is only an estimate based on the computed Z and Y. 
%
% Hence, this factor is limited in the range of 0.3 to 2 to ensure the step 
% does not change too rapidly.
%% Exercise 5
% Objective: Compare Euler to your Adaptive Euler method.
% 
% Details: Consider the IVP from exercise 3.
% 
% (a) Use Euler method to approximate the solution from |t=0| to |t=0.75| with 
% |h=0.025|.

te = 0:0.025:0.75;
ye = euler(f,0,te);
%% 
% (b) Use your Adaptive Euler method to approximate the solution from |t=0| 
% to |t=0.75| with initial |h=0.025|.

adaptive = adap_euler(f,0,0.75,0,0.025);
%% 
% (c) Plot both approximations together with the exact solution.

plot(te, ye, adaptive.t, adaptive.y, te, sin(te.^2))
legend('euler', 'adaptive euler', 'exact',"Location","eastoutside")
%% Exercise 6
% Objective: Problems with Numerical Methods.
% 
% Details: Consider the IVP from exercise 3 (and 5).
% 
% (a) From the two approximations calculated in exercise 5, which one is closer 
% to the actual solution (done in 3.b)? Explain why.

% The one using the adaptive method is closer to the actual solution. 
% The reason is that using a step size of 0.025 will result in an estimated
% error of far larger than 1e-8 near 0.075. 
adaptive.t(end)-adaptive.t(end-1)
% In fact, last step of the adaptive solution is far smaller than 0.025. 
% It is 0.00027.
%% 
% (b) Plot the exact solution (from exercise 3.b), the Euler's approximation 
% (from exercise 3.a) and the adaptive Euler's approximation (from exercise 5) 
% from |t=0| to |t=1.5|.

te = 0:0.025:1.5;
ye = euler(f,0,te);
adaptive = adap_euler(f,0,1.5,0,0.025);
plot(te, ye, adaptive.t, adaptive.y, te, sin(te.^2));
legend('euler', 'adaptive euler', 'exact', "Location","eastoutside")
%% 
% (c) Notice how the exact solution and the approximations become very different. 
% Why is that? Write your answer as a comment.

% Taking a closer inspection of the ODE in exercise 3, y' is real if and only if t=0 or |y|<=1. 
% Note that for the exact solution, y=sin(t^2), y(1)=1. However, if the
% numerical method steps over y=1 by even a tiny bit, y' becomes a complex 
% number. MATLAB gives an error when plotting because of this.

% Plotting the imaginary parts of the numerical solutions confirms this.
plot(adaptive.t, imag(adaptive.y))
title('imagnary part of adaptive euler solution')
plot(te, imag(ye))
title('imagnary part of non-adaptive euler solution')

% As the imaginary part of y is very large, when squared it is capable of
% significantly altering the value of y', which is used in euler's method. 
% Hence, it causes the approximation to deviate from the exact solution.
