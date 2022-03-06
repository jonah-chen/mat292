function out = solvesystem(f,g,t0,tN,x0,h)
    out.t = t0:h:tN; % vector of all t
    imax = size(out.t);
    out.x = zeros(2,imax(2));
    out.x(1,1) = x0(1);
    out.x(2,1) = x0(2);

    for i=1:imax(2)-1
        x1 = out.x(1,i) + h*f(out.t(i),out.x(1,i),out.x(2,i));
        x2 = out.x(2,i) + h*g(out.t(i),out.x(1,i),out.x(2,i));

        out.x(1,i+1) =  out.x(1,i) + h/2*(f(out.t(i),out.x(1,i),out.x(2,i))+f(out.t(i+1),x1,x2));
        out.x(2,i+1) =  out.x(2,i) + h/2*(g(out.t(i),out.x(1,i),out.x(2,i))+g(out.t(i+1),x1,x2));
    end
end
