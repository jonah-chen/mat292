function out = heun(odefun, t0, tN, y0, h)
    out.t = t0:h:tN; % vector of all t
    imax = size(out.t);
    out.y = zeros(imax);
    out.y(1) = y0;
    imax = imax(2);
    for i=1:imax-1
        y1 = out.y(i) + h*odefun(out.t(i),out.y(i));
        out.y(i+1) =  out.y(i) + h/2*(odefun(out.t(i),out.y(i))+odefun(out.t(i+1),y1));
    end
end