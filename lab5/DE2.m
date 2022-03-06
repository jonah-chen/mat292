
function out = DE2(P,Q,G,t0,tN,y0,y1,h)
    out.x = t0:h:tN;
    out.y = zeros(size(out.x));
    out.y(1) = y0;
    out.y(2) = y0+h*y1;
    for i = 2:length(out.x)-1
        x = out.x(i);
        y0 = out.y(i-1);
        y1 = out.y(i);
        p = P(x);
        q = Q(x);
        g = G(x);

        out.y(i+1) = g*h^2 - y0 + h*p*y0 + 2*y1 - h*p*y1 - h^2*q*y1;
    end
end
