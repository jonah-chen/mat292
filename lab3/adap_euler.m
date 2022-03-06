function out = adap_euler(odefun,t0,tN,y0,h)
    % Declare constants
    tol = 1e-8;
    
    % declare output
    out.t = [t0];
    out.y = [y0];
    
    while (out.t(end) < tN)
        nsuccess = true;
        while (nsuccess)
            % Make Y: 1 euler step size h
            Y = out.y(end) + h*odefun(out.t(end),out.y(end));

            % Make Z: 2 euler steps size h/2
            Z = out.y(end) + h/2*odefun(out.t(end),out.y(end));
            Z = Z + h/2*odefun(out.t(end)+h/2, Z);
            
            % Make D
            D = Z-Y;
            
            if (abs(D)<tol) % accept solution
                out.t = [out.t, out.t(end)+h];
                out.y = [out.y, Z-D];
                nsuccess = false;
            end

            % update h
            h = 0.9*h*min(max(tol/abs(D),0.3),2);
        end
    end
end