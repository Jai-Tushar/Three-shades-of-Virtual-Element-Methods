function pde = LEdata

% Data for linear  elliptic pde with variable coefficients
% Author : Jai Tushar, BITS - Pilani


pde = struct('exactu',@exactu,'Dux',@Dux,'Duy',@Duy,'f',@f);
         
    function u = exactu(p)
        u = sin(p(:,1)*pi).*sin(p(:,2)*pi);
    end
    
    function s = Dux(p)
        s = pi*cos(p(:,1)*pi).*sin(p(:,2)*pi);
    end

    function t = Duy(p)
        t = pi*sin(p(:,1)*pi).*cos(p(:,2)*pi);
    end

    function l = f(p)
        l = 2*pi^2*sin(pi*p(:,1)).*sin(pi*p(:,2));
    end
end