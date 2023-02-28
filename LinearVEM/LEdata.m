function pde = LEdata

% Data for linear  elliptic pde with variable coefficients
% Author : Jai Tushar, BITS - Pilani


pde = struct('exactu',@exactu,'Dux',@Dux,'Duy',@Duy,'f',@f);
         
    function u = exactu(p)
        % Equivalent projectors Example
        u = sin(2*p(:,1)+0.5).*cos(p(:,2)+0.3)+log(1+p(:,1).*p(:,2));
    end
    
    function s = Dux(p)
        % Equivalent projectors example
        s = (p(:,2)./(p(:,1).*p(:,2)+1))+2*cos(2*p(:,1)+0.5).*cos(p(:,2)+0.3);
    end

    function t = Duy(p)
        % Equivalent projectors example
        t = (p(:,1)./(p(:,1).*p(:,2)+1))-sin(2*p(:,1)+0.5).*sin(p(:,2)+0.3);
    end

    function l = f(p)
        % Equivalent projectors example
        l = (((p(:,1).^2+p(:,2).^2)./(p(:,1).*p(:,2)+1).^2)+5*sin(2*p(:,1)+0.5).*cos(p(:,2)+0.3)) + exactu(p);
    end
end