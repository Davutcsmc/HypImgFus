function [coeffM,paydaNorm]=funcCoeffMat(resR)
if isequal(resR,0.5)
    coeffM = 1./[ sqrt(2)/2 sqrt(2)/2 ; sqrt(2)/2  sqrt(2)/2  ];
elseif isequal(resR,1/3)
%     sig = (1/(2*(2.7725887)/ratio^2))^0.5;
%     sig = 1.2740;
%     coeffM = fspecial('gaussian',[3 3],sig);
    coeffM = [  0.0885 0.1205 0.0885 ;
                0.1205 0.1639 0.1205 ;
                0.0885 0.1205 0.0885 ];
elseif isequal(resR,0.25)
    coeffM = 1./[  3*sqrt(2)/2 sqrt(10)/2 sqrt(10)/2 3*sqrt(2)/2 ;
                   sqrt(10)/2  sqrt(2)/2  sqrt(2)/2  sqrt(10)/2 ;
                   sqrt(10)/2  sqrt(2)/2  sqrt(2)/2  sqrt(10)/2 ;
                   3*sqrt(2)/2 sqrt(10)/2 sqrt(10)/2 3*sqrt(2)/2 ];
elseif isequal(resR,0.125)
    coeffM =1./[sqrt(98) sqrt(74) sqrt(58) sqrt(50) sqrt(50) sqrt(58) sqrt(74) sqrt(98) ;
                sqrt(74) sqrt(50) sqrt(34) sqrt(26) sqrt(26) sqrt(34) sqrt(50) sqrt(74) ;
                sqrt(58) sqrt(34) sqrt(18) sqrt(10) sqrt(10) sqrt(18) sqrt(34) sqrt(58) ;
                sqrt(50) sqrt(26) sqrt(10) sqrt(02) sqrt(02) sqrt(10) sqrt(26) sqrt(50) ;
                sqrt(50) sqrt(26) sqrt(10) sqrt(02) sqrt(02) sqrt(10) sqrt(26) sqrt(50) ;
                sqrt(58) sqrt(34) sqrt(18) sqrt(10) sqrt(10) sqrt(18) sqrt(34) sqrt(58) ;
                sqrt(74) sqrt(50) sqrt(34) sqrt(26) sqrt(26) sqrt(34) sqrt(50) sqrt(74) ;
                sqrt(98) sqrt(74) sqrt(58) sqrt(50) sqrt(50) sqrt(58) sqrt(74) sqrt(98) ];
end
paydaNorm = sum(sum(coeffM));

end

