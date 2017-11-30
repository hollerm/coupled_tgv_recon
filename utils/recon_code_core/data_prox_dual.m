
function y = data_prox_dual(x,x0,ld,sig,type);


if strcmp(type,'gaussian')%L2 data fidelity

	y = (x - sig*x0) / (1+(sig/ld));

elseif strcmp(type,'poisson') %Poisson log-likelyhood

    y = x - (( x - ld + sqrt( abs(x-ld).^2 + 4*ld*sig*x0 ) )/ 2);
	
else
	error('Unknown noise type');
end
