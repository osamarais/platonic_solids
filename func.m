% This is just here to evaluate the peaks function in the form that the
% Nelder Mead code asks for

function z = func(x)
z = peaks(x(1),x(2));
end