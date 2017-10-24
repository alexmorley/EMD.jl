#function to calculate instantaneous frequencies
function IF(C, t, x=0, order=4, window=0)

	n,N = size(C)
	n = length(t)

	if x == 0
		x = t+1e-6
	end

	if window == 0
		window = ones(length(x))
	end

	nx = length(x)
	Phi = zeros(nx,N)

	for i =1:N
		H1 = Spline(C[:,i], t, order)

		Phi[:,i] = 1/(2pi)*(H1(x).*H1(x,1,true) - H1(x, 1).*H1(x,0,true))./(H1(x).^2 + H1(x,0,true).^2)
		Phi[:,i] = Phi[:,i].*window
	end

	return Phi

end
