"""
    function IF(C, t, x=0, order=4, window=0)
Calculate instantaneous frequencies
"""
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
        H1 = Spline1D(t, C[:,i], k = order)
        Phi[:,i] = abs(hilbert(H1[x]))
        Phi[:,i] = Phi[:,i].*window
    end

	return Phi

end
