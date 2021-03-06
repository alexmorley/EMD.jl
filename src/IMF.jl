#Function to calculate intrinsic mode functions
#N is the maximum number of modes to find
"""
    IMF(y, t; toldev=0.01, tolzero = 0.01, order=4, N=5, window=0)
Calculate the intrinsic mode functions of the sequence y along timespan t.
"""
function IMF(y, t; toldev=0.01, tolzero = 0.01, order=4, N=5, window=0)
	if window==0
		window = ones(length(t))
	end

	n = length(y)
	f = zeros(n,N)
	tempy = copy(y)
	eps = 0.00001

	n_modes = 0;

	for i = 1:N
		avg = zeros(n,1)+1
		sd = 2*toldev

		while(mean(abs.(avg))>tolzero && sd > toldev)

			# Interpolate a spline through the maxima and minima
			max_ar, min_ar, tmax, tmin = findExtrema(tempy, t)

			# # Don't use too high an order to interpolate, restrict it if there are not many extrema
			# p_max = min(order, length(max_ar)-2)
			# p_min = min(order, length(min_ar)-2)

			# # Make even order
			# p_max = Int(2*floor(p_max/2))
			# p_min = Int(2*floor(p_min/2))

			# # At least linear
			# p_max = max(p_max,2)
			# p_min = max(p_min,2)
			p_max = (length(max_ar) >= order) ? order : 4
			p_min = (length(min_ar) >= order) ? order : 4

			while true
				p_min = p_min < length(min_ar) ? (p_min; break) : (p_min - 1)
				p_max = p_max < length(max_ar) ? (p_max; break) : (p_max - 1)
			end

			S1 = Spline1D(tmax, max_ar, k = p_max)
			S2 = Spline1D(tmin, min_ar, k = p_min)

			# Find mean of envelope
			avg = (S1(t) + S2(t)) / 2
			avg = avg.*window

			tempy = tempy-avg

			sd = mean( (avg.^2)./((y-f[:,i]).^2 + eps) )

			f[:,i] = f[:,i] + avg

			#println(sd)
		end

		tempy = copy(f[:,i])

		# Check to see if it's worth continuing or if the remainder is monotone
		c = mean(abs.(tempy))
		d = diff(tempy)
		n_modes = n_modes + 1
		if all(d+c*tolzero .> 0) || all(d-c*tolzero .< 0)
			break
		end
	end

	C = zeros(n,N)
	C[:,1] = y - f[:,1]
	for i = 2:n_modes
		C[:,i] = f[:,i-1]-f[:,i]
	end

	return C
end
