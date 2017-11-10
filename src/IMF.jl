abstract type TimeSeriesModel end

mutable struct EMD <: TimeSeriesModel
    k::Int
    C::Array{Float64,2}
    var::Array{Float64,1}
    maxorder::Int
    tols=Tuple{Float64,Float64}
    function new(k, maxorder=4, tols=(0.01,0.01))
        EMD(k, Array{Float64,2}(0,0), Array{Float64,1}(0), maxorder, tols)
    end
end

#Function to calculate intrinsic mode functions
#N is the maximum number of modes to find
#"""
#    IMF(y, t; toldev=0.01, tolzero = 0.01, order=4, N=5, window=0)
#Calculate the intrinsic mode functions of the sequence y along timespan t.
#"""

function fit!(model::EMD, y)
    t = ####
    model.C,model.var 
end

function IMF!(y, t; toldev=0.01, tolzero = 0.01, order=4, N=5, window=ones(length(t)))
	n = length(y)
	f = zeros(n,N)
	tempy = copy(y)
	eps = 0.00001
    model
	n_modes = 0;
    var = zeros(N)

	for i = 1:N
		avg = zeros(n,1)+1
        var[i] = 2*toldev

		while(mean(abs.(avg))>tolzero && sd > toldev)

			# Interpolate a spline through the maxima and minima
			max_ar, min_ar, tmax, tmin = findExtrema(tempy, t)

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

            var[i] = sqrt(mean( (avg.^2)./((y-f[:,i]).^2 + eps) ))

			f[:,i] = f[:,i] + avg
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
    
	return C,var
end
