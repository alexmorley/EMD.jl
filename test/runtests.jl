# Testing Script for EMD module
using Base.Test
using EMD

println("EMDUtil test started")
#test case for findExtrema function
y = Float64[1,3,1,2]
t = Float64[1,2,3,4]
max, min, tmax, tmin = findExtrema(y, t)

@test max == [1,3,2]

@test tmax == [1,2,4]

@test min == [1,1,2]

@test tmin == [1,3,4]
