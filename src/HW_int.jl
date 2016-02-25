
module HW_int

	# question 1 a)

# Demand function equation

function dem_func(p)
q = 2p^-0.5
retun q
end


# Plotting the demand function
# The solution corresponds to the area between the lines p = 1, p = 4, q = 1 and q = 2

using Gadfly
plot(
q,0,10,
Geom.hline(size=0.25mm, color="black"),
yintercept=[q(1),q(4)],
Stat.xticks(ticks=[0,1,4]),
Stat.yticks(ticks=[0, q(1), q(4)]),
Guide.ylabel("Quantity (q)"),
Guide.xlabel("Price (p)"),
Guide.title("Demand Function")

# The integration of q(p) = 2*(p^(-0.5)) on Wolfram Alpha from 1 to 4 gives 4
# (I could not find the command line for doing an integration in Julia)




	# question 1 b)
	# here are the packages I used

	using FastGaussQuadrature
	using Roots
	using Sobol
	using PyPlot
	using Distributions

	# here are some functions I defined for usage
	# in several sub questions

	# demand function
	function dem_func(p)
	q = 2p^-0.5
	return q
	end

	# gauss-legendre adjustment factors for map change
	# The suggested transformation is (not sure if this LaTex will be readible):
	#"$$ \\int_a^b f(x)\\,dx = \\frac{b-a}{2} \\int_{-1}^1 f\\left(\\frac{b-a}{2}x + \\frac{a+b}{2}\\right)dx $$\n"
	# with f(x) = q(p) , a = 1 and b = 4, we get :
	#"$$ \\int_1^4 q(p)\\,dp = \\frac{4-1}{2} \\int_{-1}^1 q\\left(\\frac{4-1}{2}p + \\frac{1+4}{2}\\right)dp $$\n"
	#"$$ \\int_1^4 q(p)\\,dp = \\frac{3}{2} \\int_{-1}^1 q\\left(\\frac{3}{2}p + \\frac{5}{2}\\right)dp $$\n"
	#"$$ \\int_1^4 q(p)\\,dp = \\frac{3}{2} \\int_{-1}^1 q\\left(\\frac{3}{2}p + \\frac{5}{2}\\right)dp $$\n"



	# eqm condition for question 2
	# this is the equilibrium condition: total demand = supply,


	# weighted sum for integration from the slides.



# QUESTION B. a.

function question_1b(n)

nodes, weights = gausslegendre(n)
function GL(n)
GL(n) = gausslegendre(n)
end
plot(weights, 3*((1.5nodes+2.5)^(-2))) # expression of the demand function when the domain is [-1,1]
result_1b = GL(n)
println(result_1b)

end

# QUESTION B. b.

function question_1c(n)

function
MCIntegrate(dem_func, 1, 4, n) #lower and upper bounds: 1 and 4; n points
end
result_c = MCIntegrate(dem_func, 1, 4, n)
println(result_1c)

end

# QUESTION B. c.
	function question_1d(n)

println("could not find it")

	end












	function question_2a(n)

# "Supply = Demand" implies "S = D + X"

function
f(log_theta_1, log_theta_2, p)=(exp(log_theta_1))*(1/p)+(exp(log_theta_2))*(p^(-0.5)) - 2
end

Sigma = [0.02 0.01; 0.01 0.01] # Defines the variances and covariances for log_theta_1 and log_theta_2
GH=gausshermite(10)
grid=kron(GH,GH)
result_1 = mean(p[,grid])
result_2 = var(p[,grid])
println(result_1,result_2)

	end


	function question_2b(n)

	function MCIntegrate(f, 1, 4, n)
	result_2b = MCIntegrate(f, 1, 4, n)
	println(result_2b)


	end


	# function to run all questions
	function runall(n=10)
		println("running all questions of HW-integration:")
		println("results of question 1:")
		question_1b(n)	# make sure your function prints some kind of result!
		question_1c(n)
		question_1d(n)
		println("")
		println("results of question 2:")
		q2 = question_2a(n)
		println(q2)
		q2b = question_2b(n)
		println(q2b)
		println("end of HW-integration")
	end

end
