
module HW_int


	# question 1 b) 

	using FastGaussQuadrature
	using Roots
	using Sobol
	using PyPlot
	using Distributions

	# demand function
	q(p) = 2*(p.^-0.5)

	# gauss-legendre adjustment factors for map change
	ba2(lb,ub) = (ub-lb)/2
	ab2(lb,ub) = (lb+ub)/2

	# eqm condition for question 2
	# this is the equilibrium condition: total demand = supply, 
	# i.e. domestic + export demand = 2
	function dd(p,t1,t2)
	    exp(t1)/p .+ exp(t2) .* p^-0.5 - 2
	end

	function fzero_wrap(f,lb,ub) 
		try
			fzero(f,lb,ub)
		catch
			println("no eqm price found")
			return(NaN)
		end
	end

	function question_1b(n)

		gl = gausslegendre(n);

		# bounds of integration

		# function to integrate:
		# CS_4 = \int_0^4 q(p) dp - 4
		# CS_1 = \int_0^1 q(p) dp - 2
		# CS_4 - CS_1 = \int_0^4 q(p) dp - 4 - \int_0^1 q(p) dp + 2 
		#             = \int_1^4 q(p) dp - 2 

		# bounds on integration
		a = 1
		# b = pstar = 4
		pstar = 4

		# equilibrium quantity
		# qstar = q(pstar)

		# use transformation formula to map into [-1,1]
		pts  = ba2(a,pstar).*gl[1] .+ ab2(a,pstar)	# integration points
		vals = q(pts)  # function values at those points
 		Integ = ba2(a,pstar) * (gl[2]' * vals ) 
		Integ -= 2
		Integ = Integ[1]

		# plot
		figure()
		plot(pts,vals,"o-")
		title("Gauss Laguerre")

		println("estimated change in CS using $n gauss legendre nodes is $(Integ)")
		println("i.e. an error of $(round(abs(100*(2-Integ))/2,5)) percent")
		println("")
	end

	# question 1 c)

	function question_1c(n)

		# function to integrate:
		# CS_4 = \int_0^4 q(p) dp - 4
		# CS_1 = \int_0^1 q(p) dp - 2
		# CS_4 - CS_1 = \int_0^4 q(p) dp - 4 - \int_0^1 q(p) dp + 2 
		#             = \int_1^4 q(p) dp - 2 


		# get n random numbers from [1,4]
		pts = rand(n)*3 + 1	
		vals = q(pts)

		# integrate: Monte carlo is defined for the "hypercube" [0,1]
		# we need to adjust the "volume" of this cube to be 3
		Integ = 3*mean(vals) - 2
		
		# plot
		figure()
		plot(pts,vals,"o")
		axhline(mean(vals),color="red")
		title("Monte Carlo")

		println("estimated change in CS using $n monte carlo nodes is $Integ)")
		println("i.e. an error of $(round(abs(100*(2-Integ))/2,5)) percent")
		println("")

	end

	function question_1d(n)

		# CS1: p=4
		# ========

		pstar=4.0
		pqstar = pstar * q(pstar)

		s = SobolSeq(1,[1],[pstar])  # 1-dimensional sobol sequence in [1,pstar]
		pts = zeros(n)
		for i in 1:n
			pts[i] = next(s)[1]
		end

		vals = q(pts)

		# integrate: Monte carlo is defined for the hypercube [0,1]
		# we need to extend the length of that interval to be 3
		Integ = 3*mean(vals) - 2
		
		# plot
		figure()
		plot(pts,vals,"o")
		axhline(mean(vals),color="red")
		title("Quasi Monte Carlo")

		println("estimated change in CS using $n Quasi monte carlo nodes is $Integ)")
		println("i.e. an error of $(round(abs(100*(2-Integ))/2,5)) percent")
		println("")

	end

	# question 2

	function question_2a(n)

		gh = gausshermite(n)

		Sigma = hcat([0.02, 0.01],[0.01,0.01])
		Omega = chol(Sigma,Val{:U})

		mu = [0.0;0.0]

		# kronecker product of grids and weights
		gr = hcat(kron(ones(n),gh[1]),kron(gh[1],ones(n)))
		wt = kron(gh[2],gh[2]) / pi	# watch out for the pi!

		# make adjustment for correlation in shocks
		grids = Omega * gr'	 + zeros(2,n*n)   # zeros here would be a matrix with mu

		# find eqm price at each combination of theta1,theta2
		pstar = zeros(n*n)
		for i in 1:length(pstar)
			pstar[i] = fzero(x->dd(x,grids[1,i],grids[2,i]),0.001,15)
		end

		# plot
		figure()
		plot(grids[1,:],grids[2,:],"o")
		title("Question 2a: Gauss hermite theta grid")
		
		EP = dot(wt,pstar)
		VAR = dot(wt, (pstar .- EP).^2 )

		return Dict("E[p]"=>EP, "Var[p]"=>VAR)

	end

	function question_2b(n)

		# for fairness, let's also create n^2 points as in 2a
		n = n^2

		Sigma = hcat([0.02, 0.01],[0.01,0.01])
		M = MvNormal(Sigma)	# create mean zero joint normal distribution
		pts = rand(M,n)	# just draw from it randomly

		# find eqm price at each combination of theta1,theta2
		pstar = zeros(n)
		for i in 1:length(pstar)
			pstar[i] = fzero(x->dd(x,pts[1,i],pts[2,i]),0.01,10)
		end
		
		# plot
		figure()
		plot(pts[1,:],pts[2,:],"o")
		title("Question 2b: Monte Carlo theta grid")

		EP = mean(pstar)

		VAR = mean( (pstar .- EP).^2 )
		return Dict("E[p]"=>EP, "Var[p]"=>VAR)

	end

	function question_2bonus(n)

		# for fairness, let's also create n^2 points as in 2a
		n = n^2

		s = SobolSeq(2,[1,1],[4,4])  # 2-dimensional sobol sequence 
		pts = hcat([next(s) for i=1:n])
		println("here's the sobol sequence")
		println(pts)

		# find eqm price at each combination of theta1,theta2
		pstar = Float64[]
		for i in 1:n
			tmp = fzero_wrap(x->dd(x,pts[i][1],pts[i][2]),0.0001,20)
			if isnan(tmp)
				# nothing
			else
				push!(pstar,tmp)
			end
		end
		println("it's very hard to find the eqm price for those values")
		
		# plot
		figure()
		plot([pts[i][1] for i in 1:n],[pts[i][2] for i in 1:n],"o")
		title("Question 2bonus: Quasi Monte Carlo theta grid")

		EP = mean(pstar)

		VAR = mean( (pstar .- EP).^2 )
		return Dict("E[p]"=>EP, "Var[p]"=>VAR)

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
		println("")
		println("bonus question: Quasi monte carlo:")
		q2bo = question_2bonus(n)
		println(q2bo)
		println("end of HW-integration")
	end

end

