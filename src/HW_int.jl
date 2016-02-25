
module HW_int


	# question 1 b) 
	# here are the packages I used

	using FastGaussQuadrature
	using Roots
	using Sobol
	using PyPlot
	using Distributions

	# here are some functions I defined for useage 
	# in several sub questions

	# demand function

	# gauss-legendre adjustment factors for map change

	# eqm condition for question 2
	# this is the equilibrium condition: total demand = supply, 


	# weighted sum for integration from the slides.



	function question_1b(n)

	end


	function question_1c(n)

	end


	function question_1d(n)

	end


	function question_2a(n)

	end


	function question_2b(n)

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

