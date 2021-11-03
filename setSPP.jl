# --------------------------------------------------------------------------- #

# Setting an ip model of SPP
function setSPP(C, A)
  m, n = size(A)
  model = Model(GLPK.Optimizer)
  #model = Model(with_optimizer(GLPK.Optimizer))
  #ip = Model(with_optimizer(solverSelected))
  @variable(model, x[1:n], Bin)
  @objective(model, Max, sum(C[j]x[j] for j = 1:n))
  @constraint(model, cte[i=1:m], sum(A[i,j] * x[j] for j=1:n) <= 1)
  optimize!(model)
  println("z = ",objective_value(model)) # affichage de la valeur optimale
  println("x = ",value.(model[:x]))
  #return objective_value(model)
end
