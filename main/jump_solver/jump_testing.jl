using JuMP, Pajarito, HiGHS, Hypatia
import SCS
import LinearAlgebra
import Random

freq = 3.6;
M    = Int(1e2); # num of points in fine grid
N    = Int(8);  # num of measurements to take
 
function getLambdaEoptMethod(freq,M,N)
  # will need as input A,f,M,N
  eye  = Matrix{Float64}(LinearAlgebra.I, 2,2);
  tau  = LinRange(0,1,Int(M)+1);
  tau  = tau[1:Int(M)]

  # formulate E-optimality as SDP
  eOpt = Model(
      optimizer_with_attributes(
          Pajarito.Optimizer,
          "oa_solver" => optimizer_with_attributes(
              HiGHS.Optimizer,
              MOI.Silent() => true,
              "mip_feasibility_tolerance" => 1e-8,
              "mip_rel_gap" => 1e-6,
          ),
          "conic_solver" =>
              optimizer_with_attributes(Hypatia.Optimizer, 
                                        MOI.Silent() => true),
      )
  )
  X    = [cos.(2pi*freq*tau) sin.(2pi*freq*tau)];
  
  @variable(eOpt, t);
  @variable(eOpt, mu[i=1:Int(M)],Bin);
  @constraint(eOpt, sum(mu)==N);
  @constraint(
        eOpt,
        X'*LinearAlgebra.diagm(0=>mu)*X - (t.*eye) >=0,
        PSDCone(),
  );
  @objective(eOpt, Max, t)

  optimize!(eOpt)
  return eOpt
end
