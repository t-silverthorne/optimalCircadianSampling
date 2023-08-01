using JuMP, Pajarito, HiGHS, Hypatia
import SCS
import LinearAlgebra
import Random
#%# 
M    = Int(1e2); # num of points in fine grid
N    = Int(24);  # num of measurements to take
eye  = Matrix{Float64}(LinearAlgebra.I, 2,2);

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
set_attribute(model, "time_limit", 60)
X    = randn(M,2); # TODO make this actually the design matrix
@variable(eOpt, t);
@variable(eOpt, mu[i=1:Int(M)],Bin);
@constraint(eOpt, sum(mu)==N);
@constraint(
      eOpt,
      X'*LinearAlgebra.diagm(0=>mu)*X - (t.*eye) >=0,
      PSDCone(),
);
@objective(eOpt, Max, t)

# solve
optimize!(eOpt)
