using DelimitedFiles: readdlm

# Loader functions for controller matrices stored in the package's data/
# directory. The feedforward-generator matrices are computed by
# scripts/compute_feedforward.jl and consumed by CascadeFFDyadBot; the LQG
# controller matrices are computed by scripts/tune_lqg.jl and consumed by
# LQGControlledDyadBot.
const __DATA_DIR = joinpath(@__DIR__, "..", "data")

ff_A() = readdlm(joinpath(__DATA_DIR, "ff_A.csv"), ',')
ff_B() = readdlm(joinpath(__DATA_DIR, "ff_B.csv"), ',')
ff_C() = readdlm(joinpath(__DATA_DIR, "ff_C.csv"), ',')
ff_D() = readdlm(joinpath(__DATA_DIR, "ff_D.csv"), ',')
ff_nx() = parse(Int, strip(read(joinpath(__DATA_DIR, "ff_nx.txt"), String)))

lqg_A() = readdlm(joinpath(__DATA_DIR, "lqg_A.csv"), ',')
lqg_B() = readdlm(joinpath(__DATA_DIR, "lqg_B.csv"), ',')
lqg_C() = readdlm(joinpath(__DATA_DIR, "lqg_C.csv"), ',')
lqg_D() = readdlm(joinpath(__DATA_DIR, "lqg_D.csv"), ',')
lqg_nx() = parse(Int, strip(read(joinpath(__DATA_DIR, "lqg_nx.txt"), String)))
