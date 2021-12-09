import Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()

using RigidBodyDynamics
using LinearAlgebra
using MeshCatMechanisms
using MeshCat
using StaticArrays
using SparseArrays
using ForwardDiff
const FD = ForwardDiff

include("utils.jl")
include("constraint.jl")
include("opt-setup.jl")
include("leg.jl")
include("control.jl")
include("vis.jl")

#Initial Conditions
rb = [0.0; 0.0; 0.5]  # start 0.5 meters above the ground
Qb = [1.0; 0; 0; 0]
q0 = 30*pi/180;
q1 = 120*(pi/180)
q2 = 150*(pi/180)
q3 = -120*(pi/180)

Tf = 1
h = 0.01
thist = 0:h:Tf
N = length(thist)

pb = rb + rotate(Qb, l_cb)  # position vector from world frame to *JOINTS* 0 and 2

Qb0 = Expq([0, q0, 0])  # quaternion from base to link 0
Q0 = L(Qb)*Qb0  # quaternion from world frame to link 0
r0 = pb + rotate(Q0, l_c0)

Q01 = Expq([0, q1, 0])  # quaternion from link 0 to link 1
Q1 = L(Q0)*Q01  # quaternion from world frame to link 1
r1 = pb + rotate(Q0, l0) + rotate(Q1, l_c1)  

Qb2 = Expq([0, q2, 0])  # quaternion from base to link 2
Q2 = L(Qb)*Qb2  # quaternion from world frame to link 2
r2 = pb + rotate(Q2, l_c2)

Q23 = Expq([0, q3, 0])  # quaternion from base to link 2
Q3 = L(Q2)*Q23  # quaternion from world frame to link 2
r3 = pb + rotate(Q2, l2) + rotate(Q3, l_c3)  

Qb .= Qb/norm(Qb)
Q0 .= Q0/norm(Q0)
Q1 .= Q1/norm(Q1)
Q2 .= Q2/norm(Q2)
Q3 .= Q3/norm(Q3)

q_0 = [rb; Qb; r0; Q0; r1; Q1; r2; Q2; r3; Q3]  # initial state

# make gravity zero for first two timesteps
ghist = repeat([0 0 g], N)'
ghist[:, 1] = [0 0 0]
ghist[:, 2] = [0 0 0]

const n_c = size(con(q_0))[1]  # number of constraint rows, 24
const n_q = size(q_0)[1]  # number of q rows, 35

#Solve with IPOPT
n_nlp = n_q + n_c + n_c  # size of decision variables, 83 
m_nlp = n_q + n_c + n_c  # size of constraint! output, 83 
nlp_prob = ProblemMOI(n_nlp,m_nlp, idx_ineq=nonnegative_constraint_indices);

#Initial conditions
qhist = zeros(n_q,N)
qhist[:,1] .= q_0
qhist[:,2] .= q_0  # this may need to be fixed

λhist = zeros(n_c,N-1)
shist = zeros(n_c,N-1)

global F = u_f([0 0 0 0 0])
global F_prev = u_f([0 0 0 0 0])

global k = 0

for kk = 2:(N-1)
    print("Simulation ", kk/(N-1)*100, " % complete \n")
    flush(stdout)
    
    if kk/(N-1)*100 > 41
        break
    end
    
    global k = kk
    
    if k == 1 || k == 2
        global F = u_f([0 0 0 0 0])  # enforce no input for first two timesteps
    else
        global F = u_f([0 0 0 0 0]*1e-4)
    end

    z_guess = [qhist[:,k]; zeros(n_q); ones(n_q)]
    z_sol = ipopt_solve(z_guess, nlp_prob, print=0);
    qhist[:,k+1] .= z_sol[1:n_q]
    λhist[:,k] .= z_sol[n_q + 1:n_q + n_c]
    shist[:,k] .= z_sol[n_q + 1 + n_c:n_q + n_c + n_c]
    
    global F_prev = F
    e = constraint_check(z_sol, 1e-6)
    # print("\n", e, "\n")
    #=
    if e == true
        print("\n Sim stopped due to ipopt infeasibility \n")
        break
    end
    =#
    
    
end

hopper_vis(qhist)
