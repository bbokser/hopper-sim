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
qe_target = [q0; q1; q2; q3]

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

global const n_c = size(con(q_0))[1]  # number of constraint rows, 24
global const n_q = size(q_0)[1]  # number of q rows, 35

#Specify the indicies of c (constraint output) that should be non-negative.
#The rest will be treated as equality constraints.
#This can vary depending on how you stacked up c.
n_eq = 30+5+n_c-3  # number of equality constraints, 58
nonnegative_constraint_indices = (n_eq+1:n_eq+4) # update this

#Solve with IPOPT
n_nlp = n_q + n_c + 1  # size of decision variables
m_nlp = n_q + n_c + 1  # size of constraint! output
nlp_prob = ProblemMOI(n_nlp,m_nlp, idx_ineq=nonnegative_constraint_indices);

#Initial conditions
qhist = zeros(n_q,N)
qhist[:,1] .= q_0
qhist[:,2] .= q_0  # this may need to be fixed

λhist = zeros(n_c,N-1)
shist = zeros(1,N-1)

global F = u_f([0 0 0 0 0])
global F_prev = u_f([0 0 0 0 0])

global k = 0

for kk = 2:(N-1)
    
    global k = kk
    
    if k == 1 || k == 2
        global F = u_f([0 0 0 0 0])  # enforce no input for first two timesteps
    else
        global F = u_f([0 0 0 0 0]*1e-5)
        # global F = qe_control(qe_target, qe_pos, b_orient)
    end

    z_guess = [qhist[:,k]; zeros(n_c); 1]
    z_sol = ipopt_solve(z_guess, nlp_prob, print=0);
    qhist[:,k+1] .= z_sol[1:n_q]
    λhist[:,k] .= z_sol[n_q + 1:n_q + n_c]
    shist[:,k] .= z_sol[n_q + n_c + 1]  # only one slack variable
    
    global F_prev = F
    e = constraint_check(z_sol, 1e-6)
    print("Simulation ", round(kk/(N-1)*100, digits=3), " % complete \n")
    flush(stdout)
    # print("\n", e, "\n")
    
    if e == true
        print("\n Sim stopped due to ipopt infeasibility \n")
        break
    end
    
    #=
    if kk/(N-1)*100 > 41
        break
    end
    =#
    
end

function signed_d()
    hf = zeros(N)
    for i in 1:(N-1)
        r3 = qhist[29:31, i]
        Q3 = qhist[32:35, i]
        hf[i] = r3[3] + rotate(Q3, lee-l_c3)[3] # - 0.025 # position of foot
    end
    return hf
end

ph = pl.plot(thist,signed_d(), title="signed dist from foot to ground plane")
pbz = pl.plot(thist,qhist[3,:], title="height of body")
plam = pl.plot(λhist[24,:],title="contact force")
pslack = pl.plot(shist[1, :],title="slackvar")

pl.display(ph)
pl.display(pbz)
pl.display(plam)
pl.display(pslack)

for j in 1:5
    hopper_vis(qhist)
end