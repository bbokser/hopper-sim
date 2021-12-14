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
include("properties.jl")
include("control.jl")
include("vis.jl")

#Initial Conditions
rb = [0.0; 0.0; 0.5]  # start 0.5 meters above the ground
Qb = [1.0; 0; 0; 0]
q0 = 30*pi/180;
q1 = 120*(pi/180)
q2 = 150*(pi/180)
q3 = -120*(pi/180)
a_target = [q0; q1; q2; q3]

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

global const n_c = size(con(q_0))[1]  # number of constraint function rows corresponding to lagrange multipliers
global const n_q = size(q_0)[1]  # number of state rows, 35

#Solve with IPOPT
n_n = 1 # number of normal forces, MUST UPDATE MANUALLY
n_b = 2 # number of friction forces (x and y), MUST UPDATE MANUALLY
n_λf = 1  # size of λf
n_s = 2  # number of complementary slackness constraints, MUST UPDATE MANUALLY
n_nlp = n_q + n_c + n_s + n_b + n_λf + n_n  # size of decision variables
n_c_ineq = 0  # no. of ineq constraints corresponding to λ (NOT λf), MUST UPDATE MANUALLY
n_c_eq = n_c-n_c_ineq
n_ineq = 4  # total number of inequality constraints, MUST UPDATE MANUALLY
n_eq = 30+5+n_c+n_b  # number of equality constraints
m_nlp = n_eq + n_ineq  # size of constraint! output

#Specify the indicies of c (constraint output) that should be non-negative.
#The rest will be treated as equality constraints.
#This can vary depending on how you stacked up c.
nonnegative_constraint_indices = (n_eq+1:m_nlp)
nlp_prob = ProblemMOI(n_nlp,m_nlp, idx_ineq=nonnegative_constraint_indices);

#Initial conditions
qhist = zeros(n_q,N)
qhist[:,1] .= q_0
qhist[:,2] .= q_0  # this may need to be fixed

λhist = zeros(n_c,N-1)
shist = zeros(n_s,N-1)

bhist = zeros(n_b,N-1)
λfhist = zeros(n_b,N-1)
nhist = zeros(n_n, N-1)

global F = u_f([0 0 0 0 0])
global F_prev = u_f([0 0 0 0 0])

global k = 0

for kk = 2:(N-1)
    
    global k = kk

    a = a_joint(qhist[:, k])
    a_prev = a_joint(qhist[:, k-1])

    if k == 1 || k == 2  # enforce no input for first two timesteps
        global F = u_f([0 0 0 0 0])  
    else
        global F = u_f([0 0 0 0 0])*1e-4
        # global F = a_control(a_target, a, a_vel(a, a_prev, h))
    end

    z_guess = [qhist[:,k]; zeros(n_c); ones(n_s); zeros(n_b); ones(n_b); zeros(n_n)]
    z_sol = ipopt_solve(z_guess, nlp_prob, max_iter=1000, print=0);
    qhist[:,k+1] .= z_sol[1:n_q]
    λhist[:,k] .= z_sol[n_q + 1:n_q + n_c]
    shist[:,k] .= z_sol[n_q + n_c + 1:n_q + n_c + n_s]
    bhist[:,k] .= z_sol[n_q + n_c + n_s + 1:n_q + n_c + n_s + n_b]
    λfhist[:,k] .= z_sol[n_q + n_c + n_s + n_b + 1:n_q + n_c + n_s + n_b + n_λf]
    nhist[:,k] .= z_sol[n_q + n_c + n_s + n_b + n_λf + 1:end]

    global F_prev = F

    e = constraint_check(z_sol, 1e-6)  # print("\n", e, "\n")
    # f e == true; break; end
    print("Simulation ", round(kk/(N-1)*100, digits=3), " % complete \n")
    # flush(stdout)
    
    if kk/(N-1)*100 > 32; break; end
    
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

function angle_look()
    an = zeros(N)
    for i in 1:(N-1)
        Q0 = qhist[11:14, i]
        Q1 = qhist[18:21, i]
        an[i] = pi - anglesolve(L(Q0)'*Q1)
    end
    return an
end

function angle_y_look()
    an = zeros(N)
    for i in 1:(N-1)
        Q0 = qhist[11:14, i]
        Q1 = qhist[18:21, i]
        an[i] = angle_y(Q0, Q1)
    end
    return an
end

plot = false

if plot == true
    ph = pl.plot(thist,signed_d(), title="signed dist from foot to ground plane")
    pbz = pl.plot(thist,qhist[3,:], title="height of body")
    plam = pl.plot(λhist[n_c,:],title="contact force")
    # pslack = pl.plot(shist,title="slackvar")
    panbt = pl.plot(thist, angle_y_look().*180/pi,title="angle_y b/t 0 and 1")
    pb = pl.plot(thist[:, 1:end-1], bhist',title="angle_y b/t 0 and 1")
    pl.display(ph)
    pl.display(pbz)
    pl.display(plam)
    # pl.display(pslack)
    pl.display(panbt)
    pl.display(pb)
end

urdf = true  # for now manually change this

if urdf == false
    anim_init = geom_init
    anim = geom_vis
else
    anim_init = urdf_init
    anim = urdf_vis
end

mvis = anim_init()
print("\n Visualization starting in 5 seconds \n")
sleep(5)
for j in 1:5
    print("\n Visualization starting now, replay #", j, "\n")
    anim(mvis, qhist)
end
print("\n Visualization ended \n")