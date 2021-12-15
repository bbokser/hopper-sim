function con(q)
    # joint constraint function
    # c(q_0) should be all zeros
    c_ = zeros(eltype(q), 24)
    
    rb = q[1:3]
    Qb = q[4:7]
    r0 = q[8:10]
    Q0 = q[11:14]
    r1 = q[15:17]
    Q1 = q[18:21]
    r2 = q[22:24]
    Q2 = q[25:28]
    r3 = q[29:31]
    Q3 = q[32:35]
    
    rf = r3 + rotate(Q3, lee-l_c3) # position of foot

    pb = rb + rotate(Qb, l_cb) # position vector from world frame to *JOINTS* 0 and 2
    c_[1:3] = pb - r0 - rotate(Q0, -l_c0)
    c_[4:5] = [0 1 0 0; 0 0 0 1]*L(Qb)'*Q0  # y axis rotation constraint: set x & z = 0
    c_[6:8] = r0 + rotate(Q0,  l0-l_c0) - r1 - rotate(Q1, -l_c1)
    c_[9:10] = [0 1 0 0; 0 0 0 1]*L(Q0)'*Q1
    c_[11:13] = pb - r2 - rotate(Q2, -l_c2)
    c_[14:15] = [0 1 0 0; 0 0 0 1]*L(Qb)'*Q2
    c_[16:18] = r2 + rotate(Q2, l2 - l_c2) - r3 - rotate(Q3, -l_c3)
    c_[19:20] = [0 1 0 0; 0 0 0 1]*L(Q2)'*Q3 
    c_[21:23] = r1 + rotate(Q1, l1 - l_c1) - r3 - rotate(Q3, lc-l_c3)
    # Prevent base from rotating except in z axis
    # c_[24:25] = [0 1 0 0; 0 0 1 0]*Qb
    #c_[24] = rb[1]  # constrain base x axis
    #c_[25] = rb[2]  # constrain base y axis
    # Prevent foot-floor interpenetration
    c_[24] = rf[3] - 0.025  # subtract radius of foot
    return c_
end

function Dc(q)
    ForwardDiff.jacobian(dq->con(dq),q)*Ḡ(q)
end

function del(m, I, r1, r2, r3, Q1, Q2, Q3, grav)
    [(1/h)*m*(r2-r1) - (1/h)*m*(r3-r2) - m*grav*h;
     (2.0/h)*G(Q2)'*L(Q1)*H*I*H'*L(Q1)'*Q2 + (2.0/h)*G(Q2)'*T*R(Q3)'*H*I*H'*L(Q2)'*Q3]
end
    
function DEL(q_1,q_2,q_3,λ,F1,F2,grav)
    
    rb_1 = q_1[1:3]
    Qb_1 = q_1[4:7]
    r0_1 = q_1[8:10]
    Q0_1 = q_1[11:14]
    r1_1 = q_1[15:17]
    Q1_1 = q_1[18:21]
    r2_1 = q_1[22:24]
    Q2_1 = q_1[25:28]
    r3_1 = q_1[29:31]
    Q3_1 = q_1[32:35]
    
    rb_2 = q_2[1:3]
    Qb_2 = q_2[4:7]
    r0_2 = q_2[8:10]
    Q0_2 = q_2[11:14]
    r1_2 = q_2[15:17]
    Q1_2 = q_2[18:21]
    r2_2 = q_2[22:24]
    Q2_2 = q_2[25:28]
    r3_2 = q_2[29:31]
    Q3_2 = q_2[32:35]
    
    rb_3 = q_3[1:3]
    Qb_3 = q_3[4:7]
    r0_3 = q_3[8:10]
    Q0_3 = q_3[11:14]
    r1_3 = q_3[15:17]
    Q1_3 = q_3[18:21]
    r2_3 = q_3[22:24]
    Q2_3 = q_3[25:28]
    r3_3 = q_3[29:31]
    Q3_3 = q_3[32:35]

    del1 = [del(mb, Ib, rb_1, rb_2, rb_3, Qb_1, Qb_2, Qb_3, grav);
            del(m0, I0, r0_1, r0_2, r0_3, Q0_1, Q0_2, Q0_3, grav);
            del(m1, I1, r1_1, r1_2, r1_3, Q1_1, Q1_2, Q1_3, grav);
            del(m2, I2, r2_1, r2_2, r2_3, Q2_1, Q2_2, Q2_3, grav);
            del(m3, I3, r3_1, r3_2, r3_3, Q3_1, Q3_2, Q3_3, grav)] 
    #if F1[5] != F2[5]
    #        print(F1[5], ", ", F2[5], "\n")
    #end
    return del1 + (h/2.0)*F1 + (h/2.0)*F2 + reshape(h*Dc(q_2)'*λ, 30)
end

#Objective and constraint functions for IPOPT

function objective(z)
    s = z[n_q+n_c+1:n_q+n_c+n_s]
    λj = z[n_q+1:n_q+n_c-1]
    α = 1e-3
    # regularizer (Tikhonov)
    return α*λj'*λj + sum(s) #Minimize slacks associated with complementarity conditions
    # return sum(s)
end

function constraint!(c,z)
    qn = z[1:n_q]
    λ = z[n_q+1:n_q+n_c]
    s = z[n_q+n_c+1]

    # nonlinear DEL equation                (30 equality constraint)   c1
    # quaternion norm squared - 1 = 0       (5 equality constraints)   c2-c6
    # joint constraints                     (23 equality constraints)  c7
    
    # signed distance of foot from ground   (1 inequality constraint)  c8
    # collision complementarity foot        (1 inequality constraint)  c9

    ϕ = con(qn)[n_c]
    n = λ[n_c]  # normal force 1x1 (z axis only, assumes ground is flat)

    # equality constraints
    c1 = DEL(qhist[:,k-1], qhist[:,k], qn, λ, F_prev, F, ghist[:, k])  # 30x1
    c2 = norm(qn[4:7])^2 - 1  # 1x1
    c3 = norm(qn[11:14])^2 - 1
    c4 = norm(qn[18:21])^2 - 1
    c5 = norm(qn[25:28])^2 - 1
    c6 = norm(qn[32:35])^2 - 1
    c7 = con(qn)[1:n_c-1]
    # @show size([c1; c2; c3; c4; c5; c6; c7])
    # inequality constraints
    c8 = ϕ  # signed distance
    c9 = s - n*ϕ  # relaxed complementarity (signed dist) 1x1
    c .= [c1; c2; c3; c4; c5; c6; c7; c8; c9]

    return nothing
end

function primal_bounds(n)
    #Enforce simple bound constraints on the decision variables (e.g. positivity) here
    # x_l ≤ [q; λ; s] ≤ x_u
    x_l = -Inf*ones(n)  # 60
    
    # x_l[n_q+n_c] = 0  # normal force corresponding to signed dist constr
    x_l[n_q+n_c+1] = 0  # slack variable
    
    x_u = Inf*ones(n)
    
    #x_l[1:2] = zeros(2)  # constrain body movement to z-axis only
    #x_u[1:2] = zeros(2)  # constrain body movement to z-axis only
    return x_l, x_u
end

function constraint_check(z, n_tol)
    qn = z[1:n_q]
    λ = z[n_q+1:n_q+n_c]
    s = z[n_q+n_c+1]

    ϕ = con(qn)[n_c]
    n = λ[n_c]  # normal force 1x1 (z axis only, assumes ground is flat)

    # equality constraints
    c1 = DEL(qhist[:,k-1], qhist[:,k], qn, λ, F_prev, F, ghist[:, k])  # 30x1
    c2 = norm(qn[4:7])^2 - 1  # 1x1
    c3 = norm(qn[11:14])^2 - 1
    c4 = norm(qn[18:21])^2 - 1
    c5 = norm(qn[25:28])^2 - 1
    c6 = norm(qn[32:35])^2 - 1
    c7 = con(qn)[1:n_c-1]
    
    # inequality constraints
    c8 = ϕ  # signed distance
    c9 = s - n*ϕ  # relaxed complementarity (signed dist) 1x1
    
    A = [c1; c2; c3; c4; c5; c6; c7]
    B = [c8; c9]

    if !isapprox(A, zeros(size(A)[1]); atol=n_tol) # , rtol=0)  # 58
        e = 1
        #print("\n", A, "\n")
        print("\n A \n")
        print(findall(A .< -ones(size(A)[1])*n_tol), " is less than 0 \n")
        print(findall(A .> ones(size(A)[1])*n_tol), " is greater than 0 \n")
        print("\n Sim stopped due to ipopt infeasibility \n")
    elseif B < -ones(size(B)[1])*n_tol  #25
        e = 1
        #print("\n", B, "\n")
        print("\n B \n")
        print(findall(B .< -ones(size(B)[1])*n_tol))  # 25
        print("\n Sim stopped due to ipopt infeasibility \n")
    elseif (pi-angle_y(Q0, Q1)) < 18*pi/180  # constrain relative angle between links 0 and 1
        e = 1
        print("RoM lower limit surpassed")
    elseif 166*pi/180 < (pi-angle_y(Q0, Q1))
        e = 1
        print("RoM upper limit surpassed")
    else
        e = 0
    end
    
    return e
end