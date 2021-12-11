function hat(ω)
    # skew-symmetric
    return [0 -ω[3] ω[2];
            ω[3] 0 -ω[1];
            -ω[2] ω[1] 0]
end

function L(Q)
    [Q[1] -Q[2:4]'; Q[2:4] Q[1]*I + hat(Q[2:4])]
end

function R(Q)
    [Q[1] -Q[2:4]'; Q[2:4] Q[1]*I - hat(Q[2:4])]
end

H = [zeros(1,3); I];

T = Diagonal([1.0; -1.0; -1.0; -1.0])

function G(Q)
    return L(Q)*H  # 4x3
end

function Ḡ(q)
    Qb = q[4:7]
    Q0 = q[11:14]
    Q1 = q[18:21]
    Q2 = q[25:28]
    Q3 = q[32:35]
    
    # 35x30
    return blockdiag(sparse(I, 3, 3), sparse(G(Qb)), 
                     sparse(I, 3, 3), sparse(G(Q0)), 
                     sparse(I, 3, 3), sparse(G(Q1)),
                     sparse(I, 3, 3), sparse(G(Q2)),
                     sparse(I, 3, 3), sparse(G(Q3)),)
end

function Expq(ϕ)
    
    # The quaternion exponential map ϕ → q 
    Q = zeros(4)
    θ = norm(ϕ)
    Q = [cos(θ/2); 0.5*ϕ*sinc(θ/(2*pi))]
    
    return Q
end

function rotate(Q, p)
    # Rotate a position vector p by a quaternion Q
    return H'L(Q)*R(Q)'*H*p
end

function anglesolve(Q)
    # convert quaternion to angle
    2 * atan(norm(Q[2:4]),Q[1])
end

function anglebt(Q1, Q2)
    # get smallest angle b/t quaternions
    acos((Q1'*Q2) / (norm(Q1)*norm(Q2)))
end

function a_joint(q)
    # convert quaternions to relative joint angles b/t links

    Qb = q[4:7]
    Q0 = q[11:14]
    Q1 = q[18:21]
    Q2 = q[25:28]
    Q3 = q[32:35]
    
    Qb0 = L(Qb)'*Q0
    Q01 = L(Q0)'*Q1
    Qb2 = L(Qb)'*Q2
    Q23 = L(Q2)'*Q3

    a0 = -anglesolve(Qb0)
    a1 = -anglesolve(Q01)
    a2 = -anglesolve(Qb2)
    a3 = anglesolve(Q23)

    return [a0; a1; a2; a3]

end

function a_vel(a, a_prev, dt)
    # get joint velocities from current and previous joint angle
    return (a .- a_prev)/dt
end
#=
function a_urdf(a)
    # adjust angles for urdf animation
    a0 = a[1] +30*(pi/180)
    a1 = a[2] +120*(pi/180)
    a2 = a[3] +150*(pi/180)
    a3 = a[4] -120*(pi/180)
    
    return [a0, a2, a1, a3]
end
=#