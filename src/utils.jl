function hat(ω)
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
    q = zeros(4)
    θ = norm(ϕ)
    q = [cos(θ/2); 0.5*ϕ*sinc(θ/(2*pi))]
    
    return q
end

function rotate(Q, p)
    # Rotate a position vector p by a quaternion Q
    return H'L(Q)*R(Q)'*H*p
end