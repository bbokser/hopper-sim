#Torque input at joints

#=
uhist = repeat([0 0 0 0 0]*1e-4, N)'
uhist[:, 1] = [0 0 0 0 0]
uhist[:, 2] = [0 0 0 0 0]
# 1 -> joint0
# 2 -> joint1
# 3 -> joint2
# 4 -> joint3
# 5 -> joint4 (parallel constraint)
#Corresponding F
Fhist = zeros(30,N)  # wrench [xyz force, xyz torque]
for k = 1:N
    Fhist[:,k] = [zeros(4); -uhist[1, k]-uhist[3, k]; 0; # body
                zeros(4); uhist[1, k]-uhist[2, k]; 0; # link0
                zeros(4); uhist[2, k]+uhist[5, k]; 0; # link1
                zeros(4); uhist[3, k]-uhist[4, k]; 0; # link2
                zeros(4); uhist[4, k]-uhist[5, k]; 0] # link3
end             
=#

function u_f(u)
    # convert torque input to wrench

    # 1 -> joint0
    # 2 -> joint1
    # 3 -> joint2
    # 4 -> joint3
    # 5 -> joint4 (parallel constraint)
    #Corresponding F
    Fk = [zeros(4); -u[1]-u[3]; 0; # body
        zeros(4); u[1]-u[2]; 0; # link0
        zeros(4); u[2]+u[5]; 0; # link1
        zeros(4); u[3]-u[4]; 0; # link2
        zeros(4); u[4]-u[5]; 0] # link3

    return Fk
end


function qcontrol(q_target, q_pos,  b_orient)
    
    Fk = u_f(u)

    return Fk
end
