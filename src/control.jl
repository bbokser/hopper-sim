function u_f(u)
    # convert torque input to wrench
    # u[1] -> joint0
    # u[2] -> joint1
    # u[3] -> joint2
    # u[4] -> joint3
    # u[5] -> joint4 (parallel constraint)

    #Corresponding wrench "F" for each link
    Fk = [zeros(4); -u[1]-u[3]; 0; # body
        zeros(4); u[1]-u[2]; 0; # link0
        zeros(4); u[2]+u[5]; 0; # link1
        zeros(4); u[3]-u[4]; 0; # link2
        zeros(4); u[4]-u[5]; 0] # link3

    return Fk
end

function a_control(a_target, a_pos, a_vel)
    # a_target: joint angle target
    # a_pos: joint angles
    kp = 1
    kd = 0.05

    u = kp*(a_target-a_pos) - kd*(-a_vel)
    Fk = u_f(u)  # convert torque input to wrench

    return Fk
end
