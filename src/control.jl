function u_f(u)
    # convert torque input to wrench
    # u[1] -> joint0
    # u[2] -> joint1
    # u[3] -> joint2
    # u[4] -> joint3
    # u[5] -> joint4 (parallel constraint)

    u5 = 0  # we're not outputting torque there anyway

    #Corresponding wrench "F" for each link
    Fk = [zeros(4); -u[1]-u[3]; 0; # body
          zeros(4); u[1]-u[2]; 0; # link0
          zeros(4); u[2]+u5; 0; # link1
          zeros(4); u[3]-u[4]; 0; # link2
          zeros(4); u[4]-u5; 0] # link3

    return Fk
end

function a_control(a_target, a_pos, a_vel)
    # a_target: joint angle target
    # a_pos: joint angles
    kp = 0.000007
    kd = copy(kp)*0.05
    print("a_pos = ", a_pos.*180/pi, "\n")
    
    B = Diagonal([1; 0; 1; 0])  # actuator selection matrix
    u = kp*(a_pos-a_target) + kd*(a_vel)
    Fk = u_f(B*u)  # convert torque input to wrench
    return Fk
end

function max_to_min(q)
    # convert max coords state to min coords
    rb = q[1:3]
    Qb = q[4:7]
    a = a_joint(q)
    return [rb; Qb; a]
end

# Kinematics
function kin_ee_min(q_min)
    # forward kinematics of the end effector from body pose and joint positions
    rb = q_min[1:3]
    Qb = q_min[4:7]
    a = q_min[8:11]
    q2 = a[3]
    q3 = a[4]
    
    pb = rb + rotate(Qb, l_cb)  # position vector from world frame to *JOINTS* 0 and 2
    
    Qb2 = Expq([0, q2, 0])  # quaternion from base to link 2
    Q2 = L(Qb)*Qb2  # quaternion from world frame to link 2
    
    Q23 = Expq([0, q3, 0])  # quaternion from link 2 to link 3
    Q3 = L(Q2)*Q23  # quaternion from world frame to link 2
    ree = pb + rotate(Q2, l2) + rotate(Q3, lee)  

    return ree
end

# End effector Jacobian
function J_min(q_min)
    # ree = kin_ee(q_min);
    jac = ForwardDiff.jacobian(dq->kin_ee_min(dq), q_min)
    return jac
end