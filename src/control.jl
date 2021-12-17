# Control functions

function u_f(u, q)
    # convert torque input to wrench
    Qb = q[4:7]
    Q0 = q[11:14]
    Q2 = q[25:28]

    u0 = u[1]  # joint0
    u2 = u[2]  # joint2

    τb = [0; -u0-u2; 0]
    τ0 = [0; u0;     0]
    τ2 = [0; u2;     0]

    fb = f_applied(τb, l_cb)  # force in body frame
    f0 = f_applied(τ0, l_c0)
    f2 = f_applied(τ2, l_c2)

    #Corresponding wrench "F" for each link
    Fk = [rotate(Qb, fb); τb;        # body  
          rotate(Q0, f0); τ0;        # link0
          zeros(3);       zeros(3);  # link1
          rotate(Q2, f2); τ2;        # link2
          zeros(3);       zeros(3)]  # link3

    return Fk
end

function f_applied(τ, r)
    cross(τ, r)/(norm(r)^2)
end

function a_control(a_target, a_pos, a_vel, q)
    # a_target: joint angle target
    # a_pos: joint angles
    kp = 5
    kd = copy(kp)*0.03
    # print("a_pos = ", a_pos.*180/pi, "\n")
    u = kp*(a_pos-a_target) + kd*(a_vel)
    # B = Diagonal([1; 0; 1; 0])  # actuator selection matrix
    Fk = u_f(u, q)  # convert torque input to wrench
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