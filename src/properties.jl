# link lengths
const l0 = [0.1; 0; 0]
const l1 = [0.3; 0; 0]
const l2 = [0.3; 0; 0]
const l3 = 0.1
const l4 = 0.2
const l5 = 0.0205
const lc = [l3; 0; 0]
const lee = [l3 + l4; 0; l5] # sqrt((l3 + l4)^2 + l5^2)

# CoM locations
const l_cb = [0; 0.004; 0]
const l_c0 = [0.0125108364230515; 0.00117191218927888; 0]
const l_c1 = [0.149359714867044; 0; 0]
const l_c2 = [0.0469412900551914; 0; 0]
const l_c3 = [0.113177000131857; 0; -0.015332867880069]

# link masses
const mb = 7  # kg
const m0 = 0.24644240965487
const m1 = 0.0707939028219395
const m2 = 0.276735496985514
const m3 = 0.130824780046739
# const m = Diagonal([m0, m1, m2, m3])
    
# gravity, obviously
const g = 9.807

# mass moment of inertia in axis of rotation
const Ib = Array([0.0024241 5.252E-06 2.0733E-19; 
                  5.252E-06 0.0044176 -3.1153E-19; 
                  2.0733E-19 -3.1153E-19 0.0022481])

const I0 = Array([3.83120149546952E-05 1.46925714738609E-05 -8.60106401672571E-06;
                  1.46925714738609E-05 0.000172067745507247 1.0427260925207E-06;
                  -8.60106401672571E-06 1.0427260925207E-06 0.00014745218068435])

const I1 = Array([3.06999775886187E-06 7.91090301514898E-12 -1.43705963146176E-12;
                  7.91090301514898E-12 0.000147960574744097 1.30742394049546E-11;
                  -1.43705963146176E-12 1.30742394049546E-11 0.000147884231885009])

const I2 = Array([3.43038397803592E-05 -2.90339844227483E-07 6.18680397558952E-06;
                  -2.90339844227483E-07 0.000302324068012293 2.25016327583562E-08;
                  6.18680397558952E-06 2.25016327583562E-08 0.00028292376778719])

const I3 = Array([1.76996970020568E-05 -5.3695427116208E-07 7.62350214406387E-07;
                  -5.3695427116208E-07 0.000164188445564489 -2.77843753828047E-07;
                  7.62350214406387E-07 -2.77843753828047E-07 0.000160656046697151])

M̄ = [mb*I(3) zeros(3, 27);
    zeros(3,3) Ib zeros(3, 24);
    zeros(3,6) m0*I(3) zeros(3,21);
    zeros(3,9) I0 zeros(3,18);
    zeros(3,12) m1*I(3) zeros(3,15);
    zeros(3,15) I1 zeros(3, 12);
    zeros(3,18) m2*I(3) zeros(3, 9);
    zeros(3, 21) I2 zeros(3, 6);
    zeros(3, 24) m3*I(3) zeros(3, 3);
    zeros(3, 27) I3]

# Kinematics
function kin_ee(q)
    # forward kinematics of the end effector
    rb = q[1:3]
    Qb = q[4:7]
    r3 = q[29:31]
    Q2 = q[25:28]
    Q3 = q[32:35]
    
    # not the most efficient way, but will have to do for now
    pb = rb + rotate(Qb, l_cb)  # position vector from world frame to *JOINTS* 0 and 2
    
    ree = pb + rotate(Q2, l2) + rotate(Q3, lee)  
    
    # ree = r3 + rotate(Q3, lee-l3) # TODO: this is probably better, try

    return ree  # 3x1
end

function J(q)
    # ree = kin_ee(q_min);
    jac = ForwardDiff.jacobian(dq->kin_ee(dq), q)  # 3x35
    return jac
end