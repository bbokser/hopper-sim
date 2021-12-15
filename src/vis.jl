function geom_init()
    vis = Visualizer()
    render(vis)
    delete!(vis)
    return vis
end

function geom_vis(vis, qhist, dt, speed=1)

    green_material = MeshPhongMaterial(color=pl.RGBA(0, 1, 0, 0.8))
    red_material = MeshPhongMaterial(color=pl.RGBA(1, 0, 0, 0.8))

    cylinderb = Cylinder(Point(0, -0.02, l_cb[2]), Point(0, 0.02, l_cb[2]), 0.04)  #float(l_cb)
    cylinder0 = Cylinder(Point(-l_c0[1], -l_c0[2], -l_c0[3]), Point((l0-l_c0)[1], (l0-l_c0)[2], (l0-l_c0)[3]), 0.008)
    cylinder1 = Cylinder(Point(-l_c1[1], -l_c1[2], -l_c1[3]), Point((l1-l_c1)[1],(l1-l_c1)[2],(l1-l_c1)[2]), 0.008)
    cylinder2 = Cylinder(Point(-l_c2[1], -l_c2[2], -l_c2[3]), Point((l2-l_c2)[1], (l2-l_c2)[2], (l2-l_c2)[3]), 0.008)
    cylinder3 = Cylinder(Point(-l_c3[1], -l_c3[2], -l_c3[3]), Point((lee-l_c3)[1], (lee-l_c3)[2],  (lee-l_c3)[3]), 0.008)

    setobject!(vis["cylinderb"],cylinderb,red_material)
    setobject!(vis["cylinder0"],cylinder0,green_material)
    setobject!(vis["cylinder1"],cylinder1,green_material)
    setobject!(vis["cylinder2"],cylinder2,green_material)
    setobject!(vis["cylinder3"],cylinder3,green_material)

    for k = 1:N
        
        q = qhist[:, k]
        
        # set position and attitude
        positionb = Translation(q[1:3]...)
        attitudeb = LinearMap(UnitQuaternion(q[4:7]))
        position0 = Translation(q[8:10]...)
        attitude0 = LinearMap(UnitQuaternion(q[11:14]))
        position1 = Translation(q[15:17]...)
        attitude1 = LinearMap(UnitQuaternion(q[18:21]))
        position2 = Translation(q[22:24]...)
        attitude2 = LinearMap(UnitQuaternion(q[25:28]))
        position3 = Translation(q[29:31]...)
        attitude3 = LinearMap(UnitQuaternion(q[32:35]))

        settransform!(vis["cylinderb"], compose(positionb,attitudeb))
        settransform!(vis["cylinder0"], compose(position0,attitude0))
        settransform!(vis["cylinder1"], compose(position1,attitude1))
        settransform!(vis["cylinder2"], compose(position2,attitude2))
        settransform!(vis["cylinder3"], compose(position3,attitude3))
        sleep(dt/speed)
    end

    return nothing
end

function urdf_init()
    curdir = pwd()
    urdfpath = joinpath(curdir, "res/flyhopper_robot/urdf/flyhopper_robot_base_call.urdf")
    leg = parse_urdf(urdfpath, floating=true)
    # state = MechanismState(leg)
    # mvis = MechanismVisualizer(doublependulum, Skeleton(randomize_colors=true, inertias=false));
    mvis = MechanismVisualizer(leg, URDFVisuals(urdfpath));

    render(mvis)
    return mvis
end

function urdf_vis(mvis, qhist, dt, speed=1)

    for k = 1:N
        q = copy(qhist[:, k])
        pb = q[1:3]
        Qb = q[4:7]
        Q0 = q[11:14]
        Q1 = q[18:21]
        Q2 = q[25:28]
        Q3 = q[32:35]
           
        # convert quaternions to joint angles
        a0 = -angle_y(Qb, Q0) +30*(pi/180)
        a1 = -angle_y(Q0, Q1) +120*(pi/180)
        a2 = -angle_y(Qb, Q2) +150*(pi/180)
        a3 = -angle_y(Q2, Q3) -120*(pi/180)
    
        #-0.5235987755982988, -2.6179938779914944, -2.0943951023931953, 2.0943951023931953
        q_array = vcat(Qb, pb, [a0, a2, a1, a3])
        # @show [a0, a2, a1, a3]
        set_configuration!(mvis, q_array)
        sleep(dt/speed)
    end
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