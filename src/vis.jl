using MeshCat
using GeometryBasics
using CoordinateTransformations
using Rotations
import Plots as pl

vis = Visualizer()
render(vis)

# Zoom in--the robot falls through floor since there is no contact

delete!(vis)

function hopper_vis(qhist)
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
        sleep(0.1)
    end

    return nothing
end
