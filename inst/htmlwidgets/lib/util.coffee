### UTIL ###

cart2polar = (point) ->
    return {"r":Math.sqrt(point.x**2 + point.y**2), "angle":Math.atan2(point.y, point.x)}
polar2cart = (ppoint) ->
    return {"x":Math.cos(ppoint.angle) * ppoint.r, "y":Math.sin(ppoint.angle) * ppoint.r}
deg2rad = (deg) ->
    return deg * Math.PI/180
rad2deg = (rad) ->
    return rad * 180/Math.PI

hexagon_polar = (angle, radius=1) ->
    delta = 2*Math.PI/6 # the angle between every consecutive vertex of the hexagon
    return Math.cos(delta/2)/Math.cos((math.mod(angle, delta)) - delta/2) * radius # we don't use the % operator here, because it behaves "wrongly" for negative angles
hexagon_edge = (angle) ->
    # this function returns the edge number made by and angle on the hexagon
    delta = 2*Math.PI/6
    return math.mod(Math.floor(angle/delta), 6)
hexagon_edgeclip = (angle, edgeid, padding=0) ->
    delta = 2*Math.PI/6
    minangle = edgeid * delta + padding
    maxangle = (edgeid + 1) * delta - padding

    if math.mod(angle, 2*Math.PI) > maxangle
        return maxangle
    else if math.mod(angle, 2*Math.PI) < minangle
        return minangle
    else
        return angle
hexagon_path = (angle1, angle2) ->
    angles = [angle1]
    for angle in [0..math.pi*2-0.001] by math.pi/3
        if between_circular(angle, angle1, angle2)
            angles.push(angle)
    angles.push(angle2)
    coords = (polar2cart({angle:angle, r:hexagon_polar(angle, 1)}) for angle in angles)

    return coords
hexagon_edge_alignment = (edge) ->
    if edge == 1
        va = "-100%"
    else if edge == 4
        va = "0%"
    else
        va = "-50%"

    if edge == 0 || edge == 5
        ha = "start"
    else if edge == 1 || edge == 4
        ha = "middle"
    else
        ha = "end"

    return [va, ha]

class Binner
    constructor: (@nbins) ->
        @delta = math.pi*2/@nbins
        @bins = (@delta/2 + i * @delta for i in [0..@nbins-1])

    bin: (angles) ->
        binned = (0 for i in [0..@nbins-1])
        for angle in angles
            binid = math.mod(math.floor((angle + @delta/2)/@delta), @nbins)
            binned[binid] += 1

        return binned

    group: (barycoords) ->
        grouped = ([] for i in [0..@nbins-1])
        for row in barycoords
            binid = math.mod(math.floor((row.angle + @delta/2)/@delta), @nbins)
            grouped[binid].push(row)

        return grouped

cal_circular_mean = (angles) ->
	return math.atan2(math.mean(math.sin(angles)), math.mean(math.cos(angles)))

between_circular = (angle, angle1, angle2) ->
    return math.mod((angle-angle1),math.pi*2) < (angle2 - angle1)

difference_circular = (angle1, angle2) ->
    return math.atan2(math.sin(angle2-angle1), math.cos(angle2-angle1))

bbox_overlap = (rect1, rect2, widthpadding=0, heightpadding=0) ->
    !(rect1.x > rect2.x + rect2.width + widthpadding || rect2.x > rect1.x + rect1.width + widthpadding || rect1.y > rect2.y + rect2.height + heightpadding || rect2.y > rect1.y + rect1.height + heightpadding)

# does the angle increase?, especially tricky around 0
angles_increase = (angle1, angle2) ->
    math.mod(angle2-angle1, math.pi*2) < math.pi

posneg = (value) ->
    if value == 0
        return 1
    return value/math.abs(value)



#from http://stackoverflow.com/questions/9043805/test-if-two-lines-intersect-javascript-function
lineIntersect = (x1,y1,x2,y2,x3,y3,x4,y4) ->
    x=((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
    y=((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
    if (isNaN(x)||isNaN(y))
        return false
    else
        if (x1>=x2)
            if (!(x2<=x&&x<=x1))
                return false
        else
            if (!(x1<=x&&x<=x2))
                return false

        if (y1>=y2)
            if (!(y2<=y&&y<=y1))
                return false
        else
            if (!(y1<=y&&y<=y2))
                return false

        if (x3>=x4)
            if (!(x4<=x&&x<=x3))
                return false
        else
            if (!(x3<=x&&x<=x4))
                return false

        if (y3>=y4)
            if (!(y4<=y&&y<=y3))
                return false
        else
            if (!(y3<=y&&y<=y4))
                return false
    true
