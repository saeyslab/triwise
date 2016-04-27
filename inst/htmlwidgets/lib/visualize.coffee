## VISUALIZE ##

drawHexagonGrid = (ax, rmax, scale, rdelta=1, rmin=0, anglebase=0, color="#444444", alpha=0.7, lw=1) ->
    grid = ax.append("g")

    hexagon_points = []
    for angle in [anglebase..anglebase+Math.PI*2] by Math.PI/6*2
        hexagon_points.push(polar2cart({"angle":angle, "r":1}))

    for r in [rmin..rmax] by rdelta
        grid.append("polygon")
            .attr("points", (scale(point.x) * r + "," + scale(point.y) * r for point in hexagon_points).join(" "))
            .style("fill", "none")
            .style("stroke", color)
            .style("opacity", alpha)
            .style("stroke-width", lw)
        # hexagon = Polygon(np.array(hexagon_points), fill=None, zorder=10, lw=lw, color=color, alpha=alpha)
        # ax.add_artist(hexagon)

    return grid

drawCircleGrid = (ax, scale, radii=null, rmax=1, rdelta=1, rmin=0, anglebase=0, color="#444444", alpha=0.7, lw=1) ->
    grid = ax.append("g")

    if radii == null
        radii = (r for r in [rmin..rmax] by rdelta)

    grid.selectAll("circle")
        .data(radii)
        .enter()
        .append("circle")
        .attr("r", (d) -> scale(d))
        .style("fill", "none")
        .style("stroke", color)
        .style("opacity", alpha)
        .style("stroke-width", lw)
    return grid


drawDirections = (ax, rmax, scale, labels, anglebase=0, color="black", padding=0.1, lw=1) ->
    directions = ax.append("g")
        .classed("directions", true)
    angles = (angle for angle in [anglebase..anglebase+Math.PI*2-0.001] by Math.PI/3*2)
    for [angle, label] in _.zip(angles, labels)
        point = polar2cart({"angle":angle, "r":rmax})

        directions.append("line")
            .attr({
                x1:0,
                y1:0,
                x2:scale(point.x),
                y2:scale(point.y)
                })
            .style("stroke-width", lw)
            .style("stroke", color)

        point = polar2cart({"angle":angle, "r":rmax*(padding+1)})
        if point.y < 0
            va = "0%"
        else if point.y == 0
            va = "-30%"
        else
            va = "-75%"

        if scale(-0.1)<scale(point.x)<scale(0.1)
            ha = "middle"
        else if point.x < 0
            ha = "end"
        else
            ha = "start"

        directions.append("text")
            .attr("x",  scale(point.x))
            .attr("y",  scale(point.y))
            .text(label)
            .style("baseline-shift", va)
            .style("text-anchor", ha)
            .attr("class", "direction_label")
            .style("font-weight", "bold")
    return directions

repositionDirections = (directions, bboxes) ->
    

class Roseplot
    constructor: (@ax, @w, @h, @barycoords, @labels, @binner, @colorDirection) ->
        @rmax = 1
        @scale = (x) => x * (@h-0)/(@rmax*2)

        @grid = drawCircleGrid(@ax, @scale, (math.sqrt(surface) for surface in [0..@rmax] by 0.2))
        @directions = drawDirections(@ax, @rmax, @scale, @labels)

        @roses = new Roses(@ax, @scale, @binner.bins, @colorDirection)

    updateGoi: () ->
        @angles = []
        for row in @barycoords
            if row.goi and row.diffexp
                @angles.push(row.angle)

        total = @angles.length
        @binned = @binner.bin(@angles)
        @binned = (math.sqrt(count / total) for count in @binned)

        @roses.updateRadii(@binned)


class Pvalplot
    constructor: (@ax, @w, @h, @scores, @rmax, @labels) ->
        @scale = (x) => x * (@h-0)/(@rmax*2)

        @grid = drawCircleGrid(@ax, @scale, (-math.log10(pval) for pval in [0.1, 0.01, 0.001, 0.0001, 0.00001]))
        @directions = drawDirections(@ax, @rmax, @scale, @labels, 0, "black", 0.025)

        alldotsdata = @generateDotsData()
        @alldots = new Dots(@ax, alldotsdata, @scale)

        tip = d3.tip().html((d)->d.gsetname).attr('class', 'd3-tip').offset([-10, 0])
        @alldots.dots.call(tip)
        @alldots.dots.selectAll("g")
            .on("mouseover", tip.show)
            .on("mouseout", tip.hide)
            .classed("pval", true)

    generateDotsData: (filterflags= {}) ->
        dotsdata = []

        sizeScale = (d) -> (1-d.redundancy)**2 * 5+1

        for row in @scores
            console.log(-row.logqval_unidir)
            if row.logqval_unidir != undefined
                r = math.min(@rmax, -row.logqval_unidir)
                console.log(-row.logqval_unidir)
                angle = row.angle
                newpoint = polar2cart({"angle":angle, "r":r})

                classes = []

                passed = true
                if filterflags.significant and r < -math.log10(0.1)
                    passed = false

                if passed
                    dotsdata.push({
                        "x":newpoint.x,
                        "y":newpoint.y,
                        "angle":angle,
                        "r":r,
                        "gsetid":row.gsetid
                        "gsetname":row.name,
                        "classes":classes,
                        "size":sizeScale(row)
                    })

        return dotsdata


class Dotplot
    constructor: (@ax, @w, @h, @barycoords, @rmax, @labels, @Glabels) ->
        padding = 20
        @scale = (x) => x * (@h-padding)/(@rmax*2)

        @grid = drawHexagonGrid(@ax, @rmax, @scale)
        @directions = drawDirections(@ax, @rmax, @scale, @labels)
        @alldotsdata = @generateDotsData({}, {"diffexp":true})
        @alldots = new Dots(@ax, @alldotsdata, @scale)

        goidotsdata = @generateDotsData({"goi":true}, {"diffexp":true, "goi":true})
        @goidots = new Dots(@ax, goidotsdata, @scale)

        @selectinfo = {
            angle1:0,
            angle2:math.pi/3,
            rmin:4,
            rmax:null
        }

        @pins = @ax.insert("g", ":first-child")
            .classed("pins", true)

        @rings = @ax.insert("g", ":first-child")
            .classed("rings", true)

        @ghover = @ax.append("circle")
            .classed("ghover", true)
            .attr("r", 5)
            .style("visibility", "hidden")

        @logpval_scale = d3.scale.linear()
            .domain(_.range(0, -5, -0.25))
            .range(["#fff5f0","#ffede5","#fee5d8","#fed9c9","#fdcab5","#fcbba1","#fcab8f","#fc9b7c","#fc8a6a","#fb7a5a","#fb694a","#f6583e","#f14432","#e83429","#d92523","#ca181d","#bc141a","#ac1117","#980c13","#7e0610"])
            .clamp(true)

        @opts = {
            visual: {
                goi: {
                    "fill":("#FBB4AE")
                },
                goidiffexp:{
                    "fill":("#E41A1C")
                },
                direction_label:{
                    "font-weight":"bold",
                    "font-size":"1em"
                },
                pin_label:{
                    "font-size":"0.8em"
                }
            }
        }
        @updateVisual()

        tip = d3.
            tip().
            html((d)=>
                @Glabels[d.gid]
            )
            .attr('class', 'd3-tip').offset([-10, 0])

        @alldots.dots.call(tip)
        @goidots.dots.call(tip)

        @alldots.dots.selectAll("g")
            .on("click", (d) =>
                console.log(d)
                if d.gid in @Gpin
                    @updateGpin([], [d.gid])
                else
                    @updateGpin([d.gid])
            )
            .on("mouseover", (d) => @updateHover(d.gid);tip.show(d))
            .on("mouseout", (d) => @updateHover(undefined);tip.hide(d))


        @goidots.dots.selectAll("g")
            .on("click", (d) =>
                console.log(d)
                if d.gid in @Gpin
                    @updateGpin([], [d.gid])
                else
                    @updateGpin([d.gid])
            )
            .on("mouseover", (d) => @updateHover(d.gid);tip.show(d))
            .on("mouseout", (d) => @updateHover(undefined);tip.hide(d))

    updateHover: (gid) ->
        if gid?
            point = @alldotsdata[gid]
            @ghover
                .attr("cx", @scale(point.x))
                .attr("cy", @scale(point.y))
                .style("visibility", "visible")
                .style("pointer-events", "none")
        else
            @ghover
                .style("visibility", "hidden")

    updateRings: (logpvals, deltangle=math.pi/24) ->
        ringsdata = []
        for [angle1, logpval] in _.zip(_.range(-deltangle/2, math.pi*2-deltangle/2, deltangle), logpvals)

            coords = hexagon_path(angle1, angle1+deltangle)

            ringsdata.push({
                angle1:angle1,
                angle2:angle1+deltangle,
                coords:coords,
                logpval:logpval
            })

        r1 = @rmax * 1.02
        r2 = @rmax * 1.08

        @rings.html("")

        @rings.selectAll("polygon")
            .data(ringsdata)
            .enter()
            .append("polygon")
            .attr("points", (d) => ("#{@scale(r1 * point.x)},#{@scale(r1 * point.y)}" for point in d.coords).join(" ") + " " + ("#{@scale(r2 * point.x)},#{@scale(r2 * point.y)}" for point in d.coords.reverse()).join(" "))
            .style("fill", (d) => @logpval_scale(d.logpval))

    generateDotsData: (filterflags= {}, classflags={}) ->
        dotsdata = []

        for point in @barycoords
            if filterflags["goi"] and !(point.goi)
                continue

            gid = point.gid
            ppoint = cart2polar({"x":point.x, "y":point.y})
            newr = math.min(ppoint.r, hexagon_polar(ppoint.angle, @rmax))
            newpoint = polar2cart({"angle":ppoint.angle, "r":newr})

            # determine class
            classes = ["gene"]
            if classflags.goi and (point.goi)
                classes.push("goi")
            else
                classes.push("nogoi")
            if classflags.diffexp and point.diffexp
                classes.push("diffexp")
            else
                classes.push("nodiffexp")
            classes = classes.join(" ")

            dotsdata.push({
                "x":newpoint.x,
                "y":newpoint.y,
                "angle":ppoint.angle,
                "r":newr,
                "gid":gid,
                "diffexp":point.diffexp,
                "classes": classes
            })

        return dotsdata

    updateGoi: () ->
        goidotsdata = @generateDotsData({"goi":true}, {"diffexp":true, "goi":true})
        @goidots.updateData(goidotsdata)
        @updateVisual()

    initGpin: (@Gpin)->
        if @pins?
            @pins.html("")

        @pindata = []
        @updateGpin(@Gpin)

    updateGpin: (add=[], remove=[]) ->
        padding = 1.1

        for newgid in add
            dotdata = @alldotsdata[newgid]
            angle = dotdata.angle
            point = polar2cart({r:hexagon_polar(angle, @rmax * padding), angle:angle})

            edge = hexagon_edge(angle)
            [va, ha] = hexagon_edge_alignment(edge)

            @pindata.push({
                labelx:point.x,
                labely:point.y,
                labelangle:math.mod(angle, math.pi*2),
                x:dotdata.x,
                y:dotdata.y,
                angle:angle,
                gid:newgid,
                label:@Glabels[newgid],
                va:va,
                ha:ha,
                edge:edge
            })

            @Gpin.push(newgid)

        for oldgid in remove
            @pindata.splice(_.findIndex(@pindata, {gid:oldgid}), 1)
            @Gpin.splice(@Gpin.indexOf(oldgid))
        @pindata = _.sortBy(@pindata, (d) -> math.mod(d.angle, math.pi*2)) # sort by angle

        pins = @pins.selectAll("g")
            .data(@pindata, (d)-> d.gid)

        pins.exit()
            .remove()

        pinsenter = pins.enter()
            .append("g")

        pinsenter
            .append("line")

        pinsenter
            .append("text")
            .text((d) -> d.label)
            .attr("x", (d) => @scale(d.labelx))
            .attr("y", (d) => @scale(d.labely))
            .style("baseline-shift", (d) -> d.va)
            .style("text-anchor", (d) -> d.ha)
            .attr("class", "pin_label")
            .on('mouseover', (d)=> @updateHover(d.gid))
            .on('mouseout', (d)=> @updateHover(undefined))
            .on('click', (d)=>@updateGpin([], [d.gid]))

        pins.order() # this changes the order of the elements in the DOM to match the ordering of the data (important for optimizeGpin)

        shrink=false
        if remove.length > 0
            shrink=true
        @optimizeGpin(padding, shrink)

        # add lines
        pins.selectAll("line")
            .attr("x1", (d) => @scale(d.x))
            .attr("y1", (d) => @scale(d.y))
            .attr("x2", (d) => @scale(d.labelx))
            .attr("y2", (d) => @scale(d.labely))
            .style("fill", "none")
            .style("stroke", "#000")
            .style("opacity", 0.6)
            .style("stroke-width", 1)

    optimizeGpin: (padding, shrink=false) ->
        @updateVisual()

        # forcedata = which genes are next to eachother
        forcedata = []
        for i in [0..@pindata.length-1]
            # only add if next label is on the same edge (if not, no pushing from that direction will happen)
            #if hexagon_edge(@pindata[i].labelangle) == hexagon_edge(@pindata[math.mod(i+1, @pindata.length)].labelangle)
                next = math.mod(i+1, @pindata.length)
                forcedata.push([i, next])

        for row in forcedata
            ""
            #console.log(Glabels[@pindata[row[0]].gid] + "->" + Glabels[@pindata[row[1]].gid])
        #console.log(forcedata)

        labelbboxes = (label.getBBox() for label in @directions.selectAll("text")[0])
        for i in [0..1000]
            pinbboxes = (label.getBBox() for label in @pins.selectAll("text")[0])
            moved = 0
            deltangle = (0 for j in [1..@pindata.length])

            if (shrink==true) and (i == 0)
                for j in [0..@pindata.length-1]
                    delta =  difference_circular(@pindata[j].labelangle, @pindata[j].angle)
                    deltangle[j] = posneg(delta) * math.min(0.1, math.abs(delta))
                moved += 1
            else
                for forcerow in forcedata
                    [j, k] = forcerow
                    overlap = bbox_overlap(pinbboxes[j], pinbboxes[k], 1, 1)
                    if overlap# || !angles_increase(@pindata[forcerow[0]].labelangle, @pindata[forcerow[1]].labelangle)
                        #weight = 1 - (@pindata[forcerow[0]].labelangle - @pindata[forcerow[1]].labelangle) * 2 # allow bigger changes if the difference is larger
                        #console.log(weight)
                        #weight=1
                        #if @pindata[forcerow[0]].edge > 2
                        force = 0.01
                        deltangle[j] -= force
                        deltangle[k] += force
                        #if forcerow[0] == goi || forcerow[1] == goi
                            #console.log(force)
                            #console.log((@pindata[forcerow[0]].labelangle - @pindata[forcerow[1]].labelangle))

                        moved += 1
            # force to origin
            #for j in [0..@pindata.length-1]
            #    delta = @pindata[j].labelangle - @pindata[j].angle
            #    deltangle[j] -=
            #console.log(deltangle[goi])

            # forces calculated, now update the true angles and x,y positions
            for j in [0..@pindata.length-1]
                #newangle = hexagon_edgeclip(@pindata[j].labelangle + deltangle[j], @pindata[j].edge, 0.05)
                newangle = @pindata[j].labelangle + deltangle[j]
                newpoint = polar2cart({r:hexagon_polar(newangle, @rmax * padding), angle:newangle})

                edge = hexagon_edge(newangle)
                [va, ha] = hexagon_edge_alignment(edge)

                @pindata[j].labelx = newpoint.x
                @pindata[j].labely = newpoint.y
                @pindata[j].va = va
                @pindata[j].ha = ha
                @pindata[j].labelangle = newangle

            # update text
            @pins.selectAll("g")
                .data(@pindata)
                .select("text")
                .attr("x", (d) => @scale(d.labelx))
                .attr("y", (d) => @scale(d.labely))
                .style("baseline-shift", (d) -> d.va)
                .style("text-anchor", (d) -> d.ha)

            if moved == 0
                console.log("Automatic pin positioning converged after " + i + " iterations")
                break

        if moved > 0
            console.log("Automatic pin positioning not converged! " + moved + " overlaps")

    updateOptions: (newopts) ->
        if "visual" of newopts
            $.extend(@opts.visual, newopts.visual)
            @updateVisual()

    updateVisual: () ->
        @goidots.dots.selectAll("g")
            .style(@opts.visual.goi)
        @goidots.dots.selectAll("g.diffexp")
            .style(@opts.visual.goidiffexp)
        @directions.selectAll("text")
            .style(@opts.visual.direction_label)
        @pins.selectAll("text")
            .style(@opts.visual.pin_label)


class Dots
    constructor: (@ax, @dotsdata, @scale) ->
        @dots = @ax.append("g")
            .classed("dots", true)

        individualdots = @dots
            .selectAll("g")
            .data(@dotsdata)
            .enter()
            .append("g") # add a group here, so we can scale it with origin center of the circle (=> scaling through css instead of the r parameter, will also allow custom shapes etc)
            .attr("class", (d) -> d.classes)
            .attr("transform", (d) => "translate(" + @scale(d.x) + ", " + @scale(d.y) + ")")
            .append("circle")
            .attr("cx", 0)
            .attr("cy", 0)
            .attr("r", 1)
            .attr("transform", (d) ->
                if d.size?
                    return("scale(" + d.size + ")")
                else
                    return("scale(1)")
            )

    updateData: (@dotsdata) ->
        dataselection = @dots
            .selectAll("g")
            .data(@dotsdata)

        dataselection
            .enter()
            .append("g")
            .append("circle")
            .attr("cx", 0)
            .attr("cy", 0)
            .attr("r", 1)

        dataselection
            .exit()
            .remove()

        dataselection
            .attr("class", (d) -> d.classes)
            .attr("transform", (d) => "translate(" + @scale(d.x) + ", " + @scale(d.y) + ")")

class Roses
    constructor: (@ax, @scale, @bins, @colorDirection) ->
        @roseData = (
            {
                "angle1":@bins[math.mod(i-1, @bins.length)],
                "angle2":@bins[i]
            } for i in [0..@bins.length-1]
        )

        @wedges = @ax.append("g")
            .classed("roses", true)
            .selectAll("path")
            .data(@roseData)
            .enter()
            .append("path")
            .attr("d", (d) ->
                r = 1 # scaling of the radius is done through transformations (easier to animate)
                x1 = math.cos(d.angle1) * r
                y1 = math.sin(d.angle1) * r
                x2 = math.cos(d.angle2) * r
                y2 = math.sin(d.angle2) * r

                arc = (d.angle2 - d.angle1) / math.pi * 180

                return "M0,0 L#{x1},#{y1} A#{r},#{r} #{arc} 0,1 #{x2},#{y2} z"
            )
            .style("fill", (d) -> colorDirection(cal_circular_mean([d.angle1, d.angle2])))
    updateRadii: (radii) ->
        for roseDatum, i in @roseData
            roseDatum["radius"] = radii[i]

        @wedges
            .data(@roseData)

        @wedges
            .transition()
            .duration(1000)
            .attr("transform", (d) =>
                r = d.radius
                scaling = @scale(r) + 0.01
                return "scale(#{scaling}, #{scaling})"
            )
