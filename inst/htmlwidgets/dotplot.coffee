w=0
h=0
HTMLWidgets.widget({
    name: 'dotplot',
    type: 'output',
    initialize: (el, width, height) ->
        w = width
        h = height
    ,
    renderValue: (el, data, instance) ->
        window.data = data

        Eoi = data.Eoi
        Gdiffexp = data.Gdiffexp

        G = Eoi.columns
        C = Eoi.index
        labels = C

        Glabels = data.Glabels

        Goi = data.Goi
        Gpin = data.Gpin
        logpvals = data.logpvals

        options = data.options

        rmax = data.rmax

        ## add searchbar
        d3.select("div.dotplot")
            .append("div")
            .style("position", "relative")
            .style("left", "10px")
            .style("top", "50px")
            .classed("awesomplete", true)
            .append("input")
            .attr("id", "genesearch")
            .attr("placeholder", "search gene...")

        genesearch_input = document.getElementById("genesearch")
        genesearch = new Awesomplete(genesearch_input, {list:_.zip(Glabels, G), minChars:1, filter:Awesomplete.FILTER_STARTSWITH, autoFirst:true})

        ## add dotplot

        ## transform to barycentric coordinates
        anglebase = 0
        transformation = [
            [math.cos(0+anglebase),math.cos(math.pi*2/3+anglebase),math.cos(math.pi*2/3+anglebase)],
            [math.sin(0+anglebase), math.sin(math.pi*2/3+anglebase), math.sin(-math.pi*2/3+anglebase)]
        ]

        barycoords = transform_barycentric(Eoi.data, transformation)
        barycoords = barycoords.map((point, gid) ->
            ppoint = cart2polar({x:point[0], y:point[1]})
            {"x":point[0], "y":point[1], "r":ppoint.r, "angle":ppoint.angle,  "gid":gid, "goi":false, "diffexp":Gdiffexp[gid]}
        )

        for row in barycoords
            if row.gid in Goi
                row.goi = true
            else
                row.goi = false

        ## dot plot
        scale = (x) -> x * (h-0)/(rmax*2)

        ax = d3.select("div.dotplot")
            .append("svg")
            .attr("width", w)
            .attr("height", h)
            .append("g")
            .attr("transform", "translate(" + w/2 + "," + h/2 + ")")

        dotplot = new Dotplot(ax, w, h, barycoords, rmax, labels, Glabels)
        dotplot.updateGoi()
        dotplot.initGpin(Gpin)

        console.log(data)
        if options.plotLocalEnrichment[0]
            dotplot.updateRings(logpvals)

        # svg = document.querySelector( "svg" )
        # svgData = new XMLSerializer().serializeToString( svg )
        #
        # canvas = document.createElement( "canvas" )
        # svgSize = svg.getBoundingClientRect();
        # canvas.width = svgSize.width;
        # canvas.height = svgSize.height;
        # ctx = canvas.getContext( "2d" );
        #
        # img = document.createElement("img");
        # img.setAttribute( "src", "data:image/svg+xml;base64," + btoa( svgData ) );
        #
        # img.onload = () ->
        #     ctx.drawImage( img, 0, 0 );
        #     console.log( canvas.toDataURL( "image/png" ) )
        #     window.open(canvas.toDataURL( "image/png" ));

        window.dotplot = dotplot

        genesearch_input.addEventListener("awesomplete-selectcomplete", (e) ->
            dotplot.updateHover(this.value)
            if this.value not in dotplot.Gpin
                dotplot.updateGpin([this.value])
            this.value = Glabels[this.value]
        )

        insetwidth = w/3
        insetheight = h/3

        ax = d3.select("div.dotplot")
            .select("svg")
            .append("g")
            .attr("transform", "translate(" + insetwidth + "," + (h-insetheight) + ")")
        demo = new Demo(dotplot, ax)

    ,
    resize: (el, width, height, instance) ->

});
