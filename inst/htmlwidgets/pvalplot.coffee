w = 0
h = 0

HTMLWidgets.widget({
    name: 'pvalplot',
    type: 'output',
    initialize: (el, width, height) ->
        w = width
        h = height
    ,
    renderValue: (el, data, instance) ->
        console.log(data)
        window.data = data

        scores = data.scores
        labels = data.labels
        gsetlabels = data.gsetlabels

        ## dot plot
        rmax = 5

        ax = d3.select("div.pvalplot")
            .append("svg")
            .attr("width", w)
            .attr("height", h)
            .append("g")
            .attr("transform", "translate(" + w/2 + "," + h/2 + ")")

        pvalplot = new Pvalplot(ax, w, h, scores, rmax, labels, gsetlabels)

    ,
    resize: (el, width, height, instance) ->

});
