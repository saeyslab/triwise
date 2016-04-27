var h, w,
  indexOf = [].indexOf || function(item) { for (var i = 0, l = this.length; i < l; i++) { if (i in this && this[i] === item) return i; } return -1; };

w = 0;

h = 0;

HTMLWidgets.widget({
  name: 'dotplot',
  type: 'output',
  initialize: function(el, width, height) {
    w = width;
    return h = height;
  },
  renderValue: function(el, data, instance) {
    var C, Eoi, G, Gdiffexp, Glabels, Goi, Gpin, anglebase, ax, barycoords, dotplot, genesearch, genesearch_input, i, labels, len, logpvals, ref, rmax, row, scale, transformation;
    window.data = data;
    Eoi = data.Eoi;
    Gdiffexp = data.Gdiffexp;
    G = Eoi.columns;
    C = Eoi.index;
    labels = C;
    Glabels = data.Glabels;
    Goi = data.Goi;
    Gpin = data.Gpin;
    logpvals = data.logpvals;
    d3.select("div.dotplot").append("div").style("position", "relative").style("left", "10px").style("top", "50px").classed("awesomplete", true).append("input").attr("id", "genesearch").attr("placeholder", "search gene...");
    genesearch_input = document.getElementById("genesearch");
    genesearch = new Awesomplete(genesearch_input, {
      list: _.zip(Glabels, G),
      minChars: 1,
      filter: Awesomplete.FILTER_STARTSWITH,
      autoFirst: true
    });
    anglebase = 0;
    transformation = [[math.cos(0 + anglebase), math.cos(math.pi * 2 / 3 + anglebase), math.cos(math.pi * 2 / 3 + anglebase)], [math.sin(0 + anglebase), math.sin(math.pi * 2 / 3 + anglebase), math.sin(-math.pi * 2 / 3 + anglebase)]];
    barycoords = transform_barycentric(Eoi.data, transformation);
    barycoords = barycoords.map(function(point, gid) {
      var ppoint;
      ppoint = cart2polar({
        x: point[0],
        y: point[1]
      });
      return {
        "x": point[0],
        "y": point[1],
        "r": ppoint.r,
        "angle": ppoint.angle,
        "gid": gid,
        "goi": false,
        "diffexp": Gdiffexp[gid]
      };
    });
    for (i = 0, len = barycoords.length; i < len; i++) {
      row = barycoords[i];
      if (ref = row.gid, indexOf.call(Goi, ref) >= 0) {
        row.goi = true;
      } else {
        row.goi = false;
      }
    }
    rmax = 4;
    scale = function(x) {
      return x * (h - 0) / (rmax * 2);
    };
    ax = d3.select("div.dotplot").append("svg").attr("width", w).attr("height", h).append("g").attr("transform", "translate(" + w / 2 + "," + h / 2 + ")");
    dotplot = new Dotplot(ax, w, h, barycoords, rmax, labels, Glabels);
    dotplot.updateGoi();
    dotplot.initGpin(Gpin);
    if (data.plotLocalEnrichment) {
      dotplot.updateRings(logpvals);
    }
    window.dotplot = dotplot;
    return genesearch_input.addEventListener("awesomplete-selectcomplete", function(e) {
      var ref1;
      dotplot.updateHover(this.value);
      if (ref1 = this.value, indexOf.call(dotplot.Gpin, ref1) < 0) {
        dotplot.updateGpin([this.value]);
      }
      return this.value = Glabels[this.value];
    });
  },
  resize: function(el, width, height, instance) {}
});
