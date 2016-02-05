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
    var C, Eoi, G, Gdiffexp, Glabels, Goi, Gpin, anglebase, ax, barycoords, canvas, ctx, dotplot, i, img, labels, len, logpvals, ref, rmax, row, scale, svg, svgData, svgSize, transformation;
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
    dotplot.updateRings(logpvals);
    svg = document.querySelector("svg");
    svgData = new XMLSerializer().serializeToString(svg);
    canvas = document.createElement("canvas");
    svgSize = svg.getBoundingClientRect();
    canvas.width = svgSize.width;
    canvas.height = svgSize.height;
    ctx = canvas.getContext("2d");
    img = document.createElement("img");
    img.setAttribute("src", "data:image/svg+xml;base64," + btoa(svgData));
    return img.onload = function() {
      ctx.drawImage(img, 0, 0);
      return console.log(canvas.toDataURL("image/png"));
    };
  },
  resize: function(el, width, height, instance) {}
});
