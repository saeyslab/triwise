var Dotplot, Dots, Pvalplot, Roseplot, Roses, drawCircleGrid, drawDirections, drawHexagonGrid,
  indexOf = [].indexOf || function(item) { for (var i = 0, l = this.length; i < l; i++) { if (i in this && this[i] === item) return i; } return -1; };

drawHexagonGrid = function(ax, rmax, scale, rdelta, rmin, anglebase, color, alpha, lw) {
  var angle, grid, hexagon_points, l, m, point, r, ref, ref1, ref2, ref3, ref4, ref5;
  if (rdelta == null) {
    rdelta = 1;
  }
  if (rmin == null) {
    rmin = 0;
  }
  if (anglebase == null) {
    anglebase = 0;
  }
  if (color == null) {
    color = "#444444";
  }
  if (alpha == null) {
    alpha = 0.7;
  }
  if (lw == null) {
    lw = 1;
  }
  grid = ax.append("g");
  hexagon_points = [];
  for (angle = l = ref = anglebase, ref1 = anglebase + Math.PI * 2, ref2 = Math.PI / 6 * 2; ref2 > 0 ? l <= ref1 : l >= ref1; angle = l += ref2) {
    hexagon_points.push(polar2cart({
      "angle": angle,
      "r": 1
    }));
  }
  for (r = m = ref3 = rmin, ref4 = rmax, ref5 = rdelta; ref5 > 0 ? m <= ref4 : m >= ref4; r = m += ref5) {
    grid.append("polygon").attr("points", ((function() {
      var len, n, results;
      results = [];
      for (n = 0, len = hexagon_points.length; n < len; n++) {
        point = hexagon_points[n];
        results.push(scale(point.x) * r + "," + scale(point.y) * r);
      }
      return results;
    })()).join(" ")).style("fill", "none").style("stroke", color).style("opacity", alpha).style("stroke-width", lw);
  }
  return grid;
};

drawCircleGrid = function(ax, scale, radii, rmax, rdelta, rmin, anglebase, color, alpha, lw) {
  var grid, r;
  if (radii == null) {
    radii = null;
  }
  if (rmax == null) {
    rmax = 1;
  }
  if (rdelta == null) {
    rdelta = 1;
  }
  if (rmin == null) {
    rmin = 0;
  }
  if (anglebase == null) {
    anglebase = 0;
  }
  if (color == null) {
    color = "#444444";
  }
  if (alpha == null) {
    alpha = 0.7;
  }
  if (lw == null) {
    lw = 1;
  }
  grid = ax.append("g");
  if (radii === null) {
    radii = (function() {
      var l, ref, ref1, ref2, results;
      results = [];
      for (r = l = ref = rmin, ref1 = rmax, ref2 = rdelta; ref2 > 0 ? l <= ref1 : l >= ref1; r = l += ref2) {
        results.push(r);
      }
      return results;
    })();
  }
  grid.selectAll("circle").data(radii).enter().append("circle").attr("r", function(d) {
    return scale(d);
  }).style("fill", "none").style("stroke", color).style("opacity", alpha).style("stroke-width", lw);
  return grid;
};

drawDirections = function(ax, rmax, scale, labels, anglebase, color, padding, lw) {
  var angle, angles, directions, ha, l, label, len, point, ref, ref1, ref2, va;
  if (anglebase == null) {
    anglebase = 0;
  }
  if (color == null) {
    color = "black";
  }
  if (padding == null) {
    padding = 0.1;
  }
  if (lw == null) {
    lw = 1;
  }
  directions = ax.append("g").classed("directions", true);
  angles = (function() {
    var l, ref, ref1, ref2, results;
    results = [];
    for (angle = l = ref = anglebase, ref1 = anglebase + Math.PI * 2 - 0.001, ref2 = Math.PI / 3 * 2; ref2 > 0 ? l <= ref1 : l >= ref1; angle = l += ref2) {
      results.push(angle);
    }
    return results;
  })();
  ref = _.zip(angles, labels);
  for (l = 0, len = ref.length; l < len; l++) {
    ref1 = ref[l], angle = ref1[0], label = ref1[1];
    point = polar2cart({
      "angle": angle,
      "r": rmax
    });
    directions.append("line").attr({
      x1: 0,
      y1: 0,
      x2: scale(point.x),
      y2: scale(point.y)
    }).style("stroke-width", lw).style("stroke", color);
    point = polar2cart({
      "angle": angle,
      "r": rmax * (padding + 1)
    });
    if (point.y < 0) {
      va = "0%";
    } else if (point.y === 0) {
      va = "-30%";
    } else {
      va = "-75%";
    }
    if ((scale(-0.1) < (ref2 = scale(point.x)) && ref2 < scale(0.1))) {
      ha = "middle";
    } else if (point.x < 0) {
      ha = "end";
    } else {
      ha = "start";
    }
    directions.append("text").attr("x", scale(point.x)).attr("y", scale(point.y)).text(label).style("baseline-shift", va).style("text-anchor", ha).attr("class", "direction_label").style("font-weight", "bold");
  }
  return directions;
};

repositionDirections(directions, bboxes)(function() {
  return console.log("nope");
});

Roseplot = (function() {
  function Roseplot(ax1, w, h, barycoords, labels1, binner, colorDirection1) {
    var surface;
    this.ax = ax1;
    this.w = w;
    this.h = h;
    this.barycoords = barycoords;
    this.labels = labels1;
    this.binner = binner;
    this.colorDirection = colorDirection1;
    this.rmax = 1;
    this.scale = (function(_this) {
      return function(x) {
        return x * (_this.h - 0) / (_this.rmax * 2);
      };
    })(this);
    this.grid = drawCircleGrid(this.ax, this.scale, (function() {
      var l, ref, results;
      results = [];
      for (surface = l = 0, ref = this.rmax; l <= ref; surface = l += 0.2) {
        results.push(math.sqrt(surface));
      }
      return results;
    }).call(this));
    this.directions = drawDirections(this.ax, this.rmax, this.scale, this.labels);
    this.roses = new Roses(this.ax, this.scale, this.binner.bins, this.colorDirection);
  }

  Roseplot.prototype.updateGoi = function() {
    var count, l, len, ref, row, total;
    this.angles = [];
    ref = this.barycoords;
    for (l = 0, len = ref.length; l < len; l++) {
      row = ref[l];
      if (row.goi && row.diffexp) {
        this.angles.push(row.angle);
      }
    }
    total = this.angles.length;
    this.binned = this.binner.bin(this.angles);
    this.binned = (function() {
      var len1, m, ref1, results;
      ref1 = this.binned;
      results = [];
      for (m = 0, len1 = ref1.length; m < len1; m++) {
        count = ref1[m];
        results.push(math.sqrt(count / total));
      }
      return results;
    }).call(this);
    return this.roses.updateRadii(this.binned);
  };

  return Roseplot;

})();

Pvalplot = (function() {
  function Pvalplot(ax1, w, h, scores, rmax1, labels1) {
    var alldotsdata, pval, tip;
    this.ax = ax1;
    this.w = w;
    this.h = h;
    this.scores = scores;
    this.rmax = rmax1;
    this.labels = labels1;
    this.scale = (function(_this) {
      return function(x) {
        return x * (_this.h - 0) / (_this.rmax * 2);
      };
    })(this);
    this.grid = drawCircleGrid(this.ax, this.scale, (function() {
      var l, len, ref, results;
      ref = [0.1, 0.01, 0.001, 0.0001, 0.00001];
      results = [];
      for (l = 0, len = ref.length; l < len; l++) {
        pval = ref[l];
        results.push(-math.log10(pval));
      }
      return results;
    })());
    this.directions = drawDirections(this.ax, this.rmax, this.scale, this.labels, 0, "black", 0.025);
    alldotsdata = this.generateDotsData();
    this.alldots = new Dots(this.ax, alldotsdata, this.scale);
    tip = d3.tip().html(function(d) {
      return d.gsetname;
    }).attr('class', 'd3-tip').offset([-10, 0]);
    this.alldots.dots.call(tip);
    this.alldots.dots.selectAll("g").on("mouseover", tip.show).on("mouseout", tip.hide).classed("pval", true);
  }

  Pvalplot.prototype.generateDotsData = function(filterflags) {
    var angle, classes, dotsdata, l, len, newpoint, passed, r, ref, row, sizeScale;
    if (filterflags == null) {
      filterflags = {};
    }
    dotsdata = [];
    sizeScale = function(d) {
      return Math.pow(1 - d.redundancy, 2) * 5 + 1;
    };
    ref = this.scores;
    for (l = 0, len = ref.length; l < len; l++) {
      row = ref[l];
      console.log(-row.logqval_unidir);
      if (row.logqval_unidir !== void 0) {
        r = math.min(this.rmax, -row.logqval_unidir);
        console.log(-row.logqval_unidir);
        angle = row.angle;
        newpoint = polar2cart({
          "angle": angle,
          "r": r
        });
        classes = [];
        passed = true;
        if (filterflags.significant && r < -math.log10(0.1)) {
          passed = false;
        }
        if (passed) {
          dotsdata.push({
            "x": newpoint.x,
            "y": newpoint.y,
            "angle": angle,
            "r": r,
            "gsetid": row.gsetid,
            "gsetname": row.name,
            "classes": classes,
            "size": sizeScale(row)
          });
        }
      }
    }
    return dotsdata;
  };

  return Pvalplot;

})();

Dotplot = (function() {
  function Dotplot(ax1, w, h, barycoords, rmax1, labels1, Glabels) {
    var goidotsdata, padding, tip;
    this.ax = ax1;
    this.w = w;
    this.h = h;
    this.barycoords = barycoords;
    this.rmax = rmax1;
    this.labels = labels1;
    this.Glabels = Glabels;
    padding = 20;
    this.scale = (function(_this) {
      return function(x) {
        return x * (_this.h - padding) / (_this.rmax * 2);
      };
    })(this);
    this.grid = drawHexagonGrid(this.ax, this.rmax, this.scale);
    this.directions = drawDirections(this.ax, this.rmax, this.scale, this.labels);
    this.alldotsdata = this.generateDotsData({}, {
      "diffexp": true
    });
    this.alldots = new Dots(this.ax, this.alldotsdata, this.scale);
    goidotsdata = this.generateDotsData({
      "goi": true
    }, {
      "diffexp": true,
      "goi": true
    });
    this.goidots = new Dots(this.ax, goidotsdata, this.scale);
    this.selectinfo = {
      angle1: 0,
      angle2: math.pi / 3,
      rmin: 4,
      rmax: null
    };
    this.pins = this.ax.insert("g", ":first-child").classed("pins", true);
    this.rings = this.ax.insert("g", ":first-child").classed("rings", true);
    this.ghover = this.ax.append("circle").classed("ghover", true).attr("r", 5).style("visibility", "hidden");
    this.logpval_scale = d3.scale.linear().domain(_.range(0, -5, -0.25)).range(["#fff5f0", "#ffede5", "#fee5d8", "#fed9c9", "#fdcab5", "#fcbba1", "#fcab8f", "#fc9b7c", "#fc8a6a", "#fb7a5a", "#fb694a", "#f6583e", "#f14432", "#e83429", "#d92523", "#ca181d", "#bc141a", "#ac1117", "#980c13", "#7e0610"]).clamp(true);
    this.opts = {
      visual: {
        goi: {
          "fill": "#FBB4AE"
        },
        goidiffexp: {
          "fill": "#E41A1C"
        },
        direction_label: {
          "font-weight": "bold",
          "font-size": "1em"
        },
        pin_label: {
          "font-size": "0.8em"
        }
      }
    };
    this.updateVisual();
    tip = d3.tip().html((function(_this) {
      return function(d) {
        return _this.Glabels[d.gid];
      };
    })(this)).attr('class', 'd3-tip').offset([-10, 0]);
    this.alldots.dots.call(tip);
    this.alldots.dots.selectAll("g").on("mouseover", tip.show).on("mouseout", tip.hide);
    this.goidots.dots.call(tip);
    this.goidots.dots.selectAll("g").on("mouseover", tip.show).on("mouseout", tip.hide);
    this.alldots.dots.selectAll("g").on("click", (function(_this) {
      return function(d) {
        var ref;
        console.log(d);
        if (ref = d.gid, indexOf.call(_this.Gpin, ref) >= 0) {
          return _this.updateGpin([], [d.gid]);
        } else {
          return _this.updateGpin([d.gid]);
        }
      };
    })(this));
    this.goidots.dots.selectAll("g").on("click", (function(_this) {
      return function(d) {
        var ref;
        console.log(d);
        if (ref = d.gid, indexOf.call(_this.Gpin, ref) >= 0) {
          return _this.updateGpin([], [d.gid]);
        } else {
          return _this.updateGpin([d.gid]);
        }
      };
    })(this));
  }

  Dotplot.prototype.updateHover = function(gid) {
    var point;
    if (gid != null) {
      point = this.alldotsdata[gid];
      return this.ghover.attr("cx", this.scale(point.x)).attr("cy", this.scale(point.y)).style("visibility", "visible");
    } else {
      return this.ghover.style("visibility", "hidden");
    }
  };

  Dotplot.prototype.updateRings = function(logpvals, deltangle) {
    var angle1, coords, l, len, logpval, r1, r2, ref, ref1, ringsdata;
    if (deltangle == null) {
      deltangle = math.pi / 24;
    }
    ringsdata = [];
    ref = _.zip(_.range(-deltangle / 2, math.pi * 2 - deltangle / 2, deltangle), logpvals);
    for (l = 0, len = ref.length; l < len; l++) {
      ref1 = ref[l], angle1 = ref1[0], logpval = ref1[1];
      coords = hexagon_path(angle1, angle1 + deltangle);
      ringsdata.push({
        angle1: angle1,
        angle2: angle1 + deltangle,
        coords: coords,
        logpval: logpval
      });
    }
    r1 = this.rmax * 1.02;
    r2 = this.rmax * 1.08;
    this.rings.html("");
    return this.rings.selectAll("polygon").data(ringsdata).enter().append("polygon").attr("points", (function(_this) {
      return function(d) {
        var point;
        return ((function() {
          var len1, m, ref2, results;
          ref2 = d.coords;
          results = [];
          for (m = 0, len1 = ref2.length; m < len1; m++) {
            point = ref2[m];
            results.push((this.scale(r1 * point.x)) + "," + (this.scale(r1 * point.y)));
          }
          return results;
        }).call(_this)).join(" ") + " " + ((function() {
          var len1, m, ref2, results;
          ref2 = d.coords.reverse();
          results = [];
          for (m = 0, len1 = ref2.length; m < len1; m++) {
            point = ref2[m];
            results.push((this.scale(r2 * point.x)) + "," + (this.scale(r2 * point.y)));
          }
          return results;
        }).call(_this)).join(" ");
      };
    })(this)).style("fill", (function(_this) {
      return function(d) {
        return _this.logpval_scale(d.logpval);
      };
    })(this));
  };

  Dotplot.prototype.generateDotsData = function(filterflags, classflags) {
    var classes, dotsdata, gid, l, len, newpoint, newr, point, ppoint, ref;
    if (filterflags == null) {
      filterflags = {};
    }
    if (classflags == null) {
      classflags = {};
    }
    dotsdata = [];
    ref = this.barycoords;
    for (l = 0, len = ref.length; l < len; l++) {
      point = ref[l];
      if (filterflags["goi"] && !point.goi) {
        continue;
      }
      gid = point.gid;
      ppoint = cart2polar({
        "x": point.x,
        "y": point.y
      });
      newr = math.min(ppoint.r, hexagon_polar(ppoint.angle, this.rmax));
      newpoint = polar2cart({
        "angle": ppoint.angle,
        "r": newr
      });
      classes = ["gene"];
      if (classflags.goi && point.goi) {
        classes.push("goi");
      } else {
        classes.push("nogoi");
      }
      if (classflags.diffexp && point.diffexp) {
        classes.push("diffexp");
      } else {
        classes.push("nodiffexp");
      }
      classes = classes.join(" ");
      dotsdata.push({
        "x": newpoint.x,
        "y": newpoint.y,
        "angle": ppoint.angle,
        "r": newr,
        "gid": gid,
        "diffexp": point.diffexp,
        "classes": classes
      });
    }
    return dotsdata;
  };

  Dotplot.prototype.updateGoi = function() {
    var goidotsdata;
    goidotsdata = this.generateDotsData({
      "goi": true
    }, {
      "diffexp": true,
      "goi": true
    });
    this.goidots.updateData(goidotsdata);
    return this.updateVisual();
  };

  Dotplot.prototype.initGpin = function(Gpin) {
    this.Gpin = Gpin;
    if (this.pins != null) {
      this.pins.html("");
    }
    this.pindata = [];
    return this.updateGpin(this.Gpin);
  };

  Dotplot.prototype.updateGpin = function(add, remove) {
    var angle, dotdata, edge, ha, l, len, len1, m, newgid, oldgid, padding, pins, pinsenter, point, ref, shrink, va;
    if (add == null) {
      add = [];
    }
    if (remove == null) {
      remove = [];
    }
    padding = 1.1;
    for (l = 0, len = add.length; l < len; l++) {
      newgid = add[l];
      dotdata = this.alldotsdata[newgid];
      angle = dotdata.angle;
      point = polar2cart({
        r: hexagon_polar(angle, this.rmax * padding),
        angle: angle
      });
      edge = hexagon_edge(angle);
      ref = hexagon_edge_alignment(edge), va = ref[0], ha = ref[1];
      this.pindata.push({
        labelx: point.x,
        labely: point.y,
        labelangle: math.mod(angle, math.pi * 2),
        x: dotdata.x,
        y: dotdata.y,
        angle: angle,
        gid: newgid,
        label: this.Glabels[newgid],
        va: va,
        ha: ha,
        edge: edge
      });
      this.Gpin.push(newgid);
    }
    for (m = 0, len1 = remove.length; m < len1; m++) {
      oldgid = remove[m];
      this.pindata.splice(_.findIndex(this.pindata, {
        gid: oldgid
      }), 1);
      this.Gpin.splice(this.Gpin.indexOf(oldgid));
    }
    this.pindata = _.sortBy(this.pindata, function(d) {
      return math.mod(d.angle, math.pi * 2);
    });
    pins = this.pins.selectAll("g").data(this.pindata, function(d) {
      return d.gid;
    });
    pins.exit().remove();
    pinsenter = pins.enter().append("g");
    pinsenter.append("line");
    pinsenter.append("text").text(function(d) {
      return d.label;
    }).attr("x", (function(_this) {
      return function(d) {
        return _this.scale(d.labelx);
      };
    })(this)).attr("y", (function(_this) {
      return function(d) {
        return _this.scale(d.labely);
      };
    })(this)).style("baseline-shift", function(d) {
      return d.va;
    }).style("text-anchor", function(d) {
      return d.ha;
    }).attr("class", "pin_label").on('mouseover', (function(_this) {
      return function(d) {
        return _this.updateHover(d.gid);
      };
    })(this)).on('mouseout', (function(_this) {
      return function(d) {
        return _this.updateHover(void 0);
      };
    })(this));
    pins.order();
    shrink = false;
    if (remove.length > 0) {
      shrink = true;
    }
    this.optimizeGpin(padding, shrink);
    return pins.selectAll("line").attr("x1", (function(_this) {
      return function(d) {
        return _this.scale(d.x);
      };
    })(this)).attr("y1", (function(_this) {
      return function(d) {
        return _this.scale(d.y);
      };
    })(this)).attr("x2", (function(_this) {
      return function(d) {
        return _this.scale(d.labelx);
      };
    })(this)).attr("y2", (function(_this) {
      return function(d) {
        return _this.scale(d.labely);
      };
    })(this)).style("fill", "none").style("stroke", "#000").style("opacity", 0.6).style("stroke-width", 1);
  };

  Dotplot.prototype.optimizeGpin = function(padding, shrink) {
    var delta, deltangle, edge, force, forcedata, forcerow, ha, i, j, k, l, label, labelbboxes, len, len1, m, moved, n, newangle, newpoint, next, o, overlap, p, pinbboxes, q, ref, ref1, ref2, ref3, row, va;
    if (shrink == null) {
      shrink = false;
    }
    this.updateVisual();
    forcedata = [];
    for (i = l = 0, ref = this.pindata.length - 1; 0 <= ref ? l <= ref : l >= ref; i = 0 <= ref ? ++l : --l) {
      next = math.mod(i + 1, this.pindata.length);
      forcedata.push([i, next]);
    }
    for (m = 0, len = forcedata.length; m < len; m++) {
      row = forcedata[m];
      "";
    }
    labelbboxes = (function() {
      var len1, n, ref1, results;
      ref1 = this.directions.selectAll("text")[0];
      results = [];
      for (n = 0, len1 = ref1.length; n < len1; n++) {
        label = ref1[n];
        results.push(label.getBBox());
      }
      return results;
    }).call(this);
    for (i = n = 0; n <= 1000; i = ++n) {
      pinbboxes = (function() {
        var len1, o, ref1, results;
        ref1 = this.pins.selectAll("text")[0];
        results = [];
        for (o = 0, len1 = ref1.length; o < len1; o++) {
          label = ref1[o];
          results.push(label.getBBox());
        }
        return results;
      }).call(this);
      moved = 0;
      deltangle = (function() {
        var o, ref1, results;
        results = [];
        for (j = o = 1, ref1 = this.pindata.length; 1 <= ref1 ? o <= ref1 : o >= ref1; j = 1 <= ref1 ? ++o : --o) {
          results.push(0);
        }
        return results;
      }).call(this);
      if ((shrink === true) && (i === 0)) {
        for (j = o = 0, ref1 = this.pindata.length - 1; 0 <= ref1 ? o <= ref1 : o >= ref1; j = 0 <= ref1 ? ++o : --o) {
          delta = difference_circular(this.pindata[j].labelangle, this.pindata[j].angle);
          deltangle[j] = posneg(delta) * math.min(0.1, math.abs(delta));
        }
        moved += 1;
      } else {
        for (p = 0, len1 = forcedata.length; p < len1; p++) {
          forcerow = forcedata[p];
          j = forcerow[0], k = forcerow[1];
          overlap = bbox_overlap(pinbboxes[j], pinbboxes[k], 1, 1);
          if (overlap) {
            force = 0.01;
            deltangle[j] -= force;
            deltangle[k] += force;
            moved += 1;
          }
        }
      }
      for (j = q = 0, ref2 = this.pindata.length - 1; 0 <= ref2 ? q <= ref2 : q >= ref2; j = 0 <= ref2 ? ++q : --q) {
        newangle = this.pindata[j].labelangle + deltangle[j];
        newpoint = polar2cart({
          r: hexagon_polar(newangle, this.rmax * padding),
          angle: newangle
        });
        edge = hexagon_edge(newangle);
        ref3 = hexagon_edge_alignment(edge), va = ref3[0], ha = ref3[1];
        this.pindata[j].labelx = newpoint.x;
        this.pindata[j].labely = newpoint.y;
        this.pindata[j].va = va;
        this.pindata[j].ha = ha;
        this.pindata[j].labelangle = newangle;
      }
      this.pins.selectAll("g").data(this.pindata).select("text").attr("x", (function(_this) {
        return function(d) {
          return _this.scale(d.labelx);
        };
      })(this)).attr("y", (function(_this) {
        return function(d) {
          return _this.scale(d.labely);
        };
      })(this)).style("baseline-shift", function(d) {
        return d.va;
      }).style("text-anchor", function(d) {
        return d.ha;
      });
      if (moved === 0) {
        console.log("Automatic pin positioning converged after " + i + " iterations");
        break;
      }
    }
    if (moved > 0) {
      return console.log("Automatic pin positioning not converged! " + moved + " overlaps");
    }
  };

  Dotplot.prototype.updateOptions = function(newopts) {
    if ("visual" in newopts) {
      $.extend(this.opts.visual, newopts.visual);
      return this.updateVisual();
    }
  };

  Dotplot.prototype.updateVisual = function() {
    this.goidots.dots.selectAll("g").style(this.opts.visual.goi);
    this.goidots.dots.selectAll("g.diffexp").style(this.opts.visual.goidiffexp);
    this.directions.selectAll("text").style(this.opts.visual.direction_label);
    return this.pins.selectAll("text").style(this.opts.visual.pin_label);
  };

  return Dotplot;

})();

Dots = (function() {
  function Dots(ax1, dotsdata1, scale1) {
    var individualdots;
    this.ax = ax1;
    this.dotsdata = dotsdata1;
    this.scale = scale1;
    this.dots = this.ax.append("g").classed("dots", true);
    individualdots = this.dots.selectAll("g").data(this.dotsdata).enter().append("g").attr("class", function(d) {
      return d.classes;
    }).attr("transform", (function(_this) {
      return function(d) {
        return "translate(" + _this.scale(d.x) + ", " + _this.scale(d.y) + ")";
      };
    })(this)).append("circle").attr("cx", 0).attr("cy", 0).attr("r", 1).attr("transform", function(d) {
      if (d.size != null) {
        return "scale(" + d.size + ")";
      } else {
        return "scale(1)";
      }
    });
  }

  Dots.prototype.updateData = function(dotsdata1) {
    var dataselection;
    this.dotsdata = dotsdata1;
    dataselection = this.dots.selectAll("g").data(this.dotsdata);
    dataselection.enter().append("g").append("circle").attr("cx", 0).attr("cy", 0).attr("r", 1);
    dataselection.exit().remove();
    return dataselection.attr("class", function(d) {
      return d.classes;
    }).attr("transform", (function(_this) {
      return function(d) {
        return "translate(" + _this.scale(d.x) + ", " + _this.scale(d.y) + ")";
      };
    })(this));
  };

  return Dots;

})();

Roses = (function() {
  function Roses(ax1, scale1, bins, colorDirection1) {
    var i;
    this.ax = ax1;
    this.scale = scale1;
    this.bins = bins;
    this.colorDirection = colorDirection1;
    this.roseData = (function() {
      var l, ref, results;
      results = [];
      for (i = l = 0, ref = this.bins.length - 1; 0 <= ref ? l <= ref : l >= ref; i = 0 <= ref ? ++l : --l) {
        results.push({
          "angle1": this.bins[math.mod(i - 1, this.bins.length)],
          "angle2": this.bins[i]
        });
      }
      return results;
    }).call(this);
    this.wedges = this.ax.append("g").classed("roses", true).selectAll("path").data(this.roseData).enter().append("path").attr("d", function(d) {
      var arc, r, x1, x2, y1, y2;
      r = 1;
      x1 = math.cos(d.angle1) * r;
      y1 = math.sin(d.angle1) * r;
      x2 = math.cos(d.angle2) * r;
      y2 = math.sin(d.angle2) * r;
      arc = (d.angle2 - d.angle1) / math.pi * 180;
      return "M0,0 L" + x1 + "," + y1 + " A" + r + "," + r + " " + arc + " 0,1 " + x2 + "," + y2 + " z";
    }).style("fill", function(d) {
      return colorDirection(cal_circular_mean([d.angle1, d.angle2]));
    });
  }

  Roses.prototype.updateRadii = function(radii) {
    var i, l, len, ref, roseDatum;
    ref = this.roseData;
    for (i = l = 0, len = ref.length; l < len; i = ++l) {
      roseDatum = ref[i];
      roseDatum["radius"] = radii[i];
    }
    this.wedges.data(this.roseData);
    return this.wedges.transition().duration(1000).attr("transform", (function(_this) {
      return function(d) {
        var r, scaling;
        r = d.radius;
        scaling = _this.scale(r) + 0.01;
        return "scale(" + scaling + ", " + scaling + ")";
      };
    })(this));
  };

  return Roses;

})();
