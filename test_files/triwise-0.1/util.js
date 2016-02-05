
/* UTIL */
var Binner, angles_increase, bbox_overlap, between_circular, cal_circular_mean, cart2polar, deg2rad, difference_circular, hexagon_edge, hexagon_edge_alignment, hexagon_edgeclip, hexagon_path, hexagon_polar, polar2cart, posneg, rad2deg;

cart2polar = function(point) {
  return {
    "r": Math.sqrt(Math.pow(point.x, 2) + Math.pow(point.y, 2)),
    "angle": Math.atan2(point.y, point.x)
  };
};

polar2cart = function(ppoint) {
  return {
    "x": Math.cos(ppoint.angle) * ppoint.r,
    "y": Math.sin(ppoint.angle) * ppoint.r
  };
};

deg2rad = function(deg) {
  return deg * Math.PI / 180;
};

rad2deg = function(rad) {
  return rad * 180 / Math.PI;
};

hexagon_polar = function(angle, radius) {
  var delta;
  if (radius == null) {
    radius = 1;
  }
  delta = 2 * Math.PI / 6;
  return Math.cos(delta / 2) / Math.cos((math.mod(angle, delta)) - delta / 2) * radius;
};

hexagon_edge = function(angle) {
  var delta;
  delta = 2 * Math.PI / 6;
  return math.mod(Math.floor(angle / delta), 6);
};

hexagon_edgeclip = function(angle, edgeid, padding) {
  var delta, maxangle, minangle;
  if (padding == null) {
    padding = 0;
  }
  delta = 2 * Math.PI / 6;
  minangle = edgeid * delta + padding;
  maxangle = (edgeid + 1) * delta - padding;
  if (math.mod(angle, 2 * Math.PI) > maxangle) {
    return maxangle;
  } else if (math.mod(angle, 2 * Math.PI) < minangle) {
    return minangle;
  } else {
    return angle;
  }
};

hexagon_path = function(angle1, angle2) {
  var angle, angles, coords, j, ref, ref1;
  angles = [angle1];
  for (angle = j = 0, ref = math.pi * 2 - 0.001, ref1 = math.pi / 3; ref1 > 0 ? j <= ref : j >= ref; angle = j += ref1) {
    if (between_circular(angle, angle1, angle2)) {
      angles.push(angle);
    }
  }
  angles.push(angle2);
  coords = (function() {
    var k, len, results;
    results = [];
    for (k = 0, len = angles.length; k < len; k++) {
      angle = angles[k];
      results.push(polar2cart({
        angle: angle,
        r: hexagon_polar(angle, 1)
      }));
    }
    return results;
  })();
  return coords;
};

hexagon_edge_alignment = function(edge) {
  var ha, va;
  if (edge === 1) {
    va = "-100%";
  } else if (edge === 4) {
    va = "0%";
  } else {
    va = "-50%";
  }
  if (edge === 0 || edge === 5) {
    ha = "start";
  } else if (edge === 1 || edge === 4) {
    ha = "middle";
  } else {
    ha = "end";
  }
  return [va, ha];
};

Binner = (function() {
  function Binner(nbins) {
    var i;
    this.nbins = nbins;
    this.delta = math.pi * 2 / this.nbins;
    this.bins = (function() {
      var j, ref, results;
      results = [];
      for (i = j = 0, ref = this.nbins - 1; 0 <= ref ? j <= ref : j >= ref; i = 0 <= ref ? ++j : --j) {
        results.push(this.delta / 2 + i * this.delta);
      }
      return results;
    }).call(this);
  }

  Binner.prototype.bin = function(angles) {
    var angle, binid, binned, i, j, len;
    binned = (function() {
      var j, ref, results;
      results = [];
      for (i = j = 0, ref = this.nbins - 1; 0 <= ref ? j <= ref : j >= ref; i = 0 <= ref ? ++j : --j) {
        results.push(0);
      }
      return results;
    }).call(this);
    for (j = 0, len = angles.length; j < len; j++) {
      angle = angles[j];
      binid = math.mod(math.floor((angle + this.delta / 2) / this.delta), this.nbins);
      binned[binid] += 1;
    }
    return binned;
  };

  Binner.prototype.group = function(barycoords) {
    var binid, grouped, i, j, len, row;
    grouped = (function() {
      var j, ref, results;
      results = [];
      for (i = j = 0, ref = this.nbins - 1; 0 <= ref ? j <= ref : j >= ref; i = 0 <= ref ? ++j : --j) {
        results.push([]);
      }
      return results;
    }).call(this);
    for (j = 0, len = barycoords.length; j < len; j++) {
      row = barycoords[j];
      binid = math.mod(math.floor((row.angle + this.delta / 2) / this.delta), this.nbins);
      grouped[binid].push(row);
    }
    return grouped;
  };

  return Binner;

})();

cal_circular_mean = function(angles) {
  return math.atan2(math.mean(math.sin(angles)), math.mean(math.cos(angles)));
};

between_circular = function(angle, angle1, angle2) {
  return math.mod(angle - angle1, math.pi * 2) < (angle2 - angle1);
};

difference_circular = function(angle1, angle2) {
  return math.atan2(math.sin(angle2 - angle1), math.cos(angle2 - angle1));
};

bbox_overlap = function(rect1, rect2, widthpadding, heightpadding) {
  if (widthpadding == null) {
    widthpadding = 0;
  }
  if (heightpadding == null) {
    heightpadding = 0;
  }
  return !(rect1.x > rect2.x + rect2.width + widthpadding || rect2.x > rect1.x + rect1.width + widthpadding || rect1.y > rect2.y + rect2.height + heightpadding || rect2.y > rect1.y + rect1.height + heightpadding);
};

angles_increase = function(angle1, angle2) {
  return math.mod(angle2 - angle1, math.pi * 2) < math.pi;
};

posneg = function(value) {
  if (value === 0) {
    return 1;
  }
  return value / math.abs(value);
};
