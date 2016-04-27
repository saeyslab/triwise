var Eoi, Gdiffexp, Glabels, Gpin, Loader, activateOverlay, barycoords, binner, changeTagState, colorDirection, dataLoaded, deactivateOverlay, dotplot, enrichmentable, key, loadData, loadExpressionInfo, loadMessages, loader, overlayActive, overlayCaller, roseplot, selectSource, source, sourceTypeNames, sources, tagState, tagStates, tags, update, updateExpressionInfo,
  indexOf = [].indexOf || function(item) { for (var i = 0, l = this.length; i < l; i++) { if (i in this && this[i] === item) return i; } return -1; };

Eoi = 0;

Gdiffexp = 0;

barycoords = 0;

sources = 0;

source = "go";

colorDirection = 0;

enrichmentable = 0;

dotplot = 0;

roseplot = 0;

binner = 0;

Glabels = 0;

tags = 0;

tagState = 0;

Gpin = [];

loader = 0;

key = 0;

dataLoaded = function(data) {
  var C, G, anglebase, ax, container, gsetOverview, gsetid, h, labels, rmax, scale, svg, transformation, updateSources, w;
  console.log("Data downloaded");
  window.data = data;
  Eoi = data.Eoi;
  Gdiffexp = data.Gdiffexp;
  G = Eoi.columns;
  C = Eoi.index;
  labels = C;
  Glabels = data.Glabels;
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
  barycoords = barycoords;
  w = 600;
  h = 400;
  rmax = 4;
  scale = function(x) {
    return x * (h - 0) / (rmax * 2);
  };
  svg = d3.select("#leftplot svg").html("").attr("width", w).attr("height", h);
  ax = svg.append("g").attr("transform", "translate(" + w / 2 + "," + h / 2 + ")");
  dotplot = new Dotplot(ax, w, h, barycoords, rmax, labels, Glabels);
  binner = new Binner(12);
  colorDirection = function(x) {
    return d3.hsl(x % (math.pi * 2) / math.pi * 180, 0.7, 0.65);
  };
  colorDirection = colorDirection;
  svg = d3.select("#rightplot svg").html("").attr("width", w).attr("height", h);
  ax = svg.append("g").attr("transform", "translate(" + w / 2 + "," + h / 2 + ")");
  roseplot = new Roseplot(ax, w, h, barycoords, labels, binner, colorDirection);
  sources = [];
  updateSources = function(gsetOverview) {
    var chooseSource, gsetOverviewNested, prevSourceSelectRequest, sourceButtonsEnter;
    gsetOverviewNested = d3.nest().key(function(d) {
      return d.type;
    }).entries(gsetOverview);
    chooseSource = d3.select("#choose-source").selectAll("div").data(gsetOverviewNested, function(d) {
      return d.key;
    });
    chooseSource.exit().remove();
    chooseSource.enter().append("div").classed("list-group", true).append("label").text(function(d) {
      return sourceTypeNames[d.key];
    });
    sources = chooseSource.selectAll("a").data((function(d) {
      return d.values;
    }), (function(d) {
      return d.source;
    }));
    sources.exit().remove();
    prevSourceSelectRequest = void 0;
    sourceButtonsEnter = sources.enter().append("a").classed("list-group-item", true).text(function(d) {
      return d.name;
    }).on("click", function(d, e) {
      source = d.source;
      return selectSource(source);
    });
    sourceButtonsEnter.append("button").attr("class", "btn btn-link btn-xs pull-right").append("i").attr("class", "fa fa-cog fa-lg");
    sourceButtonsEnter.append("button").attr("class", "btn btn-link btn-xs pull-right").append("i").attr("class", "fa fa-info fa-lg");
    return sourceButtonsEnter.filter(function(d) {
      return d.type === "u";
    }).append("button").attr("class", "btn btn-link btn-xs pull-right").on("click", function(d) {
      $.post("/delete_source", {
        source: d.source
      }, function(d) {
        return updateSources(d.gsetoverview);
      }, "json");
      return d3.event.stopPropagation();
    }).append("span").attr("class", "fa fa-remove fa-lg");
  };
  gsetOverview = data.gsetoverview;
  updateSources(gsetOverview);
  container = $("#upload-gsets");
  $(".fileinput-button input", container).fileupload({
    dataType: "json",
    send: function(e, d) {
      console.log("started");
      return $(".progress-bar", container).attr("class", "progress-bar").css("width", "0%").text("0%");
    },
    done: function(e, d) {
      if (d.result.status === "success") {
        $(".progress-bar", container).attr("class", "progress-bar progress-bar-success").text("Gene sets uploaded");
      } else if (d.result.status === "error") {
        $(".progress-bar", container).attr("class", "progress-bar progress-bar-danger").text("Error: " + d.result.statustext);
      }
      return updateSources(d.result.gsetoverview);
    },
    progress: function(e, data) {
      var progress;
      progress = parseInt(data.loaded / data.total * 100, 10);
      return $(".progress-bar", container).css("width", progress + "%").text(progress + "%");
    }
  });
  source = "go";
  gsetid = "GO:0006355";
  selectSource(source);
  return loader.end("loadData");
};

$(document).ready(function() {
  var colorGoi, colorGoiDiffexp, columnids, container2, directionRenderer, localRenderer, oddsRenderer, pvalRenderer;
  $('[data-toggle="tooltip"]').tooltip();
  container2 = $("#upload-expression");
  $(".fileinput-button input", container2).fileupload({
    dataType: "json",
    send: function(e, d) {
      console.log("started");
      $(".progress-bar", container2).attr("class", "progress-bar").css("width", "0%").text("0%");
      return d3.select("#expression-info").style("display", "none");
    },
    done: function(e, d) {
      if (d.result.status === "success") {
        $(".progress-bar", container2).attr("class", "progress-bar progress-bar-success").text("Expression data uploaded");
        return updateExpressionInfo(d.result);
      } else if ((d.result.status === "error") || (d.result.status === "abort")) {
        return $(".progress-bar", container2).attr("class", "progress-bar progress-bar-danger").text("Error: " + d.result.message);
      }
    },
    progress: function(e, data) {
      var progress;
      progress = parseInt(data.loaded / data.total * 100, 10);
      if (progress < 100) {
        return $(".progress-bar", container2).css("width", progress + "%").text(progress + "%");
      } else {
        return d3.select("#upload-expression .progress-bar").style("width", "100%").text("Processing ").append("i").attr("class", "fa fa-spinner fa-pulse");
      }
    }
  });
  $("button#loadData").on('click', function(e) {
    if (!$(this).hasClass("disabled")) {
      loadData();
      return $("#tabRun").tab('show');
    }
  });
  $("#permalink input").on('click', function() {
    return this.select();
  });
  $("#permalink button").on('click', function() {
    $("#permalink input").select();
    return document.execCommand("copy");
  });
  loader = new Loader(d3.select("#loadMessages"));
  if ((typeof datasetAvailable !== "undefined" && datasetAvailable !== null) && datasetAvailable) {
    loadExpressionInfo();
    loadData();
  }
  columnids = {
    name: 0,
    logqval_enrichment: 1,
    odds_enrichment: 2,
    logqval_local: 3,
    logqvals_local: 4,
    logqval_unidir: 5,
    angle: 6
  };
  directionRenderer = function(angle) {
    var arrowDim, boxdim, x1, x2, y1, y2;
    if (angle === "None") {
      return "";
    }
    boxdim = 15;
    arrowDim = 7;
    x2 = Math.cos(Number(angle)) * arrowDim + boxdim / 2;
    y2 = Math.sin(Number(angle)) * arrowDim + boxdim / 2;
    x1 = Math.cos(Math.PI + Number(angle)) * arrowDim + boxdim / 2;
    y1 = Math.sin(Math.PI + Number(angle)) * arrowDim + boxdim / 2;
    return "<svg width=15 height=15 class=direction><line x1=" + x1 + " y1=" + y1 + " x2=" + x2 + " y2=" + y2 + " stroke=black stroke-width=1 marker-end='url(#arrowHead)'/></svg>";
  };
  localRenderer = function(localqvals) {
    var angle, angle1, angle2, boxdim, color, delta, j, k, l, len, len1, len2, mid, path, paths, point, qval, radius, ref, ref1;
    boxdim = 25;
    mid = boxdim / 2;
    delta = math.pi / 3;
    paths = "";
    radius = boxdim / 2;
    angle = 0;
    for (j = 0, len = localqvals.length; j < len; j++) {
      qval = localqvals[j];
      path = "";
      angle1 = angle - delta / 2;
      angle2 = angle + delta / 2;
      ref = hexagon_path(angle1, angle2);
      for (k = 0, len1 = ref.length; k < len1; k++) {
        point = ref[k];
        path += "L " + (radius * point.x + mid) + " " + (radius * point.y + mid);
      }
      color = dotplot.logpval_scale(qval);
      paths += "<path d='M " + mid + " " + mid + " " + path + " L " + mid + " " + mid + "' fill='" + color + "' />";
      angle += delta;
    }
    ref1 = [0, math.pi / 3 * 2, math.pi + math.pi / 3];
    for (l = 0, len2 = ref1.length; l < len2; l++) {
      angle = ref1[l];
      point = polar2cart({
        "r": radius,
        "angle": angle
      });
      paths += "<path d='M " + mid + " " + mid + " l " + point.x + " " + point.y + "' stroke='#333333' style='opacity:0.3'/>";
    }
    return "<svg width=" + boxdim + " height=" + boxdim + " class=direction>" + paths + "</svg>";
  };
  pvalRenderer = function(logpval) {
    var pval;
    if (logpval === "None") {
      return "<span style='display:none'>1</span>";
    }
    if (logpval === -Infinity) {
      return 0;
    }
    pval = math.pow(10, logpval);
    if (pval < 0.001) {
      return "10" + "<sup>" + math.round(logpval, 1) + "</sup>";
    } else {
      return math.round(pval, 3);
    }
  };
  oddsRenderer = function(odds) {
    if (odds === -Infinity) {
      return "-∞";
    } else if (odds === Infinity) {
      return "+∞";
    } else if (odds === "None") {
      return "<span style='display:none'>1</span>";
    } else {
      return math.round(odds, 2);
    }
  };
  enrichmentable = $('#enrichment').DataTable({
    searching: false,
    info: false,
    lengthChange: false,
    data: [],
    columns: [
      {
        data: 'name',
        className: "overflow"
      }, {
        data: 'logqval_enrichment',
        className: 'divider',
        render: {
          display: pvalRenderer
        }
      }, {
        data: 'odds_enrichment',
        render: {
          display: oddsRenderer
        }
      }, {
        data: 'logqval_local',
        className: 'divider',
        render: {
          display: pvalRenderer
        }
      }, {
        data: 'logqvals_local',
        render: {
          display: localRenderer
        }
      }, {
        data: 'logqval_unidir',
        className: 'divider',
        render: {
          display: pvalRenderer
        }
      }, {
        data: 'angle',
        "render": {
          display: directionRenderer
        },
        "createdCell": function(td, d, rowData, row, col) {
          if (d !== "None") {
            return $(td).css('background-color', colorDirection(Number(d)).toString());
          }
        }
      }
    ],
    order: [[columnids.logqval_unidir, "asc"]],
    columnDefs: [
      {
        "type": "num",
        targets: [columnids.logqval_enrichment, columnids.logqval_unidir, columnids.odds_enrichment]
      }
    ]
  });
  enrichmentable = enrichmentable;
  $("table#enrichment tbody").on('click', 'tr', function(e) {
    var d;
    d = enrichmentable.row(this).data();
    if (!(d != null)) {
      return;
    }
    return update(source, d.gsetid);
  });
  tagState = 0;
  tags = null;
  Gpin = [];
  $("button#togglePinGenes").on("click", function(e) {
    if (tagState === tagStates.hover) {
      return changeTagState(tagStates.pin);
    } else if (tagState === tagStates.pin) {
      return changeTagState(tagStates.hover);
    }
  });
  $("button#download-svg").on("click", function() {
    loader.start("genSVG");
    return $.post("/download_dotplot", {}, function(d) {
      $("body").append("<iframe src='" + d.url + "' style='display: none;' ></iframe>");
      return loader.end("genSVG");
    }, "json");
  });
  colorGoi = function(x) {
    return d3.hsl(x * 255, 0.7, 0.15);
  };
  colorGoiDiffexp = function(x) {
    return d3.hsl(x * 255, 0.9, 0.45);
  };
  return d3.select("#color-schemes").selectAll("button").data(_.range(0, 1 - 0.0001, 0.1)).enter().append("a").attr("class", "btn btn-default").style("background-color", function(d) {
    return colorGoiDiffexp(d);
  }).style("width", "10%").attr("type", "button").on("click", function(d) {
    return dotplot.updateOptions({
      visual: {
        goi: {
          fill: colorGoi(d)
        },
        goidiffexp: {
          fill: colorGoiDiffexp(d)
        }
      }
    });
  });
});

loadData = function() {
  loader.start("loadData");
  return $.post("/load_data", $("form#expression-info").serialize(), dataLoaded, "json");
};

selectSource = function(source) {
  var prevSourceSelectRequest;
  if (typeof prevSourceSelectRequest !== "undefined" && prevSourceSelectRequest !== null) {
    prevSourceSelectRequest.abort();
  }
  loader.start("selectSource");
  sources.classed("active", false);
  sources.filter(function(d) {
    return d.source === source;
  }).classed("active", true);
  return prevSourceSelectRequest = $.ajax({
    dataType: "json",
    url: "/source_enrichment",
    data: JSON.stringify({
      expressionid: 1,
      source: source
    }),
    method: "POST",
    contentType: 'application/json;charset=UTF-8',
    error: (function(_this) {
      return function(d, status, error) {
        return console.log(error);
      };
    })(this),
    success: (function(_this) {
      return function(d) {
        var ax, h, pvalplot, rmax, scale, svg, w;
        console.log("Enrichment downloaded");
        window.enrichment = d;
        enrichmentable.clear();
        enrichmentable.rows.add(d.scores);
        enrichmentable.draw();
        if (enrichmentable.rows().data().length > 0) {
          update(source, enrichmentable.rows().data()[0].gsetid);
        }
        w = 350;
        h = 350;
        rmax = 4;
        scale = function(x) {
          return x * (h - 0) / (rmax * 2);
        };
        svg = d3.select("#pvalplot svg").html("").attr("width", w).attr("height", h);
        ax = svg.append("g").attr("transform", "translate(" + w / 2 + "," + h / 2 + ")");
        return pvalplot = new Pvalplot(ax, w, h, d.scores, rmax, ["to", "be", "added"]);
      };
    })(this),
    complete: function(d) {
      return loader.end("selectSource");
    }
  });
};

loadExpressionInfo = function() {
  return $.post("/expression_info", {}, function(d) {
    return updateExpressionInfo(d);
  }, "json");
};

updateExpressionInfo = function(expression_info) {
  var onSelectedConditionChange, options;
  d3.select("#expression-info").style("display", "inline");
  d3.select("#expression-info #dimensions").text(expression_info.nG + " genes and " + expression_info.nC + " conditions");
  options = d3.select("#select-conditions").html("").selectAll("option").data(expression_info.conditions).enter().append("option").text(function(d) {
    return d.condition + " (" + d.nsamples + ")";
  }).attr("value", function(d) {
    return d.condition;
  });
  onSelectedConditionChange = function() {
    var dropdown, nonSelected, nonSelectedOptions, option, selectedOptions;
    selectedOptions = $('#select-conditions option:selected');
    if (selectedOptions.length >= 3) {
      nonSelectedOptions = $('#select-conditions option').filter(function() {
        return !$(this).is(':selected');
      });
      nonSelected = (function() {
        var j, len, results;
        results = [];
        for (j = 0, len = nonSelectedOptions.length; j < len; j++) {
          option = nonSelectedOptions[j];
          results.push(option.value);
        }
        return results;
      })();
      dropdown = $("#select-conditions").siblings(".btn-group");
      $("input", dropdown).each(function() {
        var ref;
        if (ref = $(this).val(), indexOf.call(nonSelected, ref) >= 0) {
          return $(this).prop('disabled', true).parent('li').addClass('disabled');
        }
      });
      $('button', $('#select-conditions').parent()).addClass("btn-success");
      return $('button#loadData').removeClass("disabled");
    } else {
      dropdown = $("#select-conditions").siblings(".btn-group");
      $("input", dropdown).each(function() {
        return $(this).prop('disabled', false).parent('li').removeClass('disabled');
      });
      $('button', $('#select-conditions').parent()).removeClass("btn-success");
      return $('button#loadData').addClass("disabled");
    }
  };
  $('#select-conditions').multiselect({
    numberDisplayed: 3,
    onChange: onSelectedConditionChange,
    buttonClass: 'form-control',
    buttonWidth: '100%',
    allSelectedText: false,
    maxHeight: 300,
    enableFiltering: true
  });
  $(".multiselect-clear-filter").on("click", function(e) {
    $("#select-conditions").multiselect('deselectAll', false).multiselect('updateButtonText');
    e.stopPropagation();
    return onSelectedConditionChange();
  });
  $('#select-conditions').multiselect("select", expression_info.Coi, true);
  if (key !== expression_info.key) {
    key = expression_info.key;
    $("#permalink input").val("http://127.0.0.1:5000/dataset/" + key);
    return history.pushState({
      key: key
    }, '', "/dataset/" + key);
  }
};

$.fn.dataTable.Api.register('row().show()', function() {
  var new_row_index, page_info, page_to_display, row_position;
  page_info = this.table().page.info();
  new_row_index = this.index();
  row_position = this.table().rows()[0].indexOf(new_row_index);
  page_to_display = Math.floor(row_position / this.table().page.len());
  if (page_to_display >= 0) {
    this.table().page(page_to_display).draw(false);
  }
  return this;
});

d3.selection.prototype.moveToFront = function() {
  return this.each(function() {
    return this.parentNode.appendChild(this);
  });
};

sourceTypeNames = {
  f: "Functions",
  i: "Interactions",
  p: "Pathways",
  r: "Gene regulation",
  u: "User"
};

overlayActive = false;

overlayCaller = 0;

activateOverlay = function(parent, caller) {
  if (!overlayActive) {
    $("#overlay").css({
      top: parent.offset().top,
      width: parent.outerWidth(),
      height: parent.outerHeight()
    });
    $("#overlay").fadeIn(0);
    overlayCaller = caller;
  }
  return overlayActive = true;
};

deactivateOverlay = function(caller) {
  "";
  if (caller === overlayCaller) {
    $("#overlay").fadeOut(0);
    return overlayActive = false;
  }
};

update = function(source, gsetid) {
  loader.start("loadGset");
  return $.post("/get_gset", {
    gsetid: gsetid,
    source: source
  }, function(gsetData) {
    var Goi, bandwidth, barycoords_filtered, deltangle, g, grouped, i, j, k, l, len, len1, len2, ref, row, symbol, tagdata;
    Goi = gsetData;
    enrichmentable.$('tr.selected').removeClass('selected');
    row = enrichmentable.row(function(idx, d) {
      return d.gsetid === gsetid;
    });
    $(row.node()).addClass('selected');
    row.show();
    window.Goi = Goi;
    for (j = 0, len = barycoords.length; j < len; j++) {
      row = barycoords[j];
      if (ref = row.gid, indexOf.call(Goi, ref) >= 0) {
        row.goi = true;
      } else {
        row.goi = false;
      }
    }
    dotplot.updateGoi();
    bandwidth = math.pi / 3;
    deltangle = math.pi / 24;
    $.ajax({
      dataType: "json",
      url: "/locality",
      data: JSON.stringify({
        source: source,
        gsetid: gsetid,
        deltangle: deltangle,
        bandwidth: bandwidth
      }),
      method: "POST",
      contentType: 'application/json;charset=UTF-8',
      success: (function(_this) {
        return function(d) {
          return dotplot.updateRings(d.logpvals, deltangle);
        };
      })(this)
    });
    roseplot.updateGoi();
    barycoords_filtered = [];
    for (k = 0, len1 = barycoords.length; k < len1; k++) {
      row = barycoords[k];
      if (row.goi && row.diffexp) {
        barycoords_filtered.push(row);
      }
    }
    barycoords_filtered.sort(function(d1, d2) {
      return math.mod(d1.angle, math.pi * 2) - math.mod(d2.angle, math.pi * 2);
    });
    grouped = binner.group(barycoords_filtered);
    tagdata = [
      (function() {
        var l, ref1, results;
        results = [];
        for (i = l = 0, ref1 = binner.nbins - 1; 0 <= ref1 ? l <= ref1 : l >= ref1; i = 0 <= ref1 ? ++l : --l) {
          results.push([]);
        }
        return results;
      })()
    ];
    tags = d3.select("#genetags").style("width", "100%").style({
      "text-align": "justify",
      "text-justify": "inter-word"
    }).html("").selectAll("span").data(barycoords_filtered);
    tags.enter().append("span").classed("genetag", true).style({
      "padding-top": "2px",
      "padding-bottom": "2px",
      "padding-left": "3px",
      "padding-right": "3px"
    }).style("display", "inline-block").style("text-color", "white").style("border-radius", "4px");
    tags.exit().remove();
    tags.style("background-color", function(d) {
      return colorDirection(d.angle);
    });
    tags.append("a").text(function(d) {
      return Glabels[d.gid];
    }).on('mouseover', function(d) {
      return dotplot.updateHover(d.gid);
    }).on('mouseout', function(d) {
      return dotplot.updateHover(void 0);
    }).on('click', function(d) {
      var ref1;
      if (tagState === tagStates.hover) {
        return console.log("not implemented!");
      } else if (tagState === tagStates.pin) {
        if (ref1 = d.gid, indexOf.call(Gpin, ref1) < 0) {
          Gpin.push(d.gid);
          dotplot.updateGpin([d.gid]);
        } else {
          Gpin.splice(Gpin.indexOf(d.gid), 1);
          dotplot.updateGpin([], [d.gid]);
        }
        return changeTagState(tagStates.pin);
      }
    }).attr("role", "button").attr("tabindex", function(d, i) {
      return i;
    });
    Gpin = (function() {
      var l, len2, ref1, results;
      ref1 = ["Pole", "Pola2", "Pola1", "Pold2"];
      results = [];
      for (l = 0, len2 = ref1.length; l < len2; l++) {
        symbol = ref1[l];
        results.push(data.Glabels.indexOf(symbol));
      }
      return results;
    })();
    Gpin = [];
    for (l = 0, len2 = Goi.length; l < len2; l++) {
      g = Goi[l];
      if (Gdiffexp[g]) {
        Gpin.push(g);
      }
    }
    Gpin = _.pluck(_.sortBy(barycoords_filtered, "r"), "gid");
    if (Gpin.length > 20) {
      Gpin = Gpin.splice(Gpin.length - 10, 10);
    }
    dotplot.initGpin(Gpin);
    changeTagState(tagStates.hover);
    return loader.end("loadGset");
  }, "json");
};

loadMessages = {
  loadData: "Loading data",
  selectSource: "Loading gene set enrichment",
  loadGset: "Loading gene set",
  genSVG: "Generating SVG"
};

Loader = (function() {
  function Loader(element) {
    this.element = element;
    this.queue = [];
  }

  Loader.prototype.end = function(process) {
    this.queue.splice(this.queue.indexOf(process), 1);
    return this.update();
  };

  Loader.prototype.start = function(process) {
    this.queue.push(process);
    return this.update();
  };

  Loader.prototype.update = function() {
    var j, len, messages, process, ref;
    if (this.queue.length > 0) {
      this.element.transition(1).style("opacity", 1);
      messages = [];
      ref = this.queue;
      for (j = 0, len = ref.length; j < len; j++) {
        process = ref[j];
        if (process in loadMessages) {
          messages.push(loadMessages[process]);
        } else {
          messages.push(process);
        }
      }
      return this.element.select("span#messages").text(messages.join(", "));
    } else {
      return this.element.transition(1).style("opacity", 0);
    }
  };

  return Loader;

})();

tagStates = {
  hover: 0,
  pin: 1
};

changeTagState = function(newstate) {
  if (tagState === tagStates.hover || !(tagState != null)) {
    "";
  } else if (tagState === tagStates.pin) {
    tags.classed("pinned", false);
    d3.select("button#togglePinGenes").classed("active", false);
  }
  tagState = newstate;
  if (tagState === tagStates.hover) {
    return "";
  } else if (tagState === tagStates.pin) {
    tags.filter(function(d) {
      var ref;
      return ref = d.gid, indexOf.call(Gpin, ref) >= 0;
    }).classed("pinned", true);
    return d3.select("button#togglePinGenes").classed("active", true);
  }
};
