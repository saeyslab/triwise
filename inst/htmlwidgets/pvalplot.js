var h, w;

w = 0;

h = 0;

HTMLWidgets.widget({
  name: 'pvalplot',
  type: 'output',
  initialize: function(el, width, height) {
    w = width;
    return h = height;
  },
  renderValue: function(el, data, instance) {
    var ax, gsetlabels, labels, pvalplot, rmax, scores;
    console.log(data);
    window.data = data;
    scores = data.scores;
    labels = data.labels;
    gsetlabels = data.gsetlabels;
    rmax = 5;
    ax = d3.select("div.pvalplot").append("svg").attr("width", w).attr("height", h).append("g").attr("transform", "translate(" + w / 2 + "," + h / 2 + ")");
    return pvalplot = new Pvalplot(ax, w, h, scores, rmax, labels, gsetlabels);
  },
  resize: function(el, width, height, instance) {}
});
