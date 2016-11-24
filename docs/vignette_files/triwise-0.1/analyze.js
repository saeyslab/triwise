var transform_barycentric;

transform_barycentric = function(data, transformation) {
  var genecoords;
  genecoords = math.transpose(math.multiply(transformation, data));
  return genecoords;
};
