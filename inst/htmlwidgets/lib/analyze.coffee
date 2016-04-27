## ANALYZE ##

transform_barycentric = (data, transformation) ->
	genecoords = math.transpose(math.multiply(transformation, data))

	return(genecoords)
