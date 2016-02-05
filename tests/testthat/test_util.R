test_that("hexagonPolar gives the correct radius given an angle", {
  expect_equal(hexagonPolar(0, 1), 1)
  expect_equal(hexagonPolar(pi/3, 1), 1)
  expect_equal(hexagonPolar(pi/3+pi, 1), 1)
  expect_true(hexagonPolar(pi/6, 1) < 1)
})

test_that("betweenCircular returns true if the given angle is between angle1 and angle2", {
  expect_true(betweenCircular(1, 0.9, 1.2))
  expect_true(betweenCircular(0, -0.1, 0.1))
  expect_false(betweenCircular(1, 1.2, 0.9))
})

test_that("jaccard works as expected", {
  expect_equal(jaccard(c("1", "2", "3"), c("1", "2", "3", "4")), 0.75)
})
