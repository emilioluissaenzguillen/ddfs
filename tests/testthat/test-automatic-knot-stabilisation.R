test_that("automatic trimming is clipped to its tuning domain", {
  expect_equal(clip_trimming_alpha(-0.2), 0)
  expect_equal(clip_trimming_alpha(0.05), 0.05)
  expect_equal(clip_trimming_alpha(0.10), 0.10)
  expect_equal(clip_trimming_alpha(0.35), 0.10)
})
