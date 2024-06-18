# Test check_cmd function
test_that("check_cmd returns 0 if command is available", {
    cmd <- "ls"

    result <- check_cmd(cmd)

    expect_equal(result, 0)
})

# Test check_cmd function returns error if command is not available
test_that("check_cmd returns error if command is not available", {
    cmd <- "not_a_command_arble_garble"

    expect_error(check_cmd(cmd))
})
