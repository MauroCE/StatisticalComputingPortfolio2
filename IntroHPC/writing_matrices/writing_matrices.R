m1 = matrix(1:9, 3, 3)
m2 = matrix(10:18, 3, 3)
write(c(m1), file="matrix_file.txt", ncol = 9)
write(c(m2), file="matrix_file.txt", ncols = 9, append=TRUE)