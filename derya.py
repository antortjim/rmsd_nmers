def rmsd(x_value, y_value):
    import numpy
    ##take all the average for  each axes for all amoino acids
    mass_x_value = numpy.transpose(numpy.matrix([[numpy.average(x_value[0, :]), numpy.average(x_value[1, :]), numpy.average(x_value[2, :])]]))
    mass_y_value = numpy.transpose(numpy.matrix([[numpy.average(y_value[0, :]), numpy.average(y_value[1, :]), numpy.average(y_value[2, :])]]))
    mass_center_subtract = numpy.subtract(mass_x_value, mass_y_value)
    ##rotate the y value
    centered_y_value = y_value + mass_center_subtract

    ##find the transpose x and y value
    y_value_transpose = numpy.transpose(centered_y_value)
    x_value_transpose = numpy.transpose(x_value)

    ## R = Y*Xt
    R_value_2 = numpy.dot(centered_y_value, x_value_transpose)

    ##SVD Decompostion

    v, s, W_t = numpy.linalg.svd(R_value_2, full_matrices=True)

    w = numpy.transpose(W_t)

    V_t = numpy.transpose(v)

    u = numpy.dot(w, V_t)

    z_diag = numpy.array([1, 1, -1])

    z = numpy.diag(z_diag)

    ##E zero value
    total_sum = 0
    for i in range(0, (len(x_value[1,:]))):
        val_x = numpy.square(numpy.absolute(x_value[:, i]))
        print val_x
        val_y = numpy.square(numpy.absolute(y_value[:, i]))
        print val_y
        total_value = val_x + val_y
        sum(total_value)
        total_sum = total_sum + sum(total_value)

    ##rotate y if the det(U) = -1
    if int(numpy.linalg.det(u)) == -1:
        rmsd = numpy.sqrt(numpy.multiply(1/(len(x_value[1,:])), (total_sum - numpy.multiply(2, (s[0] + s[1] - s[2])))))
    else:
        rmsd = numpy.sqrt(numpy.multiply(1/(len(x_value[1,:]), (total_sum - numpy.multiply(2, (s[0] + s[1] + s[2])))))) 

    return rmsd
