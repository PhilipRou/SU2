include("test_head.jl")



function ran_U2_test()
    @assert is_U2(ran_U2(rand())) "Something went wrong while testing ran_U2"
    return true
end
# It was also tested that for ϵ → 0 no -id can appear anymore

ran_U2_test()

function mult_U2_test()
    coeffs_1 = ran_U2(rand())
    coeffs_2 = ran_U2(rand())
    mat_1 = coeffs2grp(coeffs_1)
    mat_2 = coeffs2grp(coeffs_2)
    prod_coeffs = coeffs_1 * coeffs_2
    prod_mat = coeffs2grp(prod_coeffs)
    @assert isapprox(prod_mat, mat_1*mat_2) "Something went wrong while testing mult_U2"
    return true
end

mult_U2_test()

function tr_test()
    coeffs = ran_U2(rand())
    mat = coeffs2grp(coeffs)
    @assert isapprox(tr(coeffs), mat[1,1] + mat[2,2])
    return true
end

tr_test()

function det_test()
    coeffs = ran_U2(rand())
    mat = coeffs2grp(coeffs)
    @assert isapprox(det(mat), det(coeffs))
    return true
end

det_test()

function proj2man_test()
    coeffs = rand() * ran_U2(rand())
    @assert is_U2(proj2man(coeffs))
    return true
end

proj2man_test()

function log_U2_test()
    bla = ran_U2(rand())
    @assert isapprox(proj2man(log_U2(bla)), grp2coeffs_U2(proj2man_mat_U2(log(coeffs2grp(bla)))))
    return true
end

log_U2_test()   

function exp_u2_test()
    bla = ran_U2(rand())
    @assert isapprox(bla, exp_u2(log_U2(bla)))
    return true
end

exp_u2_test()

function coeffs2grp_test()
    coeffs = ran_U2(rand())
    @assert is_U2(coeffs2grp(coeffs)) "Something went wrong while testing coeffs2grp"
    return true
end

coeffs2grp_test()

function grp2coeffs_U2_test()
    coeffs = ran_U2(rand())
    mat = coeffs2grp(coeffs)
    new_coeffs = grp2coeffs_U2(mat)
    # @assert isapprox(coeffs.a, new_coeffs.a) * isapprox(coeffs.b, new_coeffs.b) * isapprox(coeffs.c, new_coeffs.c) * isapprox(coeffs.d, new_coeffs.d)
    @assert is_U2(mat)
    # println(new_coeffs)
    @assert isapprox(coeffs, new_coeffs)
    return true
end

grp2coeffs_U2_test()

function adjoint_test()
    coeffs = ran_U2(rand())
    adj_coeffs = adjoint(coeffs)
    mat = coeffs2grp(coeffs)
    adj_mat = coeffs2grp(adj_coeffs)
    @assert isapprox(adj_mat, adjoint(mat)) "Something went wrong while testing adjoint"
    return true
end

adjoint_test()

function temp_gauge_U2_test()
    test_field = gaugefield_U2(N_x, N_t, true)
    coords = [[rand(1:N_x), rand(1:N_t)], [1,1], [1,N_t], [N_x,1], [N_x,N_t]]
    for coord in coords
        x, t = coord
        old_plaq = tr(plaq(test_field,x,t))
        new_plaq = tr(plaq(temp_gauge_U2(test_field),x,t))
        @assert isapprox(old_plaq, new_plaq) "temp_gauge_U2() not yet correct at x = $x, t = $t"
    end
    return true
end

temp_gauge_U2_test()





####    Hexagonal Shenanigans    ####





function temp_gauge_hex_U2_test()
    test_field = hexfield_U2(N_x, N_t, true)
    coords = [rand(half_chess_coords(N_x,N_t)), [1,1], [2,N_t], [N_x,2], [N_x,N_t], [2*rand(1:div(N_x,2)), 2*rand(1:div(N_t,2))]]
    for coord in coords
        x, t = coord
        old_plaq = tr(hexplaq(test_field,x,t))
        new_plaq = tr(hexplaq(temp_gauge_hex_U2(test_field),x,t))
        @assert isapprox(old_plaq, new_plaq) "temp_gauge_test() not yet correct at x = $x, t = $t"
    end
    return true
end

temp_gauge_hex_U2_test()
