import picos as pic
import numpy as np
from sys import argv


ROUNDING_FACTOR = 6


def parse_input_matrix_form(filepath):
    vectors = []
    rows = 0
    cols = 0
    with open(filepath, "r") as f:
        for line in f:
            rows += 1
            temp = [float(_) for _ in line.split()]
            cols = len(temp)
            vectors.append(temp)
    return vectors, rows, cols


def init_square_matrix(size):
    return [[float(0) for x in range(size)] for y in range(size)]


def solve(param, path):
    adj_matrix, size, size = parse_input_matrix_form(path)

    # Formulate the SDP
    prob = pic.Problem()
    A = pic.new_param('A', adj_matrix)
    k = pic.new_param('k', param)
    X = prob.add_variable('X', (size, size), 'symmetric')

    # Constraints
    prob.add_constraint(pic.sum([X[i, i] for i in range(size)], 'i', '0..size-1') == k)

    for i in range(size):
        prob.add_constraint(pic.sum([X[i, j]/2.0 + X[j, i]/2.0 for j in range(size)], 'j', '0..size-1') <= k * X[i, i])

    for i in range(size):
        for j in range(size):
            prob.add_constraint(X[i, j]/2.0 + X[j, i]/2.0 <= X[i, i])

    for i in range(size):
        for j in range(size):
            prob.add_constraint(X[i, j]/2.0 + X[j, i]/2.0 >= 0.0)

    for i in range(size):
        prob.add_constraint(X[i, i] <= 1.0)

    prob.add_constraint(X >> 0)

    # Objective
    prob.set_objective('max', (A | X) / 2.0)
    prob.solve(solver='cvxopt', verbose=0)

    # Output
    print prob
    print 'Type:   ' + prob.type
    print '\nStatus: ' + prob.status

    print "\n<<Input Variables>>\n"
    print 'A :'
    print np.asarray(adj_matrix)

    print '\nK :\n' + str(param)

    print "\n<<Primal Variables>>"
    X_opt = X.value
    # --------------------#
    #  objective value   #
    # --------------------#
    print '\nThe optimal value of this problem (Objective) <A, X*> is :'
    print prob.obj_value() 
    # --------------------#
    #  optimal variable  #
    # --------------------#
    print '\nX* :' 
    print X_opt

    print "<<Dual Variables>>"

    constraint_id = 0

    t = prob.get_constraint(constraint_id).dual
    print '\nt : \n' + str(t)

    constraint_id += 1

    print ""
    Y_1 = []
    for i in range(1, size+1):
        Y_1.append(prob.get_constraint(i).dual)
        constraint_id += 1
    print 'Y_1 :'
    print [round(x, ROUNDING_FACTOR) for x in Y_1]

    Z_1 = init_square_matrix(size)
    print ""
    for i in range(0, size):
        for j in range(0, size):
            Z_1[i][j] = prob.get_constraint(constraint_id).dual
            constraint_id += 1

    print "Z_1 : "
    print np.asarray(Z_1).round(ROUNDING_FACTOR)

    Z_2 = init_square_matrix(size)
    print ""
    for i in range(0, size):
        for j in range(0, size):
            Z_2[i][j] = prob.get_constraint(constraint_id).dual
            constraint_id += 1

    print "Z_2 : "
    print np.asarray(Z_2).round(ROUNDING_FACTOR)

    
    print ""
    Y_2 = []
    for i in range(1, size + 1):
        Y_2.append(prob.get_constraint(i).dual)
        constraint_id += 1
    print 'Y_2 :'
    print [round(x, ROUNDING_FACTOR) for x in Y_2]

    print "\nU-A :"
    dual_X_psd = prob.get_constraint(constraint_id).dual
    print dual_X_psd

    U_A = init_square_matrix(size)
    for i in range(0, size):
        for j in range(0, size):
            U_A[i][j] = dual_X_psd[i, j]

    print "Eigen-values of U-A"
    print sorted(np.linalg.eigvals(np.array(U_A).round(ROUNDING_FACTOR)))

    print "\nU:"
    print np.array(np.array(U_A) + np.array(adj_matrix)).round(ROUNDING_FACTOR)


def main():
    script, path, param = argv
    solve(int(param), path)


if __name__ == "__main__": main()