import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import pytest
import numpy as np
import symplectic as sy

## Example inputs
# Can be used in multiple test functions, to avoid redundancy.
lower_tri_matrix = np.array([
                            [1, 0, 0, 0],
                            [1, 1, 0, 0],
                            [1, 1, 1, 0],
                            [1, 1, 1, 1]
                            ])

matrix_4_x_3 = np.array([[1, 1, 0],
                         [1, 0, 1],
                         [0, 1, 1],
                         [1, 1, 1]
                         ])

matrix_4_x_5 = np.array([[1, 1, 0, 1, 1],
                         [0, 1, 1, 0, 1],
                         [1, 1, 1, 1, 0],
                         [0, 0, 1, 1, 1]
                         ])

symplectic_matrix = np.array([[1, 0, 1, 0],
                              [0, 1, 0, 1],
                              [0, 0, 1, 0],
                              [0, 0, 0, 1]
                              ])

mat = np.array([[1, 0, 1, 1],
                [0, 1, 0, 0],
                [1, 1, 1, 0],
                [0, 0, 0, 1]
                ])






## Tests for is_symplectic()
    #
@pytest.mark.parametrize("input, expected_result", [
    (symplectic_matrix, True),
    (mat, False),
    (matrix_4_x_3, False),
    ])
def test_is_symplectic(input, expected_result):
        assert sy.is_symplectic(input) == expected_result, "Symplectic check of matrix is not correct."
        
# def test_is_symplectic_product_coord(input, expected_result):
#     assert sy.is_symplectic_product_coord(input) == expected_result, "Symplectic check of matrix is not correct."


## Tests for ref()
#
@pytest.mark.parametrize("input, expected_result", [
    (mat, np.array([[1, 0, 1, 1],
                    [0, 1, 0, 0],
                    [0, 0, 0, 1],
                    [0, 0, 0, 0]])),
    (lower_tri_matrix, np.array([[1, 0, 0, 0],
                                  [0, 1, 0, 0],
                                  [0, 0, 1, 0],
                                  [0, 0, 0, 1]])),
    (matrix_4_x_3, np.array([[1, 1, 0],
                              [0, 1, 1],
                              [0, 0, 1],
                              [0, 0, 0]])),
    (matrix_4_x_5, np.array([[1, 1, 0, 1, 1],
                              [0, 1, 1, 0, 1],
                              [0, 0, 1, 0, 1],
                              [0, 0, 0, 1, 0]])),
    ])
def test_ref(input, expected_result):
    input_ref = sy.ref(input)
    assert np.array_equal(input_ref, expected_result), "Row echolon form of matrix is not correct."


## Tests for rref()
#
@pytest.mark.parametrize("Mat, expected_result", [
    (mat, np.array([[1, 0, 1, 0],
                    [0, 1, 0, 0],
                    [0, 0, 0, 1],
                    [0, 0, 0, 0]])),
    (lower_tri_matrix, np.array([[1, 0, 0, 0],
                                  [0, 1, 0, 0],
                                  [0, 0, 1, 0],
                                  [0, 0, 0, 1]])),
    (matrix_4_x_3, np.array([[1, 0, 0],
                              [0, 1, 0],
                              [0, 0, 1],
                              [0, 0, 0]])),
    (matrix_4_x_5, np.array([[1, 0, 0, 0, 1],
                              [0, 1, 0, 0, 0],
                              [0, 0, 1, 0, 1],
                              [0, 0, 0, 1, 0]])),
    ])
def test_rref(Mat, expected_result):
        input_rref, _ = sy.rref(Mat)
        assert np.array_equal(input_rref, expected_result), "Reduced row echolon form of matrix is not correct."


mat1 = np.array([[1, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0],
                         [0, 0, 1, 0, 0],
                         [0, 0, 0, 1, 1],
                         [0, 0, 0, 1, 1]]
                        )
mat2 = np.array([[1, 0, 0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 1, 0, 0],
                         [0, 0, 0, 0, 0, 0, 1, 0],
                         [0, 0, 0, 0, 0, 0, 0, 1]]
                        )
## test for direct_sum()
#
@pytest.mark.parametrize("args, expected_result", [
    ((np.eye(3), np.ones([2,2])), mat1),
    ((sy.identity(2), sy.zero_matrix(2, 3), sy.identity(3)), mat2)
    ])
def test_direct_sum(args, expected_result):
    result = sy.direct_sum(*args)
    print(args)
    print(result)
    print(expected_result)
    assert np.array_equal(result, expected_result), "Direct sum of matrices is not correct"
    
## test for symplectic_form()
#
@pytest.mark.parametrize("v, w, expected_result", [
    # Example 1: standard basis vectors
    (np.array([1, 0, 0, 0], dtype=np.uint8), np.array([1, 0, 1, 1], dtype=np.uint8), 1),
    
    # # Example 2: commuting Pauli operators (e.g., ZI and IZ)
    (np.array([1, 0, 0, 0], dtype=np.uint8), np.array([0, 1, 0, 1], dtype=np.uint8), 0),
    
    # # Example 3: anti-commuting pair (e.g., XZ and ZX on 2 qubits)
    (np.array([1, 1, 1, 0], dtype=np.uint8), np.array([1, 0, 0, 1], dtype=np.uint8), 0),
])
def test_symplectic_form(v, w, expected_result):
    result = sy.symplectic_form(v, w)
    assert result == expected_result, "Symplectic form of binary vectors is wrong"


