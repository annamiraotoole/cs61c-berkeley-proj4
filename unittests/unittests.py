from utils import *
from unittest import TestCase
import dumbpy

"""
For each operation, you should write tests to test  on matrices of different sizes.
Hint: use dp_mc_matrix to generate dumbpy and numc matrices with the same data and use
      cmp_dp_nc_matrix to compare the results
"""
class TestAdd(TestCase):
    def test_small_add(self):
        # TODO: YOUR CODE HERE
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(2, 2, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(2, 2, seed=1)
        
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

        # nc_result = nc_mat1 + nc_mat2
        # dp_result = dp_mat1 + dp_mat2
        # self.assertEqual(dp_result.get(0, 0), nc_result.get(0,0))
        # self.assertEqual(dp_result.get(0, 1), nc_result.get(0,1))
        # self.assertEqual(dp_result.get(1, 0), nc_result.get(1,0))
        # self.assertEqual(dp_result.get(1, 1), nc_result.get(1,1))

    def test_medium_add(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(200, 200, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(200, 200, seed=1)
        
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

 
    def test_large_add(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(2000, 2000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(2000, 2000, seed=1)
        
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
    
    def test_large_add_2(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(2015, 2002, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(2015, 2002, seed=1)

        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_verylarge_add(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(6000, 6000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(6000, 6000, seed=1)
        
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
    
    def test_gigantic_add(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(10000, 10000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(10000, 10000, seed=1)

        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
 

class TestSub(TestCase):
    num_seeds = 1
    def test_small_sub(self):
        total = 0
        for i in range(self.num_seeds):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix(2, 2, seed=0)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix(2, 2, seed=1)

            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "sub")
            self.assertTrue(is_correct)
            print_speedup(speed_up)
            total += speed_up
        print()
        print("!!!! average speed up is", total/self.num_seeds, "!!!!")
       #  nc_result = nc_mat1 - nc_mat2
       #  dp_result = dp_mat1 - dp_mat2
       #  self.assertEqual(dp_result.get(0, 0), nc_result.get(0,0))
       #  self.assertEqual(dp_result.get(0, 1), nc_result.get(0,1))
       #  self.assertEqual(dp_result.get(1, 0), nc_result.get(1,0))
       #  self.assertEqual(dp_result.get(1, 1), nc_result.get(1,1))

    def test_medium_sub(self):
        total = 0
        for i in range(self.num_seeds):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix(20, 200, seed=0)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix(20, 200, seed=1)

            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "sub")
            self.assertTrue(is_correct)
            print_speedup(speed_up)
            total += speed_up
        print()
        print("!!!! average speed up is", total/self.num_seeds, "!!!!")

    def test_medium_sub2(self):
        total = 0
        for i in range(self.num_seeds):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix(21, 203, seed=0)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix(21, 203, seed=1)

            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "sub")
            self.assertTrue(is_correct)
            print_speedup(speed_up)
            total += speed_up
        print()
        print("!!!! average speed up is", total/self.num_seeds, "!!!!")


    def test_large_sub(self):
        total = 0
        for i in range(self.num_seeds):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix(200, 200, seed=0)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix(200, 200, seed=1)

            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "sub")
            self.assertTrue(is_correct)
            print_speedup(speed_up)
            total += speed_up
        print()
        print("!!!! average speed up is", total/self.num_seeds, "!!!!")

    def test_very_large_sub(self):
        total = 0
        for i in range(self.num_seeds):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix(2000, 1501, seed=0)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix(2000, 1501, seed=1)

            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "sub")
            self.assertTrue(is_correct)
            print_speedup(speed_up)
            total += speed_up
        print()
        print("!!!! average speed up is", total/self.num_seeds, "!!!!")

class TestAbs(TestCase):
    def test_small_abs(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(2, 2, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "abs")
        self.assertTrue(is_correct)
        print()
        print_speedup(speed_up)

    def test_medium_abs(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(200, 20, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "abs")
        self.assertTrue(is_correct)
        print()
        print_speedup(speed_up)

    def test_medium_abs2(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(202, 21, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "abs")
        self.assertTrue(is_correct)
        print()
        print_speedup(speed_up)

    def test_large_abs(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(200, 200, seed=0)
        is_correct, speed_up = compute([dp_mat], [nc_mat], "abs")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

class TestNeg(TestCase):
    num_seeds = 1
    def test_small_neg(self):
        total = 0
        for i in range(self.num_seeds):
            dp_mat, nc_mat = rand_dp_nc_matrix(2, 2, seed=0)
            is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
            self.assertTrue(is_correct)
            print_speedup(speed_up)
            total += speed_up
        print()
        print("!!!! average speed up is", total/self.num_seeds, "!!!!")
    
    def test_medium_neg(self):
        total = 0
        for i in range(self.num_seeds):
            dp_mat, nc_mat = rand_dp_nc_matrix(200, 20, seed=0)
            is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
            self.assertTrue(is_correct)
            print_speedup(speed_up)
            total += speed_up
        print()
        print("!!!! average speed up is", total/self.num_seeds, "!!!!")
    
    def test_medium_neg2(self):
        total = 0
        for i in range(self.num_seeds):
            dp_mat, nc_mat = rand_dp_nc_matrix(203, 20, seed=0)
            is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
            self.assertTrue(is_correct)
            print_speedup(speed_up)
            total += speed_up
        print()
        print("!!!! average speed up is", total/self.num_seeds, "!!!!")

    def test_large_neg(self):
        total = 0
        for i in range(self.num_seeds):
            dp_mat, nc_mat = rand_dp_nc_matrix(200, 200, seed=0)
            is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
            self.assertTrue(is_correct)
            print_speedup(speed_up)
            total += speed_up
        print()
        print("!!!! average speed up is", total/self.num_seeds, "!!!!")

    def test_very_large_neg(self):
        total = 0
        for i in range(self.num_seeds):
            dp_mat, nc_mat = rand_dp_nc_matrix(2000, 1501, seed=0)
            is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
            self.assertTrue(is_correct)
            print_speedup(speed_up)
            total += speed_up
        print()
        print("!!!! average speed up is", total/self.num_seeds, "!!!!")

class TestMul(TestCase):
    def test_small_mul(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(6, 2, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(2, 9, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_mul(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(20, 200, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(200, 30, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

        dp_mat1, nc_mat1 = rand_dp_nc_matrix(200, 20, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(20, 300, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_mul2(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(21, 203, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(203, 32, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_mul(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(300, 200, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(200, 400, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

        dp_mat1, nc_mat1 = rand_dp_nc_matrix(2000, 20, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(20, 300, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_verylarge_mul(self):
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(3000, 2000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(2000, 4000, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

        dp_mat1, nc_mat1 = rand_dp_nc_matrix(2000, 200, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(200, 3000, seed=1)
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

class TestPow(TestCase):

    def test_zero_pow(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(1000, 1000, seed=0)
        is_correct, speed_up = compute([dp_mat, 0], [nc_mat, 0], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_one_pow(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(1000, 1000, seed=0)
        is_correct, speed_up = compute([dp_mat, 1], [nc_mat, 1], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_small_pow(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(20, 20, seed=0)
        is_correct, speed_up = compute([dp_mat, 5], [nc_mat, 5], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_pow(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(500, 500, seed=0)
        is_correct, speed_up = compute([dp_mat, 3], [nc_mat, 3], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_large_pow(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(1000, 1000, seed=0)
        is_correct, speed_up = compute([dp_mat, 3], [nc_mat, 3], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_smallmat_largepow_10_20(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(10, 10, seed=0)
        is_correct, speed_up = compute([dp_mat, 20], [nc_mat, 20], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_smallmat_largepow_3_50(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 3, seed=0)
        is_correct, speed_up = compute([dp_mat, 50], [nc_mat, 50], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_smallmat_largepow_3_50(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(3, 3, seed=0)
        is_correct, speed_up = compute([dp_mat, 100], [nc_mat, 100], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)


class TestGet(TestCase):
    def test_get(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(2, 2, seed=0)
        # rand_row = np.random.randint(dp_mat.shape[0])
        # rand_col = np.random.randint(dp_mat.shape[1])
        # self.assertEqual(round(dp_mat[rand_row][rand_col], decimal_places),
        #     round(nc_mat[rand_row][rand_col], decimal_places))
        self.assertEqual(dp_mat.get(0, 0), nc_mat.get(0,0))
        self.assertEqual(dp_mat.get(0, 1), nc_mat.get(0,1))

class TestSet(TestCase):
    def test_set(self):
        # TODO: YOUR CODE HERE
        dp_mat, nc_mat = rand_dp_nc_matrix(2, 2, seed=0)
        # rand_row = np.random.randint(dp_mat.shape[0])
        # rand_col = np.random.randint(dp_mat.shape[1])
        # self.assertEquals(round(dp_mat[rand_row][rand_col], decimal_places),
        #     round(nc_mat[rand_row][rand_col], decimal_places))
        dp_mat.set(0, 0, 100.1)
        nc_mat.set(0, 0, 100.1)
        self.assertEqual(dp_mat.get(0, 0), nc_mat.get(0, 0))
        self.assertEqual(nc_mat.get(0, 0), 100.1)

class TestShape(TestCase):
    def test_shape(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(2, 2, seed=0)
        self.assertTrue(dp_mat.shape == nc_mat.shape)

class TestSubscript(TestCase):
    def test_int_int(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(2, 2, seed=0)
        self.assertEqual(dp_mat[0, 0], nc_mat[0, 0])
        self.assertEqual(dp_mat[0, 1], nc_mat[0, 1])
        dp_mat, nc_mat = rand_dp_nc_matrix(10, 10, seed=0)
        self.assertEqual(dp_mat[4, 9], nc_mat[4, 9])
        self.assertEqual(dp_mat[0, 9], nc_mat[0, 9])
        self.assertEqual(dp_mat[0, 6], nc_mat[0, 6])
        self.assertEqual(dp_mat[3, 7], nc_mat[3, 7])

    def test_1D_int(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(10, 1, seed=0)
        self.assertEqual(dp_mat[0], nc_mat[0])
        self.assertEqual(dp_mat[9], nc_mat[9])
        self.assertEqual(dp_mat[5], nc_mat[5])

        dp_mat, nc_mat = rand_dp_nc_matrix(1, 10, seed=0)
        self.assertEqual(dp_mat[0], nc_mat[0])
        self.assertEqual(dp_mat[9], nc_mat[9])
        self.assertEqual(dp_mat[5], nc_mat[5])

    def test_2D_slice(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(10, 11, seed=0)
        self.assertEqual(dp_mat[2:8], nc_mat[2:8])
        self.assertEqual(dp_mat[0:1], nc_mat[0:1])
        self.assertEqual(dp_mat[5:10], nc_mat[5:10])
        self.assertEqual(dp_mat[0:10], nc_mat[0:10])

        self.assertEqual(dp_mat[0:-1], nc_mat[0:-1])
        self.assertEqual(dp_mat[-8:-5], nc_mat[-8:-5])
        self.assertEqual(dp_mat[5:-1], nc_mat[5:-1])
        self.assertEqual(dp_mat[-2:-1], nc_mat[-2:-1])

    def test_1D_slice(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(1, 10, seed=0)
        self.assertEqual(dp_mat[2:8], nc_mat[2:8])
        self.assertEqual(dp_mat[0:1], nc_mat[0:1])
        self.assertEqual(dp_mat[5:10], nc_mat[5:10])
        self.assertEqual(dp_mat[0:10], nc_mat[0:10])
    
        self.assertEqual(dp_mat[0:-1], nc_mat[0:-1])
        self.assertEqual(dp_mat[-8:-5], nc_mat[-8:-5])
        self.assertEqual(dp_mat[5:-1], nc_mat[5:-1])
        self.assertEqual(dp_mat[-2:-1], nc_mat[-2:-1])

        dp_mat, nc_mat = rand_dp_nc_matrix(10, 1, seed=0)
        self.assertEqual(dp_mat[2:8], nc_mat[2:8])
        self.assertEqual(dp_mat[0:1], nc_mat[0:1])
        self.assertEqual(dp_mat[5:10], nc_mat[5:10])
        self.assertEqual(dp_mat[0:10], nc_mat[0:10])

        self.assertEqual(dp_mat[0:-1], nc_mat[0:-1])
        self.assertEqual(dp_mat[-8:-5], nc_mat[-8:-5])
        self.assertEqual(dp_mat[5:-1], nc_mat[5:-1])
        self.assertEqual(dp_mat[-2:-1], nc_mat[-2:-1])

    def test_2D_int(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(10, 11, seed=0)
        self.assertEqual(dp_mat[0], nc_mat[0])
        self.assertEqual(dp_mat[9], nc_mat[9])
        self.assertEqual(dp_mat[5], nc_mat[5])
    
    def test_int_slice(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(10, 11, seed=0)
        self.assertEqual(dp_mat[0, 0:11], nc_mat[0, 0:11])
        self.assertEqual(dp_mat[9, 0:11], nc_mat[9, 0:11])
        self.assertEqual(dp_mat[0, 4:8], nc_mat[0, 4:8])
        self.assertEqual(dp_mat[9, 2:7], nc_mat[9, 2:7])
        self.assertEqual(dp_mat[5, 1:8], nc_mat[5, 1:8])
        self.assertEqual(dp_mat[5, 0:1], nc_mat[5, 0:1])
        self.assertEqual(dp_mat[0, 10:11], nc_mat[0, 10:11])

        self.assertEqual(dp_mat[0, 0:-1], nc_mat[0, 0:-1])
        self.assertEqual(dp_mat[5, -6:-2], nc_mat[5, -6:-2])
        self.assertEqual(dp_mat[5, 6:-2], nc_mat[5, 6:-2])
    
    def test_slice_int(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(10, 11, seed=0)
        self.assertEqual(dp_mat[0:10, 0], nc_mat[0:10, 0])
        self.assertEqual(dp_mat[0:10, 10], nc_mat[0:10, 10])
        self.assertEqual(dp_mat[4:8, 0], nc_mat[4:8, 0])
        self.assertEqual(dp_mat[2:7, 9], nc_mat[2:7, 9])
        self.assertEqual(dp_mat[1:8, 5], nc_mat[1:8, 5])
        self.assertEqual(dp_mat[0:1, 5], nc_mat[0:1, 5])
        self.assertEqual(dp_mat[9:10, 0], nc_mat[9:10, 0])
   
        self.assertEqual(dp_mat[0:-1, 0], nc_mat[0:-1, 0])
        self.assertEqual(dp_mat[-6:-2, 5], nc_mat[-6:-2, 5])
        self.assertEqual(dp_mat[6:-2, 5], nc_mat[6:-2, 5])

    def test_slice_slice(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(10, 11, seed=0)
        self.assertEqual(dp_mat[0:10, 0:11], nc_mat[0:10, 0:11])
        self.assertEqual(dp_mat[0:10, 10:11], nc_mat[0:10, 10:11])
        self.assertEqual(dp_mat[4:8, 0:1], nc_mat[4:8, 0:1])
        self.assertEqual(dp_mat[2:7, 9:10], nc_mat[2:7, 9:10])
        self.assertEqual(dp_mat[1:8, 5:8], nc_mat[1:8, 5:8])
        self.assertEqual(dp_mat[0:1, 5:10], nc_mat[0:1, 5:10])
        self.assertEqual(dp_mat[9:10, 0:1], nc_mat[9:10, 0:1])

        self.assertEqual(dp_mat[0:-1, 0:-1], nc_mat[0:-1, 0:-1])
        self.assertEqual(dp_mat[-8:-5, -7:-3], nc_mat[-8:-5, -7:-3])
        self.assertEqual(dp_mat[5:-1, 4:-2], nc_mat[5:-1, 4:-2])
        self.assertEqual(dp_mat[-2:-1, -4:-3], nc_mat[-2:-1, -4:-3])

class TestSubscriptSet(TestCase):
    def test_set_int_int(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(10, 11, seed=0)
        dp_mat[5, 10] = 101.3
        nc_mat[5, 10] = 101.3
        self.assertEqual(dp_mat[5, 10], nc_mat[5, 10])
        self.assertEqual(nc_mat[5, 10], 101.3)

        dp_mat[0, 10] = -101.9
        nc_mat[0, 10] = -101.9
        self.assertEqual(dp_mat[0, 10], nc_mat[0, 10])
        self.assertEqual(nc_mat[0, 10], -101.9)

    def test_set_int_slice(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(10, 11, seed=0)
        dp_mat[5, 10:11] = 101.4
        nc_mat[5, 10:11] = 101.4
        self.assertEqual(dp_mat[5, 10:11], nc_mat[5, 10:11])
        self.assertEqual(nc_mat[5, 10:11], 101.4)

        dp_mat, nc_mat = rand_dp_nc_matrix(11, 12, seed=0)
        dp_mat[0, 9:11] = [101, 102.2]
        nc_mat[0, 9:11] = [101, 102.2]
        self.assertEqual(dp_mat[0, 9:11], nc_mat[0, 9:11])
        dp_mat[10, 2:5] = [101, 102, 103.5]
        nc_mat[10, 2:5] = [101, 102, 103.5]
        self.assertEqual(dp_mat[10, 2:5], nc_mat[10, 2:5])
        self.assertEqual(dp_mat[0:11], nc_mat[0:11])

    def test_set_slice_int(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(11, 10, seed=0)
        dp_mat[10:11,5] = 101.8
        nc_mat[10:11, 5] = 101.8
        self.assertEqual(dp_mat[10:11, 5], nc_mat[10:11, 5])
        self.assertEqual(nc_mat[10:11, 5], 101.8)

        dp_mat, nc_mat = rand_dp_nc_matrix(11, 12, seed=0)
        dp_mat[9:11, 0] = [101, 102.3]
        nc_mat[9:11, 0] = [101, 102.3]
        self.assertEqual(dp_mat[9:11, 0], nc_mat[9:11, 0])
        dp_mat[2:5, 10] = [101, 102, 103.3]
        nc_mat[2:5, 10] = [101, 102, 103.3]
        self.assertEqual(dp_mat[2:5, 10], nc_mat[2:5, 10])
        self.assertEqual(dp_mat[0:11], nc_mat[0:11])

    def test_set_slice_slice(self):
        # TODO need to complete
        dp_mat, nc_mat = rand_dp_nc_matrix(11, 10, seed=0)
        dp_mat[10:11,5:6] = 101.4
        nc_mat[10:11, 5:6] = 101.4
        self.assertEqual(dp_mat[10:11, 5:6], nc_mat[10:11, 5:6])
        self.assertEqual(nc_mat[10:11, 5:6], 101.4)
        
        dp_mat, nc_mat = rand_dp_nc_matrix(11, 10, seed=0)
        dp_mat[10:11,5:7] = [101, 102.8]
        nc_mat[10:11, 5:7] = [101, 102.8]
        self.assertEqual(dp_mat[10:11, 5:7], nc_mat[10:11, 5:7])
        
        dp_mat, nc_mat = rand_dp_nc_matrix(11, 10, seed=0)
        dp_mat[9:11,5:6] = [101, 102.5]
        nc_mat[9:11, 5:6] = [101, 102.5]
        self.assertEqual(dp_mat[9:11, 5:6], nc_mat[9:11, 5:6])
        self.assertEqual(dp_mat[0:11], nc_mat[0:11])
        
        dp_mat, nc_mat = rand_dp_nc_matrix(11, 10, seed=0)
        dp_mat[8:11,5:7] = [[101, 102.9], [103, 104], [105, 106]]
        nc_mat[8:11, 5:7] = [[101, 102.9], [103, 104], [105, 106]]
        self.assertEqual(dp_mat[8:11, 5:7], nc_mat[8:11, 5:7])
        self.assertEqual(dp_mat[0:11], nc_mat[0:11])

        dp_mat, nc_mat = rand_dp_nc_matrix(11, 10, seed=0)
        dp_mat[0:3,8:10] = [[101.6, 102], [103, 104], [105, 106]]
        nc_mat[0:3, 8:10] = [[101.6, 102], [103, 104], [105, 106]]
        self.assertEqual(dp_mat[0:3, 8:10], nc_mat[0:3,8:10])
        self.assertEqual(dp_mat[0:11], nc_mat[0:11])
        
    def test_set_1D_int(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(11, 1, seed=0)
        dp_mat[0] = -101.5
        nc_mat[0] = -101.5
        self.assertEqual(dp_mat[0], nc_mat[0])
        self.assertEqual(nc_mat[0], -101.5)

        dp_mat, nc_mat = rand_dp_nc_matrix(1, 11, seed=0)
        dp_mat[0] = -101.4
        nc_mat[0] = -101.4
        self.assertEqual(dp_mat[0], nc_mat[0])
        self.assertEqual(nc_mat[0], -101.4)

    def test_set_1D_slice(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(11, 1, seed=0)
        dp_mat[0:1] = -101.1
        nc_mat[0:1] = -101.1
        self.assertEqual(dp_mat[0:1], nc_mat[0:1])
        self.assertEqual(nc_mat[0:1], -101.1)
        
        dp_mat[0:3] = [100.1, 101, 102]
        nc_mat[0:3] = [100.1, 101, 102]
        self.assertEqual(dp_mat[0:3], nc_mat[0:3])
        dp_mat[9:11] = [101, 102.4]
        nc_mat[9:11] = [101, 102.4]
        self.assertEqual(dp_mat[9:11], nc_mat[9:11])
        dp_mat[5:7] = [-101, 102.4]
        nc_mat[5:7] = [-101, 102.4]
        self.assertEqual(dp_mat[5:7], nc_mat[5:7])

        dp_mat, nc_mat = rand_dp_nc_matrix(1, 11, seed=0)
        dp_mat[0:1] = -101
        nc_mat[0:1] = -101
        self.assertEqual(dp_mat[0:1], nc_mat[0:1])
        self.assertEqual(nc_mat[0:1], -101)

        dp_mat[0:3] = [100, 101, 102.4]
        nc_mat[0:3] = [100, 101, 102.4]
        self.assertEqual(dp_mat[0:3], nc_mat[0:3])
        dp_mat[9:11] = [101, 102.4]
        nc_mat[9:11] = [101, 102.4]
        self.assertEqual(dp_mat[9:11], nc_mat[9:11])
        dp_mat[5:7] = [-101, 102.4]
        nc_mat[5:7] = [-101, 102.4]
        self.assertEqual(dp_mat[5:7], nc_mat[5:7])

    def test_set_2D_int(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(11, 4, seed=0)
        dp_mat[0] = [100, 101, 102, 103.4]
        nc_mat[0] = [100, 101, 102, 103.4]
        self.assertEqual(dp_mat[0], nc_mat[0])
        dp_mat[10] = [100, 101, 102, 103.8]
        nc_mat[10] = [100, 101, 102, 103.8]
        self.assertEqual(dp_mat[10], nc_mat[10])
        dp_mat[5] = [100, 101, 102.2, 103]
        nc_mat[5] = [100, 101, 102.2, 103]
        self.assertEqual(dp_mat[5], nc_mat[5])
                
    def test_set_2D_slice(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(11, 3, seed=0)
        dp_mat[0:3] = [[100, 101.2, 102], [103.8, 104, 105], [106, 107, 108]]
        nc_mat[0:3] = [[100, 101.2, 102], [103.8, 104, 105], [106, 107, 108]]
        self.assertEqual(dp_mat[0:3], nc_mat[0:3])
        dp_mat[9:11] = [[100, 101, 102.5], [103, 104, 105]]
        nc_mat[9:11] = [[100, 101, 102.5], [103, 104, 105]]
        self.assertEqual(dp_mat[9:11], nc_mat[9:11])
        dp_mat[5:7] = [[100, 101, 102], [103.6, 104, 105]]
        nc_mat[5:7] = [[100, 101, 102], [103.6, 104, 105]]
        self.assertEqual(dp_mat[5:7], nc_mat[5:7])
    
class TestSubscriptIntegration(TestCase):
    def test_slices(self):
        dp_mat, nc_mat = rand_dp_nc_matrix(11, 5, seed=0)
        dp_a = dp_mat[0:5, 0:2]
        nc_a = nc_mat[0:5, 0:2]
        dp_a[0] = [3, 4]
        nc_a[0] = [3, 4]
        self.assertEqual(dp_mat, nc_mat)
        self.assertEqual(dp_mat[0, 0:2], nc_mat[0, 0:2])
        dp_mat[4, 1] = 100;
        nc_mat[4, 1] = 100;
        self.assertEqual(dp_a[4, 1], nc_a[4, 1])

class TestMemLeak(TestCase):
    def test_multiple_operations(self):
        # expect this test to run slowly if there are memory leaks
        dp_mat1, nc_mat1 = rand_dp_nc_matrix(8000, 8000, seed=0)
        dp_mat2, nc_mat2 = rand_dp_nc_matrix(8000, 8000, seed=1)
        
        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

        is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)


