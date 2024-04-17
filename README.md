# numc

Here's what I did in project 4:

- Breakdown of work:   Theodora Worledge and I, Annamira O'Toole, were careful to breakdown the work between us as evenly as possible but ended up playing it mostly by ear.   We started working on the project two and a half weeks before the due date, hoping to get most of the work done for Task 1 and implementing the C to Python interface to the specification before Thanksgiving break.   We were able to work on it a lot over break and completed most of Task 1 and Task 2 and Task 3. However, we did notice that there was a memory leak from the matrix.c code that, according to Piazza, was unavoidable due to the initializion of the Python-C interface. So, we decided not to focus on that anymore and wanted to focus our energy on optimization. We are still wondering whether or not we could have fixed the memory leak or if it could possibly be slowing us down.

Annamira implemented most of easier matrix.c methods during the initial functionality runthrough, while Theodora implemented the fast and recursive version of the pow_matrix function so that it did repeated squaring.  Annamira did Task 2 (setting up the package).  For Task 3, Annamira implemented the interface between all matrix functions from matrix.c while Theodora took the task of dealing with the subscript and set_subscript functions, which were very tedious and difficult due to all of the complicated case breakdown. 

For Task 4, we both worked on almost everything so we won't include a careful breakdown of the changes between us, you can see the commit history on our Github repository history. We left this writeup of what we did to the very last day. 

- Whether or not Zoom calls helped, workflow, git, random stuff to take up words
- Cyberduck and spacing 

- Task 1

The important design decisions made here intitially for the matrix operations was to take a transpose of the second operand before sending it into 2 for loops instead of 3 for loops.   Additionally, we initially implemented recursive fast powers of matrices instead of the naive version for the sake of doing less work down the line.

We implemented deallocate_matrix recursively and did not enforce a correct ref_cnt for matrices that are not the initial matrix data allocated by allocate_matrix.   Furthermore, we kept an invariant where at any point, a matrix's path to the intitial origin parent should be kept at one.   We think of the data structure as a tree of pointers to root matrices, which is always of depth of only one.

- Task 2

Task 2 was pretty confusing to make sense of but overall like filling in blanks, and didn't require design decisions so we don't have anything to discuss here.

- Task 3

This did not require many design decisions either, besides the challenging tasks of subscript and set_subscript.   For those functions, fully implemented by Theodora, she did a careful case-by-case handling of possible queries to the matrix, then processing each.   We wondered whether or not there was a more efficient or less-tedious way of implementing the subscripting functions so that they took up fewer lines and were possibly recursive, but didn't end up trying to implement it that way once Theodora had a working implemention.  

- Task 4

We started out by reading through the Piazza thread and collecting ideas on potential changes for speedups.   The first strategy we addded was transposing the second operand of the matrix multiplication function so that it could be iterated through more easily.   This helps the efficiency of mul_matrix by allowing us to do a 2-layered for loop instead of a 3-layered naive implementation.   This is good because the second operand of mul_matrix must have its columns accessed, which requires a large number of separate indexing across each row, since each element of a column is stored in a different row.   By transposing ititially, we are then able to multiply a row vector times a row vector within the inner for loop of the iterations.   This also lends itself well to loading doubles into __m256d vectors later on during SIMD parallelization. 

The next optimization we added was initially allocating the space for the doubles as one large chunk in memory with 1D dimension mat->rows * mat->cols and this is to help the effect of caching of memory in the background.  By having the whole array pulled into cache memory the first time it is accessed, it saves time later on when the next rows are accessed after it.   To deal with the keeping the access format of the array the same as before, we altered the matrix struct so that it kept both the 1D memeory as "double_data", as well as storing an array caled data of dimension mat->rows, where each entry is a pointer to the double at the start of a row.  This way it can be subscribed as a 2D array the same as before. 

The next optimizations we added were to add SIMD Intel Intrinsics for parallelized operations of four 64 bit double precision numbers at a time in every operation within matrix.c.   This was able to be accomplished using the __m256d type from the Intel Intrinsics Manual.   We tried to reduce the total number of SIMD operations required by combining them when we could, such as:

In add and sub matrix, we used a repeated sequence of:
__m256d vector1 = _mm256_loadu_pd((__m256d*) (mat1->double_data + i));
__m256d vector2 = _mm256_loadu_pd((__m256d*) (mat2->double_data + i));
__m256d sum = _mm256_add_pd(vector1, vector2);

In mul matrix, we used a repeated sequence of:
vector1 = _mm256_loadu_pd((__m256d*) (m1 + i + 4));
vector2 = _mm256_loadu_pd((__m256d*) (m2 + i + 4));
sumb = _mm256_fmadd_pd (vector1, vector2, sumb);

By using _mm256_fmadd_pd instead of just multiplying and then adding seperately we are able to reduce at least one operation per loop.

In pow matrix we call mul_matrix and no direct SIMD operations were to be utilized.   However, we would ideally avoid reallocating space for the squared matrix struct each recursive level.   The workaround we tried to add is by writing a different version of the mul_matrix that only makes use of 1D arrays of doubles and passing just the array into each level instead of matrix structs.   However, we ran into issues implementing this near the end of the project and ultimately abandoned it.  

In negate matrix we use repeated sequence of:
__m256d vector = _mm256_loadu_pd((__m256d*) (mat->double_data + i));
__m256d neg_vector = _mm256_sub_pd(zeros, vector);
_mm256_storeu_pd(result->double_data + i, neg_vector);

In absolute value we reduce total SIMD operations by defining a static variable 
static double double_mask = nan("0xfffffffffffff");
__m256d mask = _mm256_set1_pd(double_mask); 
Which we use to mask the sign bit by calling a sequence of 
__m256d vector = _mm256_loadu_pd((__m256d*) (mat->double_data + i));
__m256d abs_vector = _mm256_and_pd(mask, vector);
_mm256_storeu_pd(result->double_data + i, abs_vector);


After adding SIMD parallelization to each function as much as possible, we added simple openMP multithreading by adding the macro # pragma omp parallel for the outer for loop on any iteration-intensive for loop that is not at risk of suffering from a data race with over written data settings.

After adding openOMP multithreading, we unrolled the for loops as much as possible within every method.   Instead of doing groups of four operations at a time using SIMD vectors, handling one tail case at the end for multiplying the remaining elements modulo 4.   We went futher to unroll to doing four groups of the 4-double-precsion operations at once within each body of the for loop.   Thus, we are doing 16 operations per for loop, which results in two tail cases: the remaining groups of 4 up to the total dimension floor divded by four, and then the remaining individual double multiplications for the 3, 2, or 1 doubles remaining in the very last tail case.

- Task 5

We decided to leave this task to the last step in our Project 4: Numc Optimization since we were rethinking how to improve the power function even more and more on the last day of submissions but ultimately kept our Friday submission.

